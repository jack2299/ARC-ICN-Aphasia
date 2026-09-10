# tobit_sensitivity.py
#
# Censored Tobit sensitivity analysis for the BM17 engagement-severity
# relationship, using the clean aphasia sample.
#
# Rationale: WAB-R is bounded at 100, which can cause mild ceiling effects in
# ordinary least squares. Tobit regression with upper censoring at 100 tests
# whether the OLS estimates are biased by that ceiling.
#
# Requires:
#   - ARC_03b_v3_Master_Wide.mat  (master wide table from the pipeline)
#   - Python packages: numpy, pandas, scipy, statsmodels

import numpy as np
import pandas as pd
from scipy.io import loadmat
from statsmodels.censored.tobit_model import Tobit
import statsmodels.api as sm

# =========================================================================
# EDIT ONLY THESE
# =========================================================================
MASTER_FILE = 'ARC_03b_v3_Master_Wide.mat'
OUTPUT_CSV  = 'Tobit_Sensitivity_Results.csv'
# =========================================================================

# =========================================================================
# 1. Load data and build clean aphasia sample
# =========================================================================
mat = loadmat(MASTER_FILE)
mw = mat['masterWide']

if mw.dtype.names:
    mw = mw[0, 0]

def struct_to_df(s):
    return pd.DataFrame({name: s[name].flatten() for name in s.dtype.names})

masterWide = struct_to_df(mw)

# Compute SexNumeric if not present
if 'SexNumeric' not in masterWide.columns:
    if 'Sex' in masterWide.columns:
        masterWide['SexNumeric'] = (masterWide['Sex'].astype(str).str.upper() == 'M').astype(int)
    else:
        raise ValueError("Neither SexNumeric nor Sex found in master wide file.")

restOnly = ['M2088','M2097','M2100','M2101','M2113','M2114','M2117','M2118',
            'M2122','M2126','M2129','M2131','M2135','M2140','M2141','M2142',
            'M2143','M2144','M2145','M2146','M2149','M2150','M2151','M2152',
            'M2153','M2155','M2156','M2158','M2159','M2160','M2162','M2164',
            'M2165','M2169','M2184','M2254']

masterWide['isAph']  = masterWide['GroupRole'].astype(str).str.contains('Aphasia')
masterWide['hasC2']  = ~np.isnan(masterWide['C2_ICN17_IRi'])
masterWide['isRest'] = masterWide['PatientID'].isin(restOnly)

clean = masterWide[masterWide['isAph'] & masterWide['hasC2'] & ~masterWide['isRest']]
print(f'Clean sample N = {len(clean)}')

# =========================================================================
# 2. Level 3: WAB ~ BM17 C2 IRi + Age + Sex + Days
# =========================================================================
X3 = clean[['C2_ICN17_IRi', 'Age_At_Stroke', 'SexNumeric', 'Days_Post_Stroke']].copy()
X3 = X3.dropna()
y3 = clean.loc[X3.index, 'WAB_AQ'].astype(float)
X3 = sm.add_constant(X3)

ols3 = sm.OLS(y3, X3).fit()
print('\nLevel 3 OLS summary:')
print(ols3.summary().tables[1])

tob3 = Tobit(y3, X3, cens='upper', right=100, method='bfgs')
res3 = tob3.fit(disp=0)
print('\nLevel 3 Tobit summary:')
print(res3.summary().tables[1])

# =========================================================================
# 3. Level 4: WAB ~ C2 + C1 + Age + Sex + Days
# =========================================================================
X4 = clean[['C2_ICN17_IRi', 'C1_ICN17_IRi', 'Age_At_Stroke', 'SexNumeric', 'Days_Post_Stroke']].copy()
X4 = X4.dropna()
y4 = clean.loc[X4.index, 'WAB_AQ'].astype(float)
X4 = sm.add_constant(X4)

ols4 = sm.OLS(y4, X4).fit()
print('\nLevel 4 OLS summary:')
print(ols4.summary().tables[1])

tob4 = Tobit(y4, X4, cens='upper', right=100, method='bfgs')
res4 = tob4.fit(disp=0)
print('\nLevel 4 Tobit summary:')
print(res4.summary().tables[1])

# =========================================================================
# 4. Comparison table
# =========================================================================
def compare(ols, tob):
    return (
        ols.params['C2_ICN17_IRi'], ols.pvalues['C2_ICN17_IRi'],
        tob.params['C2_ICN17_IRi'], tob.pvalues['C2_ICN17_IRi'],
    )

l3 = compare(ols3, res3)
l4 = compare(ols4, res4)

print('\n==== Comparison ====')
print('Level 3 - BM17 C2 IRi coefficient:')
print(f'  OLS:   beta = {l3[0]:.4f}, p = {l3[1]:.4f}')
print(f'  Tobit: beta = {l3[2]:.4f}, p = {l3[3]:.4f}')
print('Level 4 - BM17 C2 IRi coefficient (from C2 + C1 model):')
print(f'  OLS:   beta = {l4[0]:.4f}, p = {l4[1]:.4f}')
print(f'  Tobit: beta = {l4[2]:.4f}, p = {l4[3]:.4f}')

# Save comparison table
results = pd.DataFrame([
    {'Level': 'Level 3', 'Method': 'OLS',   'Beta': l3[0], 'p': l3[1]},
    {'Level': 'Level 3', 'Method': 'Tobit', 'Beta': l3[2], 'p': l3[3]},
    {'Level': 'Level 4', 'Method': 'OLS',   'Beta': l4[0], 'p': l4[1]},
    {'Level': 'Level 4', 'Method': 'Tobit', 'Beta': l4[2], 'p': l4[3]},
])
results.to_csv(OUTPUT_CSV, index=False)
print(f'\nSaved: {OUTPUT_CSV}')
