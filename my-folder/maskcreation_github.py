# maskcreation.py
#
# Generates the 18 BrainMap20 ICN masks from the ICN_Atlas .atl file.
# Each mask is thresholded at > 3 and saved as a NIfTI file in the output directory.
#
# Requires:
#   - ICN_Atlas toolbox (r20180306 or compatible), which provides
#     atlas_BrainMap20.atl
#
# Output:
#   - ICN01_mask.nii.gz ... ICN18_mask.nii.gz

import os
import numpy as np
import scipy.io as sio
import nibabel as nib

# =========================================================================
# EDIT ONLY THESE PATHS
# =========================================================================
ATLAS_FILE = "path/to/atlas_BrainMap20.atl"
OUTPUT_DIR = "path/to/output/ICN_masks"
# =========================================================================

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Load atlas
mat = sio.loadmat(ATLAS_FILE)

# Variable names as stored in the .atl file
icn_data = mat["atl_mapdatamatrix"]  # shape (91, 109, 91, 20)
affine = mat["atl_Ta"]               # 4x4 transformation matrix

# Generate masks for the first 18 ICNs (exclude the two artefact components)
for n in range(18):
    # Extract ICN n+1 (0-indexed here, but ICN1..ICN18 in the atlas)
    vol = icn_data[:, :, :, n]       # shape (91, 109, 91)
    mask = (vol > 3).astype(np.float32)

    # Create NIfTI image
    nii = nib.Nifti1Image(mask, affine)
    fname = os.path.join(OUTPUT_DIR, f"ICN{n+1:02d}_mask.nii.gz")
    nib.save(nii, fname)
    print(f"Saved ICN{n+1:02d} ({int(mask.sum())} voxels)")

print(f"All 18 masks saved in {OUTPUT_DIR}")