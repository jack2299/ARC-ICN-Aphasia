# diagnostic.py
#
# Prints basic properties of a NIfTI mask or image:
# file size, shape, data type, unique values, and number of non-zero voxels.
#
# Useful for sanity-checking ICN or BM masks before downstream analysis.

import os
import numpy as np
import nibabel as nib

# =========================================================================
# EDIT ONLY THESE
# =========================================================================
FOLDER = "path/to/mask/folder"
FILES = [
    "mask1.nii",
    "mask2.nii.gz",
]
# =========================================================================

for f in FILES:
    path = os.path.join(FOLDER, f)

    if not os.path.exists(path):
        print(f"\n--- {f} ---")
        print("File not found:", path)
        continue

    print(f"\n--- {f} ---")
    print("Size (bytes):", os.path.getsize(path))

    img = nib.load(path, mmap=False)
    print("Shape:", img.shape)
    print("Data type:", img.get_data_dtype())

    data = np.asanyarray(img.dataobj)
    print("Unique values:", np.unique(data))
    print("Sum of voxels > 0:", int((data > 0).sum()))
