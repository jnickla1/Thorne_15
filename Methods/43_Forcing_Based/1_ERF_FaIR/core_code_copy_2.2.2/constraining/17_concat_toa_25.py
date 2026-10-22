#!/usr/bin/env python3
import numpy as np
import glob
import re
import os

# directory containing the per-member .npy files
indir = "../../../../../output/all-headstails25/toa_imbalance"
pattern = os.path.join(indir, "all-headstails25r*cut2025_toa_imbalance_member*.npy")

files = sorted(
    glob.glob(pattern),
    key=lambda f: int(re.search(r"member(\d+)\.npy", f).group(1))
)

print(f"Found {len(files)} files")

arrays = []
for f in files:
    arr = np.load(f)
    if arr.shape != (2, 177):
        raise ValueError(f"Unexpected shape {arr.shape} in {f}")
    arrays.append(arr)

big_array = np.stack(arrays, axis=0)  # (100, 2, 175)

outfile = os.path.join(indir, "all-headstails_current_TOA_ensemble_25.npy")
np.save(outfile, big_array)
print(f"Saved {outfile} with shape {big_array.shape}")

