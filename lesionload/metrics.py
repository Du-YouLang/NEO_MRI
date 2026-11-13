"""Coverage-aware volume and axial measurements on a common physical grid."""

from fractions import Fraction

import numpy as np


def counts(n, c, h):
    n, c, h = map(int, (n, c, h))
    u = n - c
    return {
        "atlas_voxels": n, "covered_voxels": c, "hits_voxels": h, "unknown_voxels": u,
        "hits_ml": h / 1000,
        "fraction_whole": h / n if n and u == 0 else None,
        "fraction_covered": h / c if c else None,
        "fraction_bounds": [h / n, (h + u) / n] if n else None,
    }


def measure(tract, lesion, covered, z_mm, *, hcpa=False):
    """Return one summary and rows at physical axial levels, on a 1 mm grid."""
    n = tract.sum(axis=(0, 1))
    c = (tract & covered).sum(axis=(0, 1))
    h = (tract & lesion & covered).sum(axis=(0, 1))
    whole = counts(n.sum(), c.sum(), h.sum())
    whole.update(peak_fraction=None, peak_fraction_bounds=None, peak_rows=[])
    area_max = int(n.max())
    if hcpa:
        whole.update(max_slice_voxels=area_max, max_wll_voxel=None, max_wll_bounds_voxel=None)
    if not area_max:
        return whole, []
    levels = np.flatnonzero(n)
    rows = []
    lower, upper = [], []
    for k in levels:
        row = counts(n[k], c[k], h[k])
        row["z_ras_mm"] = float(z_mm[k])
        lo = Fraction(int(h[k]), int(n[k]))
        hi = Fraction(int(h[k] + n[k] - c[k]), int(n[k]))
        lower.append(lo)
        upper.append(hi)
        if hcpa:
            row["wll_bounds_voxel"] = [float(lo * area_max), float(hi * area_max)]
            row["wll_voxel"] = float(lo * area_max) if lo == hi else None
        rows.append(row)
    lo, hi = max(lower), max(upper)
    whole["peak_fraction_bounds"] = [float(lo), float(hi)]
    whole["peak_fraction"] = float(lo) if lo == hi else None
    if lo > 0:
        whole["peak_rows"] = [row for row, value in zip(rows, lower, strict=True) if value == lo]
    if hcpa:
        whole["max_wll_bounds_voxel"] = [float(lo * area_max), float(hi * area_max)]
        whole["max_wll_voxel"] = float(lo * area_max) if lo == hi else None
    return whole, rows


def fa_summary(roi, sampled, acquired):
    values = sampled[roi & np.isfinite(sampled)]
    n = len(values)
    return {"roi_voxels": int(roi.sum()), "dwi_interpolable_fov_voxels": int((roi & acquired).sum()),
            "valid_fa_voxels": n, "FA_mean": float(values.mean()) if n else None,
            "FA_median": float(np.median(values)) if n else None}


def fa_groups(tract, lesion, covered, sampled, acquired, z_mm, *, mode):
    known = tract & covered
    if mode == "bilateral":
        groups = {"covered_all": known, "inside": known & lesion, "outside": known & ~lesion}
    elif mode == "contralateral":
        groups = {"contralateral_homolog": known & ~lesion}
    else:
        groups = {"inside": known & lesion, "outside": known & ~lesion}
    whole, layers = [], []
    for group, roi in groups.items():
        whole.append(dict(group=group, **fa_summary(roi, sampled, acquired)))
        for k in np.flatnonzero(tract.any(axis=(0, 1))):
            layers.append(dict(group=group, z_ras_mm=float(z_mm[k]),
                               **fa_summary(roi[:, :, k], sampled[:, :, k], acquired[:, :, k])))
    return whole, layers
