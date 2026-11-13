# Spatial input review

Use source dataset, subject and session IDs in the config's `visit`. For a local visit, use a stable acquisition ID for `session`. The report folder name can differ from the source session. Keep case-specific receipts outside the checkout.

A review JSON uses `schema: "neo-spatial-input-review-1"`, `decision: "APPROVED"`, the same `visit`, and the `scope` listed below. Include the reviewer, decision time (`reviewed_utc`), a statement of what the reviewer approved (`basis`), and one or more named file records in `evidence`. File records have `path` and `sha256` fields. Evidence paths are relative to the review, while config paths are relative to the config. Include the exact files or structured source decisions that support the review. The workflow verifies the receipt, its evidence files and each input role. Reviewers are responsible for confirming what the original decisions cover.

## CST atlas mapping and laterality

1. Prepare the [map-label config](../configs/examples/tractography_label.example.json) with the b0, T1, 2009c template and atlas reference. The binary label and atlas reference must describe the same physical lattice. Storage-axis flips and permutations are allowed.
2. Use `python run_tractography.py audit-label --config map.local.json --output <new-mask.nii.gz>` to make a technical candidate and print its input SHA-256 map and LPS point directions. Inspect the native T1 and b0 in matched planes. On each side of a transform, mark corresponding anatomical landmarks independently of that transform. Include at least three distinct, non-collinear points per segment. The landmark file contains `visit`, `grids_sha256` (the four reference roles and their SHA-256 values), and `segments` (transform roles to arrays of `name`, `from_lps_mm`, `to_lps_mm`). Coordinates are physical LPS millimetres. Add the landmark record to the config and rerun the audit to get per-point error and its summary.
3. Following review of the errors, coverage and anatomy, create the approval for `scope: "cst_atlas_mapping"` and `branch: null`. Copy the audit's `input_sha256` and `point_directions_lps` into `input_sha256` and `point_directions` in the review, respectively. Record the maximum allowable landmark distance for each segment under `landmark_max_mm`, in physical millimetres. Choose limits for the anatomy, acquisition and intended use. Add the review record to `spatial_review` and use `map-label` with a new output path. The output is sampled once on the b0 grid. An audit label is a technical candidate pending anatomical review.
4. For real `screen` jobs the same coordinate chain defines the side of each streamline origin. Use the same reference, anatomy and landmark records and an approval scoped to `cst_laterality`. The reviewed hashes cover the reference roles, transforms and landmarks. Omit the single atlas `label` hash for this scope. The `dwi_to_t1_ras` hash uses the 4×4 matrix as row-major little-endian float64 bytes.

For Gen6 labels the sampling chain is DWI → T1 → 2009c → Gen6 in LPS. The pinned TemplateFlow H5 is named for moving Gen6 images to 2009c. Its point direction is 2009c → Gen6. The inverse used for lesion-load sampling points in the opposite direction. Template-to-template anatomical point error and inverse consistency answer different questions. An approval can be issued only after examining the actual case and application.

## Lesion load

Use `spatial_review` with `scope: "lesion_load"`. Bind `source`, `lesion`, `t1`, `mni`, `t1_to_t2_pull` and `mni_to_t1_pull` by SHA-256. For the frozen T2 path, also bind `structural_request`, `structural_receipt` and the selected `branch`. Review the branch and the exact native lesion against the approved source before issuing the receipt. For an archived explicit chain, also bind the old `review` file as provenance and set `branch: null`. Confirm what the original approval covers before reusing it. Use these values in `point_directions`.

```json
{
  "t1_to_t2_pull": ["LPS_mm", "T1w", "native_T2w"],
  "mni_to_t1_pull": ["LPS_mm", "MNI152NLin2009cAsym", "T1w"]
}
```

If using optional FA, supply `fa_spatial_review` with `scope: "lesion_load_fa"` and `branch: null`. Bind `t1`, `b0`, `fa`, `md`, `fit_mask`, `reviewed_mask_in_t1`, `flirt_dwi_to_t1`, `flirt_inverse`, `review` (the registration review), and `model_receipt` by SHA-256. Confirm the subject and acquisition, selected FLIRT method, and model and output identity in the original records. Set the directions to the following.

```json
{
  "flirt_dwi_to_t1": ["FSL_scaled_voxels", "DWI", "T1w"],
  "flirt_inverse": ["FSL_scaled_voxels", "T1w", "DWI"]
}
```

Move an input by updating its config path and retaining the same hash. If evidence files move, rebind their paths in a new receipt while retaining their content hashes. Existing reports and result directories are readable in place. Synthetic demos generate fixture-specific evidence automatically.
