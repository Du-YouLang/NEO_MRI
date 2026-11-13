"""Physical-space atlas mapping and validated native DWI grids."""
import hashlib
import json

import nibabel as nib
import numpy as np

from spatial_contract import MNI2009C, MNI6, checked, digest, locked_template_bridge, require_review, visit_identity


def load_reference(path):
    ref = nib.as_closest_canonical(nib.load(path))
    affine = ref.affine
    if len(ref.shape) != 3 or not np.isfinite(affine).all():
        raise ValueError("A finite three-dimensional corrected DWI reference is required")
    axes = affine[:3, :3] / nib.affines.voxel_sizes(affine)
    if not np.allclose(axes.T @ axes, np.eye(3), atol=1e-5, rtol=0):
        raise ValueError("The DWI grid must be orthogonal; oblique rotation is supported")
    return ref


def load_mask(path, ref):
    image = nib.as_closest_canonical(nib.load(path))
    if image.shape != ref.shape or not np.allclose(image.affine, ref.affine, atol=1e-4, rtol=0):
        raise ValueError(f"Mask does not share the DWI physical grid: {path}")
    values = np.asanyarray(image.dataobj)
    if not np.isin(values, (0, 1)).all():
        raise ValueError(f"Expected a binary hard label: {path}")
    return values.astype(bool)


def atlas_stages(anatomy, base, source_space):
    """Read consecutive LPS point maps, starting at output b0."""
    import SimpleITK as sitk

    ras = np.asarray(anatomy["dwi_to_t1_ras"], dtype=float)
    if (ras.shape != (4, 4) or not np.isfinite(ras).all()
            or not np.allclose(ras[3], [0, 0, 0, 1], atol=1e-7, rtol=0)
            or not np.allclose(ras[:3, :3].T @ ras[:3, :3], np.eye(3), atol=1e-4, rtol=0)
            or not np.isclose(np.linalg.det(ras[:3, :3]), 1, atol=1e-4, rtol=0)):
        raise ValueError("DWI-to-T1 must be a rigid physical RAS matrix")
    flip = np.diag([-1, -1, 1, 1])
    lps = flip @ ras @ flip
    pose = sitk.AffineTransform(lps[:3, :3].ravel().tolist(), lps[:3, 3].tolist())
    t1_to_standard = sitk.ReadTransform(str(checked(anatomy["mni2009c_to_t1_image_h5"], base)))
    stages = {"dwi_to_t1_ras": (["DWI", "T1w"], pose),
              "mni2009c_to_t1_image_h5": (["T1w", MNI2009C], t1_to_standard)}
    if source_space == MNI6:
        stages["mni6_to_mni2009c_image_h5"] = (
            [MNI2009C, MNI6], locked_template_bridge(anatomy["mni6_to_mni2009c_image_h5"], base))
    elif source_space != MNI2009C:
        raise ValueError(f"Unsupported atlas space: {source_space}")
    return stages


def landmark_errors(entries, stages, grids):
    """Compare transformed points to independently marked counterparts, in millimetres."""
    if set(entries) != set(stages):
        raise ValueError("Provide landmarks for every segment of the mapping")
    errors = {}
    for role, (spaces, transform) in stages.items():
        rows = entries[role]
        start = np.asarray([row["from_lps_mm"] for row in rows], dtype=float)
        end = np.asarray([row["to_lps_mm"] for row in rows], dtype=float)
        if (len(rows) < 3 or start.shape != (len(rows), 3) or end.shape != start.shape
                or not np.isfinite(start).all() or not np.isfinite(end).all()
                or np.linalg.matrix_rank(start[1:] - start[0]) < 2
                or len({row["name"] for row in rows}) != len(rows)):
            raise ValueError(f"Expected distinct, non-collinear paired landmarks for {role}")
        for space, points in zip(spaces, (start, end), strict=True):
            grid = grids[space]
            indices = np.asarray([grid.TransformPhysicalPointToContinuousIndex(tuple(p)) for p in points])
            if np.any((indices < -.5) | (indices >= np.asarray(grid.GetSize()) - .5)):
                raise ValueError(f"Landmarks outside the {space} field of view for {role}")
        actual = np.asarray([transform.TransformPoint(tuple(p)) for p in start])
        distances = np.linalg.norm(actual - end, axis=1)
        if not np.isfinite(distances).all():
            raise ValueError(f"Nonfinite mapped landmarks for {role}")
        target = grids[spaces[1]]
        actual_indices = np.asarray([target.TransformPhysicalPointToContinuousIndex(tuple(p)) for p in actual])
        if np.any((actual_indices < -.5) | (actual_indices >= np.asarray(target.GetSize()) - .5)):
            raise ValueError(f"Mapped landmarks fall outside the {spaces[1]} field of view for {role}")
        errors[role] = {"point_direction": spaces, "p95_mm": float(np.percentile(distances, 95)),
                        "max_mm": float(distances.max()),
                        "points": [dict(name=row["name"], distance_mm=float(distance))
                                   for row, distance in zip(rows, distances, strict=True)]}
    return errors


def same_physical_lattice(left, right):
    """Compare the same lattice even if stored with different axis order or signs."""
    a, b = (nib.as_closest_canonical(nib.load(path)) for path in (left, right))
    return a.shape == b.shape and np.allclose(a.affine, b.affine, atol=1e-4, rtol=0)


def compose_stages(stages):
    from structural_multimodal import compose_pull

    _, pull = next(iter(stages.values()))
    for _, transform in list(stages.values())[1:]:
        pull = compose_pull(transform, pull)
    return pull


def mapping_context(config, base, source_space, *, technical_replay=False, label=None):
    """Validate the pull and the evidence used to choose it."""
    import SimpleITK as sitk
    from t2_lesion_propagation import check_nifti_grid

    visit = visit_identity(config.get("visit"))
    anatomy = config["anatomy"]
    if source_space not in (MNI2009C, MNI6):
        raise ValueError(f"Unsupported atlas space: {source_space}")
    roles = ("reference", "t1_reference", "mni2009c_reference", "atlas_reference")
    paths = {role: checked(config[role], base) for role in roles}
    for path in paths.values():
        check_nifti_grid(path)
    grids = {space: sitk.ReadImage(str(paths[role])) for space, role in (
        ("DWI", "reference"), ("T1w", "t1_reference"), (MNI2009C, "mni2009c_reference"))}
    atlas = sitk.ReadImage(str(paths["atlas_reference"]))
    if source_space == MNI6:
        grids[MNI6] = atlas

    stages = atlas_stages(anatomy, base, source_space)
    directions = {role: spaces for role, (spaces, _) in stages.items()}
    ras_bytes = np.asarray(anatomy["dwi_to_t1_ras"], dtype="<f8").tobytes()
    bindings = {role: config[role]["sha256"].lower() for role in roles}
    if label is not None:
        bindings["label"] = label["sha256"].lower()
    bindings["dwi_to_t1_ras"] = hashlib.sha256(ras_bytes).hexdigest()
    bindings.update({role: anatomy[role]["sha256"].lower() for role in stages if role != "dwi_to_t1_ras"})
    review, errors = None, None
    if config.get("landmarks") is not None:
        points = checked(config["landmarks"], base)
        data = json.loads(points.read_text(encoding="utf-8-sig"))
        if data["visit"] != visit or data["grids_sha256"] != {role: bindings[role] for role in roles}:
            raise ValueError("Landmarks belong to a different visit or physical reference")
        errors = landmark_errors(data["segments"], stages, grids)
        bindings["landmarks"] = config["landmarks"]["sha256"].lower()
    if not technical_replay:
        if errors is None:
            raise ValueError("Independently marked anatomical landmarks are required for a reviewed map")
        if "spatial_review" not in config:
            raise ValueError("A spatial review of the landmarks and exact inputs is required")
        review = require_review(config["spatial_review"], base, visit=visit,
                                scope="cst_atlas_mapping" if label is not None else "cst_laterality",
                                inputs=bindings, directions=directions, landmark_errors=errors)
    return compose_stages(stages), grids["DWI"], paths["atlas_reference"], {
        "visit": visit, "input_sha256": bindings, "point_directions_lps": directions,
        "landmark_errors": errors, "input_review": review,
    }


def mapped_label(config, base, output, *, technical_replay=False):
    """Apply a composed pull to a binary atlas label in one interpolation."""
    import SimpleITK as sitk
    from structural_multimodal import resample
    from t2_lesion_propagation import check_nifti_grid, write_verified

    source_path = checked(config["label"], base)
    check_nifti_grid(source_path)
    source = sitk.ReadImage(str(source_path))
    if not np.isin(sitk.GetArrayViewFromImage(source), (0, 1)).all():
        raise ValueError("Threshold the source atlas to an explicit binary label first")
    if output.exists():
        raise FileExistsError(output)
    source_space = config["source_space"]
    pull, reference, atlas_path, evidence = mapping_context(
        config, base, source_space, technical_replay=technical_replay, label=config["label"])
    if not same_physical_lattice(source_path, atlas_path):
        raise ValueError("Atlas label and atlas reference do not share the same physical grid")
    result = resample(source, reference, pull, label=True)
    write_verified(sitk, result, output)
    return {"output": str(output), "sha256": digest(output),
            "voxels": int(sitk.GetArrayViewFromImage(result).sum()),
            "source_space": source_space, "interpolations": 1, **evidence,
            "status": "TECHNICAL_REPLAY_ONLY" if technical_replay else "REVIEWED_INPUTS_MAPPED",
            }
