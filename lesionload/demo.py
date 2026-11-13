"""Generate a small atlas and DWI fixture for all three report modes."""

import numpy as np
import nibabel as nib

from anatomy_review import file_record, save_new
from spatial_contract import REVIEW_SCHEMA
from structural_multimodal import sitk_module
from .spatial import FA_DIRECTIONS, T2_DIRECTIONS
from .workflow import SCHEMA, run


def demo(output):
    root = output.resolve()
    inputs = root / "inputs"
    inputs.mkdir(parents=True)
    affine = np.eye(4)
    grid = (6, 5, 5)
    source = (5, 5, 5)

    def image(name, values):
        path = inputs / f"{name}.nii.gz"
        nib.save(nib.Nifti1Image(values, affine), path)
        return file_record(path)

    lesion = np.zeros(source, dtype=np.uint8)
    lesion[0, :1, 1] = 1
    lesion[:3, 1, 2] = 1
    left = np.zeros(grid, dtype=np.uint8)
    left[:, :1, 1] = 1
    left[5, 0, 1] = 0
    left[:4, 1:3, 2] = 1
    left[5, 1:3, 2] = 1
    right = np.zeros(grid, dtype=np.uint8)
    right[2:5, 3, 1:4] = 1
    brain = image("reference", np.indices(grid).sum(axis=0).astype(np.float32) + 2)
    native = image("native_T2", np.ones(source, dtype=np.float32))
    original_lesion = image("original_lesion", lesion)
    fit = np.ones(source, dtype=np.uint8)
    fit[0, 0, 1] = 0
    fit_mask = image("fit_mask", fit)
    review = inputs / "synthetic_review.json"
    save_new(review, {"synthetic": True})
    sitk = sitk_module()
    transform = inputs / "identity.h5"
    sitk.WriteTransform(sitk.Transform(3, sitk.sitkIdentity), str(transform))
    matrix = inputs / "identity.mat"
    np.savetxt(matrix, np.eye(4))
    empty = image("empty", np.zeros(grid, dtype=np.uint8))
    banks = {}
    for method, origin in (("SMATT", "M1"), ("HCPA", "CST"), ("Boyne2025", "Whole")):
        banks[method] = {"reference": brain,
                         "tracts": {f"Left_{origin}": image(f"{method}_left", left),
                                    f"Right_{origin}": image(f"{method}_right", right),
                                    "Left_Empty": empty}}
    visit = {"dataset": "generated", "subject": "demo", "session": "01"}
    config = {"schema": SCHEMA, "synthetic": True, "visit": visit, "methods": list(banks),
              "synthetic_atlases": banks, "synthetic_gen6_to_2009c_pull": file_record(transform),
              "t2": {"source": native, "lesion": original_lesion, "t1": native, "mni": brain,
                     "t1_to_t2_pull": file_record(transform), "mni_to_t1_pull": file_record(transform),
                     "review": file_record(review)}, "affected_side": "Left",
              "fa": {"b0": native, "fa": image("FA", np.full(source, .625, dtype=np.float32)),
                     "md": image("MD", np.full(source, .001, dtype=np.float32)),
                     "fit_mask": fit_mask, "reviewed_mask_in_t1": fit_mask,
                     "flirt_dwi_to_t1": file_record(matrix), "flirt_inverse": file_record(matrix),
                     "world_forward_ras": np.eye(4).tolist(), "review": file_record(review),
                     "model_receipt": file_record(review),
                     "t1_sha256": native["sha256"], "valid_native_voxels": int(fit.sum())}}
    for key, scope, roles, directions in (
        ("spatial_review", "lesion_load", config["t2"], T2_DIRECTIONS),
        ("fa_spatial_review", "lesion_load_fa", {**config["fa"], "t1": native}, FA_DIRECTIONS),
    ):
        path = inputs / f"{scope}_review.json"
        save_new(path, {"schema": REVIEW_SCHEMA, "decision": "APPROVED", "visit": visit,
                        "scope": scope, "branch": None, "point_directions": directions,
                        "reviewer": "Synthetic fixture generator", "reviewed_utc": "2000-01-01T00:00:00Z",
                        "basis": "Synthetic inputs with identity transforms", "evidence": {"fixture": file_record(review)},
                        "input_sha256": {role: entry["sha256"] for role, entry in roles.items()
                                         if isinstance(entry, dict) and "sha256" in entry}})
        config[key] = file_record(path)
    config_path = root / "synthetic.local.json"
    save_new(config_path, config)
    run(config_path, root / "results")
    return {"config": str(config_path), "receipt": str(root / "results" / "receipt.json")}
