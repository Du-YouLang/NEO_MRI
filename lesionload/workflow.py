"""Run independent tract overlap measurements from one configured subject."""

import csv
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from anatomy_review import file_record, save_new
from structural_multimodal import compose_pull, resample, sitk_module
from t2_lesion_propagation import write_verified
from .assets import BRIDGE_SOURCE, METHODS, locked_file, real_atlas
from .bridge import invert_template_bridge
from .metrics import fa_groups, measure
from .report import render
from .spatial import (MNI6, binary, checked, image_grid, itk_grid, map_pair, native_fa, sample_fa,
                      t2_chain, target_grid)


SCHEMA = "neo-lesionload-1"


def synthetic_atlas(spec, base):
    path = checked(spec["reference"], base)
    grid = target_grid(path)
    bank = {name: binary(image_grid(checked(entry, base), grid)) for name, entry in spec["tracts"].items()}
    if not bank or "Union" in bank:
        raise ValueError("Provide synthetic tracts with distinct names; Union is derived")
    bank["Union"] = np.logical_or.reduce(list(bank.values()))
    return path, grid, bank, [file_record(path)]


def write_csv(path, rows):
    if rows:
        with path.open("x", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


def write_mask(path, values, target):
    sitk = sitk_module()
    image = sitk.GetImageFromArray(values.T.astype(np.uint8))
    image.CopyInformation(target)
    write_verified(sitk, image, path)
    return file_record(path)


def calculate_bank(bank, lesion, covered, z, method, fa, affected_side):
    rows, layers, diffusion, diffusion_layers = [], [], [], []
    for name, tract in bank.items():
        whole, axial = measure(tract, lesion, covered, z, hcpa=method == "HCPA")
        rows.append(dict(region=name, **whole))
        layers.extend(dict(region=name, **row) for row in axial)
        if fa is not None and (method == "Boyne2025" or name.startswith(("Left_", "Right_"))):
            mode = "bilateral" if method == "Boyne2025" else (
                "affected" if name.startswith(f"{affected_side}_") else "contralateral")
            fa_whole, fa_axial = fa_groups(tract, lesion, covered, *fa, z, mode=mode)
            diffusion.extend(dict(region=name, **row) for row in fa_whole)
            diffusion_layers.extend(dict(region=name, **row) for row in fa_axial)
    return {"whole": rows, "axial": layers, "fa_whole": diffusion, "fa_axial": diffusion_layers}


def run(config_path, output):
    config_path = config_path.resolve(strict=True)
    config = json.loads(config_path.read_text(encoding="utf-8-sig"))
    base = config_path.parent
    output = output.resolve()
    if config.get("schema") != SCHEMA:
        raise ValueError(f"Expected schema {SCHEMA}")
    methods = config["methods"]
    if not methods or len(set(methods)) != len(methods) or set(methods) - METHODS.keys():
        raise ValueError("Select one or more distinct supported methods")
    synthetic = config.get("synthetic") is True
    assets = None
    if not synthetic:
        if "asset_paths" in config:
            assets = {name: (base / path).resolve(strict=True) for name, path in config["asset_paths"].items()}
        else:
            assets = (base / config["asset_dir"]).resolve(strict=True)
    affected_side = config.get("affected_side")
    if "fa" in config and any(method != "Boyne2025" for method in methods) and affected_side not in ("Left", "Right"):
        raise ValueError("SMATT and HCPA FA require the reviewed affected side")
    paths, provenance, label, fov, t1_to_t2, mni_to_t1 = t2_chain(config["t2"], base, output)
    model = native_fa(config["fa"], base, paths["t1"]) if "fa" in config else None
    output.mkdir(parents=True)
    receipt = {"schema": SCHEMA, "status": "COMPUTED_REVIEW_PENDING", "synthetic": synthetic,
               "config": file_record(config_path), "t2": provenance,
               "fa": model[-1] if model else None, "methods": {},
               "dlpfc": "DLPFC 待解决", "clinical_validation": False,
               "resampling": "One nearest-neighbour sample from original T2 labels and source FOV per atlas",
               "pull_direction_lps": "atlas → MNI152NLin2009cAsym → individual T1 → original T2"}

    # Process methods on the same lattice together so the bridge only runs once.
    atlases = {}
    for method in methods:
        atlases[method] = (synthetic_atlas(config["synthetic_atlases"][method], base) if synthetic
                          else real_atlas(assets, method))
    bridge = None
    for method, (target_path, reference, bank, atlas_inputs) in atlases.items():
        target = itk_grid(reference)
        method_fov = fov
        if method == "HCPA":
            # Historical HCPA used the co-gridded label header for this source mask.
            method_fov = sitk_module().Image(fov)
            method_fov.CopyInformation(label)
        template_to_t1 = mni_to_t1
        bridge_qc = None
        if METHODS[method] == MNI6:
            if synthetic:
                bridge_path = checked(config["synthetic_gen6_to_2009c_pull"], base)
                bridge = sitk_module().ReadTransform(str(bridge_path))
                bridge_qc = {"synthetic_transform": file_record(bridge_path)}
            elif bridge is None:
                bridge_path = locked_file(assets, BRIDGE_SOURCE)
                support = np.zeros(reference.shape, dtype=bool)
                for candidate, (_, _, related, _) in atlases.items():
                    if METHODS[candidate] == MNI6:
                        for tract in related.values():
                            support |= tract
                bridge, bridge_qc = invert_template_bridge(
                    sitk_module().ReadTransform(str(bridge_path)), target, support.T, reference.affine)
                bridge_qc = dict(bridge_qc, source=file_record(bridge_path))
                receipt["gen6_bridge"] = bridge_qc
            else:
                bridge_qc = receipt["gen6_bridge"]
            template_to_t1 = compose_pull(mni_to_t1, bridge)
        lesion, covered = map_pair(label, method_fov, target, compose_pull(t1_to_t2, template_to_t1))
        if np.any(lesion & ~covered):
            raise ValueError("Mapped lesion extends outside the original T2 FOV")
        z = reference.affine[2, 3] + reference.affine[2, 2] * np.arange(reference.shape[2])
        union = np.logical_or.reduce(list(bank.values()))
        fa = sample_fa(model, template_to_t1, reference, union) if model is not None else None
        results = calculate_bank(bank, lesion, covered, z, method, fa, affected_side)
        folder = output / method
        folder.mkdir()
        output_files = {
            "lesion": write_mask(folder / "lesion.nii.gz", lesion, target),
            "source_fov": write_mask(folder / "source_fov.nii.gz", covered, target),
        }
        for name, rows in results.items():
            if rows:
                table = folder / f"{name}.csv"
                write_csv(table, rows)
                output_files[name] = file_record(table)
        figure = folder / "report.html"
        if method == "HCPA" and not synthetic:
            sitk = sitk_module()
            normalized = sitk.ReadImage(str(paths["mni"]))
            background = sitk.GetArrayFromImage(resample(
                normalized, target, sitk.Transform(3, sitk.sitkIdentity))).T
        else:
            background = np.asarray(reference.dataobj)
        render(figure, method, reference.affine, background, bank, lesion, covered, results)
        output_files["report"] = file_record(figure)
        receipt["methods"][method] = {
            "template": METHODS[method], "target_reference": file_record(target_path),
            "shape": reference.shape, "affine_ras": reference.affine.tolist(), "atlas_inputs": atlas_inputs,
            "source_fov_header": "original_T2_label" if method == "HCPA" else "original_T2_image",
            "native_positive_lesion_missing_in_target": bool(
                np.any(sitk_module().GetArrayViewFromImage(label)) and not lesion.any()),
            "gen6_bridge_qc": bridge_qc, "outputs": output_files, "summary": results["whole"],
        }
    save_new(output / "receipt.json", receipt)
    return receipt
