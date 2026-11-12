"""Portable configuration, tracking commands, screening receipts and exports."""
from base64 import b64encode
from html import escape
from io import BytesIO
import json
from pathlib import Path
import re
import subprocess
import os
import time

import nibabel as nib
import numpy as np

from .geometry import MNI6, atlas_pull, checked, digest, load_mask, load_reference
from .screening import GATE_ROLES, Screening, core_mask, read_tck_header, slab_indices

MODELS = {"Tensor_Det": "Tensor_Det", "ST_iFOD2": "iFOD2"}
PARCELS = ("M1", "PMd", "PMv", "S1", "SMA", "preSMA", "S2", "PPC")


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def write(path, data):
    with Path(path).open("x", encoding="utf-8") as stream:
        json.dump(data, stream, indent=2, allow_nan=False)
        stream.write("\n")


def record(path):
    return {"path": str(path.resolve()), "sha256": digest(path)}


def token(value):
    if not isinstance(value, str) or not re.fullmatch(r"[A-Za-z0-9_-]{1,80}", value):
        raise ValueError(f"Invalid run identifier: {value!r}")
    return value


def positive(value, name):
    if type(value) is not int or value < 1:
        raise ValueError(f"{name} must be a positive integer")
    return value


def config(path):
    path = Path(path).resolve(strict=True)
    data = read(path)
    if data.get("schema") != "neo-tractography-1" or not data.get("runs"):
        raise ValueError("A neo-tractography-1 configuration with explicit runs is required")
    seen = set()
    for run in data["runs"]:
        key = token(run["id"])
        if key in seen:
            raise ValueError(f"Duplicate run ID: {key}")
        seen.add(key)
        if run["origin"] not in PARCELS or run["side"] not in ("L", "R") or run["model"] not in MODELS:
            raise ValueError("Declare a supported origin, hemisphere and model; M2 summarizes SMA and preSMA")
        if "plic_waypoints" in run and run["model"] != "ST_iFOD2":
            raise ValueError("The reviewed PLIC+2 profile uses ST-CSD/iFOD2")
    return data, path.parent


def mm(value, name):
    if type(value) not in (int, float) or not np.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def tracking_plan(path, output):
    data, base = config(path)
    ref = load_reference(checked(data["reference"], base))
    jobs = []
    for run in data["runs"]:
        name = run["id"]
        model = checked(run["tracker_input"], base)
        seeds = checked(run["seed_mask"], base)
        support = checked(run["tracking_support"], base)
        if not load_mask(seeds, ref).any() or not load_mask(support, ref).any():
            raise ValueError(f"Empty seed or tracking support: {name}")
        attempts = positive(run.get("attempts", 200000), "attempts")
        threads = positive(run.get("threads", 8), "threads")
        rng = positive(run["rng_seed"], "rng_seed")
        argv = [
            "tckgen", str(model), str(output / f"{name}.tck"), "-algorithm", MODELS[run["model"]],
            "-seed_image", str(seeds), "-seed_unidirectional", "-output_seeds",
            str(output / f"{name}.seeds.txt"), "-mask", str(support), "-select", "0",
            "-seeds", str(attempts), "-cutoff", "0.10", "-seed_cutoff", "0.10",
            "-minlength", "5", "-maxlength", "250", "-info", "-nthreads", str(threads),
        ]
        if run["model"] == "Tensor_Det":
            argv += ["-config", "BZeroThreshold", "50"]
        jobs.append({"id": name, "argv": argv, "attempts": attempts, "rng_seed": rng,
                     "shell_b": mm(run["shell_b"], "shell_b") if run["model"] == "Tensor_Det" else None})
    return {"version": "3.0.4", "jobs": jobs}


def fit_plan(path, output):
    """Plan the reviewed OLS metrics and single-tissue CSD recipes."""
    path = Path(path).resolve(strict=True)
    data = read(path)
    if data.get("schema") != "neo-tractography-fit-1":
        raise ValueError("A neo-tractography-fit-1 configuration is required")
    ref = load_reference(checked(data["reference"], path.parent))
    mask = checked(data["fit_mask"], path.parent)
    if not load_mask(mask, ref).any():
        raise ValueError("Empty fit mask")
    threads = str(positive(data.get("threads", 4), "threads"))
    low = checked(data["tensor_shell_mif"], path.parent)
    high = checked(data["csd_shell_mif"], path.parent)
    bzero = ["-config", "BZeroThreshold", "50", "-nthreads", threads]
    output = output.resolve()
    mask_flag = ["-mask", str(mask)]
    jobs = [
        ["dwi2tensor", str(low), str(output / "tensor.mif"), *mask_flag, "-ols", "-iter", "0", *bzero],
        ["tensor2metric", str(output / "tensor.mif"), *mask_flag,
         "-fa", str(output / "FA.nii.gz"), "-adc", str(output / "MD.nii.gz"),
         "-ad", str(output / "AD.nii.gz"), "-rd", str(output / "RD.nii.gz"), *bzero],
        ["dwi2response", "tournier", str(high), str(output / "response.txt"), *mask_flag,
         "-lmax", "6", "-scratch", str(output / "response_work"), *bzero],
        ["dwi2fod", "csd", str(high), str(output / "response.txt"), str(output / "wmfod.mif"),
         *mask_flag, "-lmax", "6", *bzero],
    ]
    return {"version": "3.0.4", "jobs": jobs,
            "shell_b": {kind: mm(data[f"{kind}_shell_b"], f"{kind}_shell_b") for kind in ("tensor", "csd")},
            "inputs": {"tensor_shell_mif": data["tensor_shell_mif"], "csd_shell_mif": data["csd_shell_mif"],
                       "fit_mask": data["fit_mask"], "reference": data["reference"]}}


def mrtrix_version(executable):
    result = subprocess.run([executable, "-version"], capture_output=True, text=True, check=True, timeout=15)
    version = (result.stdout + result.stderr).strip()
    if not re.search(r"(?<![\d.])3\.0\.4(?![\d.])", version):
        raise ValueError(f"Expected MRtrix 3.0.4 for this workflow: {version}")
    return version


def mrinfo(path, *options):
    result = subprocess.run(["mrinfo", str(path), *options, "-config", "BZeroThreshold", "50"],
                            capture_output=True, text=True, check=True, timeout=30)
    return result.stdout


def check_tracker_grid(path, ref):
    """Compare the physical grid, allowing storage-axis permutations."""
    lines = mrinfo(path, "-size", "-spacing", "-transform").splitlines()
    shape = tuple(map(int, lines[0].split()))
    spacing = np.fromstring(lines[1], sep=" ")
    affine = np.loadtxt(lines[2:])
    if len(shape) != 4 or spacing.shape != (4,) or affine.shape != (4, 4):
        raise ValueError(f"Expected a four-dimensional model input with a physical grid: {path}")
    affine[:3, :3] *= spacing[:3]
    orientation = nib.orientations.io_orientation(affine)
    canonical = affine @ nib.orientations.inv_ornt_aff(orientation, shape[:3])
    spatial_shape = tuple(shape[int(i)] for i in np.argsort(orientation[:, 0]))
    if spatial_shape != ref.shape or not np.allclose(canonical, ref.affine, rtol=0, atol=.001):
        raise ValueError(f"MRtrix input and corrected b0 do not share a physical grid: {path}")
    return shape[3]


def check_gradients(path, volumes, shell, minimum_directions):
    gradients = np.loadtxt(mrinfo(path, "-dwgrad").splitlines(), ndmin=2)
    if gradients.shape != (volumes, 4) or not np.isfinite(gradients).all():
        raise ValueError(f"Gradient count or finite values disagree with the DWI input: {path}")
    b = gradients[:, 3]
    dw = b >= 50
    if ((b < 0).any() or dw.all() or dw.sum() < minimum_directions
            or not np.all(np.abs(b[dw] - shell) <= 50)):
        raise ValueError(f"Expected b0 plus one declared shell within 50 s/mm²: {path}")
    directions = gradients[dw, :3]
    if not np.allclose(np.linalg.norm(directions, axis=1), 1, atol=.01, rtol=0):
        raise ValueError(f"Nonunit diffusion directions: {path}")
    if len(np.unique(np.round(directions[:, :, None] * directions[:, None, :], 3),
                     axis=0)) < minimum_directions:
        raise ValueError(f"Too few distinct diffusion directions: {path}")


def fit(path, output):
    output = output.resolve()
    plan = fit_plan(path, output)
    executables = {"mrinfo", *(argv[0] for argv in plan["jobs"])}
    versions = {command: mrtrix_version(command) for command in sorted(executables)}
    base = path.resolve().parent
    ref = load_reference(checked(plan["inputs"]["reference"], base))
    for kind, minimum in (("tensor", 6), ("csd", 28)):
        item = checked(plan["inputs"][f"{kind}_shell_mif"], base)
        check_gradients(item, check_tracker_grid(item, ref), plan["shell_b"][kind], minimum)
    output.mkdir(parents=True)
    write(output / "plan.json", {**plan, "versions": versions})
    for argv in plan["jobs"]:
        with (output / f"{argv[0]}.log").open("xb") as log:
            subprocess.run(argv, stdout=log, stderr=subprocess.STDOUT, check=True)
    names = ("tensor.mif", "FA.nii.gz", "MD.nii.gz", "AD.nii.gz", "RD.nii.gz", "response.txt", "wmfod.mif")
    result = {"files": {name: record(output / name) for name in names}}
    write(output / "fit_receipt.json", result)
    return result


def absolute_records(value, base):
    """Rebase explicit file records when moving a generated config to another directory."""
    if isinstance(value, dict):
        if "path" in value and "sha256" in value:
            return {**value, "path": str((base / value["path"]).resolve())}
        return {key: absolute_records(item, base) for key, item in value.items()}
    if isinstance(value, list):
        return [absolute_records(item, base) for item in value]
    return value


def track(path, output):
    output = output.resolve()
    data, base = config(path)
    plan = tracking_plan(path, output)
    versions = {program: mrtrix_version(program) for program in ("mrinfo", "tckgen")}
    ref = load_reference(checked(data["reference"], base))
    for job in plan["jobs"]:
        item = Path(job["argv"][1])
        size = check_tracker_grid(item, ref)
        if job["shell_b"] is not None:
            check_gradients(item, size, job["shell_b"], 6)
        elif size != 28:
            raise ValueError(f"Expected a single-tissue SH lmax=6 FOD with 28 coefficients: {item}")
    output.mkdir(parents=True)
    write(output / "plan.json", {**plan, "observed_versions": versions})
    replay = absolute_records(data, base)
    with (output / "commands.jsonl").open("x", encoding="utf-8") as receipt:
        for job, run in zip(plan["jobs"], replay["runs"], strict=True):
            name = job["id"]
            start = time.monotonic()
            log = output / f"{name}.log"
            with log.open("xb") as stream:
                subprocess.run(job["argv"], stdout=stream, stderr=subprocess.STDOUT,
                               env={**os.environ, "MRTRIX_RNG_SEED": str(job["rng_seed"])}, check=True)
            parent = output / f"{name}.tck"
            header = read_tck_header(parent)
            progress = re.findall(rb"\[100%\]\s+(\d+) seeds,\s+(\d+) streamlines,\s+(\d+) selected",
                                  log.read_bytes())
            if not progress:
                raise ValueError(f"Missing final tracking counters: {name}")
            attempts, generated, saved = map(int, progress[-1])
            if (attempts != job["attempts"] or saved != int(header["count"])
                    or int(header["max_num_seeds"]) != attempts or header["unidirectional"] != "1"
                    or not np.isclose(float(header["threshold"]), .10)
                    or not np.isclose(float(header["init_threshold"]), .10)
                    or header.get("roi") != f"mask {job['argv'][job['argv'].index('-mask') + 1]}"):
                raise ValueError(f"Recorded tracking settings differ from the plan: {name}")
            run.update(parent=record(parent), seed_log=record(output / f"{name}.seeds.txt"),
                       attempts=attempts, saved=saved)
            receipt.write(json.dumps({"id": name, "attempts": attempts, "generated": generated,
                                      "saved": saved, "seconds": time.monotonic() - start,
                                      "parent": run["parent"], "seed_log": run["seed_log"],
                                      "log": record(log)}) + "\n")
            receipt.flush()
    target = output / "replay.local.json"
    write(target, replay)
    return {"replay_config": str(target)}


def laterality(data, base):
    if "anatomy" not in data:
        if data.get("synthetic") is not True:
            raise ValueError("An approved DWI-to-atlas physical coordinate chain is required")
        return lambda point: "L" if point[0] < 0 else "R"
    pull = atlas_pull(data["anatomy"], base, MNI6)

    def classify_side(point):
        x, y, z = map(float, point)
        return "L" if pull.TransformPoint((-x, -y, z))[0] > 0 else "R"

    return classify_side


def screen(path, output):
    data, base = config(path)
    reference = checked(data["reference"], base)
    ref = load_reference(reference)

    def mask(item):
        return load_mask(checked(item, base), ref)

    ribbon = mask(data["ribbon"])
    gate_support = mask(data["gate_support"])
    masks = {f"{role}_{side}": mask(data["gates"][f"{role}_{side}"])
             for role in GATE_ROLES for side in "LR"}
    side_at = laterality(data, base)
    slab_mm = mm(data["slab_mm"], "slab_mm")
    slabs = slab_indices(ref.shape, ref.affine, slab_mm)
    output = output.resolve()
    output.mkdir(parents=True)
    write(output / "inputs.local.json", absolute_records(data, base))
    rows = []
    for run in data["runs"]:
        attempts = positive(run["attempts"], "attempts")
        saved = run["saved"]
        if type(saved) is not int or not 0 <= saved <= attempts:
            raise ValueError(f"Invalid saved count: {run['id']}")
        parent = checked(run["parent"], base)
        seed_log = checked(run["seed_log"], base)
        plic = None
        if "plic_waypoints" in run:
            plic = {name: mask(run["plic_waypoints"][name]) for name in ("PLIC", "CP", "pons")}
        engine = Screening(ref, ribbon, mask(run["seed_mask"]), mask(run["tracking_support"]),
                           gate_support, masks, side_at, plic)
        metrics, indices, visits = engine.replay(parent, seed_log, run["side"],
                                               mm(run["ribbon_distance_mm"], "ribbon_distance_mm"), saved)
        root = output / run["id"]
        root.mkdir()
        cores = {gate: core_mask(visit, slabs) for gate, visit in visits.items()}
        np.savez_compressed(root / "selection_indices.npz", **indices)
        np.savez_compressed(root / "visits.npz", **visits)
        np.savez_compressed(root / "core_15pct.npz", **cores)
        row = {key: run[key] for key in ("id", "origin", "side", "model", "attempts", "saved", "ribbon_distance_mm")}
        row.update(parent={"path": str(parent), "sha256": run["parent"]["sha256"].lower()},
                   seed_log={"path": str(seed_log), "sha256": run["seed_log"]["sha256"].lower()},
                   metrics=metrics, gate_support=engine.gate_counts,
                   gates={gate: {"selected": len(indices[gate]), "visited_voxels": int(np.count_nonzero(visits[gate])),
                                 "core_voxels": int(cores[gate].sum())} for gate in indices},
                   files={name: digest(root / name)
                          for name in ("selection_indices.npz", "visits.npz", "core_15pct.npz")})
        rows.append(row)
    result = {"schema": "neo-tractography-result-1",
              "reference": {"path": str(reference), "sha256": data["reference"]["sha256"].lower()},
              "config_sha256": digest(output / "inputs.local.json"),
              "slab_mm": slab_mm, "core_fraction": .15, "runs": rows}
    write(output / "results.json", result)
    return result


def result_run(folder, run_id, gate):
    result = read(folder / "results.json")
    if result.get("schema") != "neo-tractography-result-1":
        raise ValueError("Unrecognized result format")
    run = next((row for row in result["runs"] if row["id"] == run_id), None)
    if run is None or gate not in run["gates"]:
        raise ValueError(f"Unavailable run or gate: {run_id} / {gate}")
    root = folder / token(run["id"])
    return result, run, root


def result_array(root, run, filename, gate):
    path = checked({"path": filename, "sha256": run["files"][filename]}, root)
    with np.load(path, allow_pickle=False) as store:
        return store[gate]


def selected_lines(parent, indices):
    if len(indices) == 0:
        return
    position = 0
    for i, line in enumerate(nib.streamlines.load(parent, lazy_load=True).streamlines):
        if i == indices[position]:
            yield line
            position += 1
            if position == len(indices):
                return
    raise ValueError("Selected indices extend past the parent tractogram")


def indices_for(root, run, gate):
    values = result_array(root, run, "selection_indices.npz", gate)
    if (values.ndim != 1 or values.dtype.kind not in "iu" or len(values) != run["gates"][gate]["selected"]
            or (len(values) and (values[0] < 0 or values[-1] >= run["saved"] or np.any(np.diff(values) <= 0)))):
        raise ValueError("Selected indices disagree with the result receipt")
    return values


def export(folder, run_id, gate, output):
    from nibabel.streamlines import LazyTractogram, TckFile
    from nibabel.streamlines.tractogram import TractogramItem

    result, run, root = result_run(folder, run_id, gate)
    ref = load_reference(checked(result["reference"], folder))
    parent = checked(run["parent"], folder)
    indices = indices_for(root, run, gate)
    visits = result_array(root, run, "visits.npz", gate)
    core = result_array(root, run, "core_15pct.npz", gate)
    if not np.array_equal(core, core_mask(visits, slab_indices(ref.shape, ref.affine, result["slab_mm"]))):
        raise ValueError("Coverage and visits disagree")
    output.mkdir(parents=True)
    np.save(output / "indices.npy", indices, allow_pickle=False)
    for filename, values in (("visits.nii.gz", visits), ("core_15pct.nii.gz", core.astype(np.uint8))):
        nib.save(nib.Nifti1Image(values, ref.affine), output / filename)
    if len(indices):
        def items():
            return (TractogramItem(line, {}, {}) for line in selected_lines(parent, indices))

        tracks = LazyTractogram.from_data_func(items)
        tracks.affine_to_rasmm = np.eye(4)
        target = output / "selected.tck"
        TckFile(tracks).save(str(target))
        if int(read_tck_header(target)["count"]) != len(indices):
            raise ValueError("Selected TCK count differs from the selected indices")
    receipt = {"run": run_id, "gate": gate, "parent_sha256": run["parent"]["sha256"],
               "counts": run["gates"][gate], "files": {f.name: digest(f) for f in output.iterdir()}}
    write(output / "export_receipt.json", receipt)
    return receipt


def report(folder, output):
    from .viewer import report_preview

    result = read(folder / "results.json")
    if result.get("schema") != "neo-tractography-result-1":
        raise ValueError("Expected a tractography results receipt")
    ref = load_reference(checked(result["reference"], folder))
    rows = []
    for run in result["runs"]:
        for gate, row in run["gates"].items():
            values = (run["id"], run["origin"], run["side"], run["model"], run["attempts"],
                      run["saved"], run["metrics"]["descending"], gate, row["selected"],
                      row["visited_voxels"], row["core_voxels"])
            rows.append("<tr>" + "".join(f"<td>{escape(str(value))}</td>" for value in values) + "</tr>")
    labels = ("运行", "起点", "侧别", "模型", "播种尝试", "保存流线", "向下流线",
              "筛选门", "候选流线", "访问体素", "15% 核心体素")
    head = "".join(f"<th scope='col'>{label}</th>" for label in labels)
    with output.open("x", encoding="utf-8") as stream:
        stream.write("""<!doctype html><html lang="zh-CN"><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>CST／CRP 追踪与筛选</title><style>
body{font:16px system-ui,sans-serif;margin:2rem;color:#172333}table{border-collapse:collapse}
th,td{padding:.55rem;border:1px solid #cad3dd;text-align:left}th{background:#e8eef5}
details{margin:.7rem 0}summary{cursor:pointer}img{display:block;width:min(100%,54rem);height:auto}
</style><h1>CST／CRP 追踪与筛选</h1>""")
        stream.write(f"<p>按原生 RAS 高度每 {result['slab_mm']:g} mm 分层，保留各层访问数达到峰值 15% 的体素。")
        stream.write("每条候选在同一体素最多贡献一次访问。</p>")
        stream.write("<p>SMA 和 preSMA 分开报告。Notebook 的 M2 入口可分别查看两组流线。</p>")
        stream.write(f"<table><thead><tr>{head}</tr></thead><tbody>{''.join(rows)}</tbody></table>")
        stream.write("<h2>各组候选在原生 b0 上的位置</h2><p>展开一组查看覆盖计数。每组显示一个解剖门。"
                     "更多门、父流线和核心图见分步 Notebook。</p>")
        for run in result["runs"]:
            gate = ("CST_PLIC_plus2_CP_Tang_pons" if "CST_PLIC_plus2_CP_Tang_pons" in run["gates"]
                    else "Tang_pons_medulla")
            run_id = token(run["id"])
            volume = result_array(folder / run_id, run, "visits.npz", gate)
            image = BytesIO()
            report_preview(ref, volume).savefig(image, format="png", dpi=75)
            embedded = b64encode(image.getvalue()).decode("ascii")
            stream.write(f"<details><summary>{escape(run_id)} · {escape(gate)}</summary>"
                         f"<img alt='原生 b0 三切面覆盖计数' loading='lazy' src='data:image/png;base64,{embedded}'></details>")
        stream.write("</html>")
    return {"output": str(output)}
