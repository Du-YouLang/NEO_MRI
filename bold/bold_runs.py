"""Run identity and sample-clock contracts for the versioned BOLD workflow."""
from dataclasses import dataclass
from pathlib import Path
import math
import re

import nibabel as nib
import numpy as np

import bold_workflow as workflow
from bold_workflow import StageError, TimingBlocked

SCHEMA = "neo-bold-stages-2"
KERNELS = (0, 2, 3)
SAFE_ID = re.compile(r"[A-Za-z0-9_-]+")
RUN_FILES = ("bold", "bold_json", "mask", "confounds", "confounds_json")


def require_id(value):
    if not isinstance(value, str) or not SAFE_ID.fullmatch(value):
        raise StageError(f"Invalid identifier: {value!r}")


def finite_number(value):
    return type(value) in (int, float) and math.isfinite(value)


def contrast_for(config, run):
    return run.get("contrast", config["analysis"]["statistics"]["contrast"])


def validate(config):
    analysis = config["analysis"]
    stats = analysis["statistics"]
    if stats["hrf_model"] != "spm + derivative":
        raise StageError("BOLD v2 requires the static SPM GLM with temporal derivative")
    if stats["noise_model"] not in ("ar1", "ols"):
        raise StageError("BOLD v2 noise_model must be ar1 or ols")
    if type(stats["compcor_components"]) is not int or stats["compcor_components"] < 0:
        raise StageError("compcor_components must be a nonnegative integer")
    if any(not finite_number(stats[k]) or stats[k] < 0 for k in ("fd_threshold_mm", "dvars_threshold")):
        raise StageError("Invalid motion spike thresholds")
    if not finite_number(stats["fdr_q"]) or not 0 < stats["fdr_q"] < 1:
        raise StageError("Invalid FDR level")
    kernels = analysis["kernels_mm"]
    if not isinstance(kernels, list) or not kernels or any(type(k) is not int or k not in KERNELS for k in kernels):
        raise StageError("Choose explicit analysis kernels from 0, 2, 3 mm")
    if len(kernels) != len(set(kernels)) or stats["smoothing_fwhm_mm"] not in kernels:
        raise StageError("Kernels must be unique and include the declared primary kernel")
    scan = analysis["double_scanner"]
    if set(scan) != {"lags_seconds", "duration_ratios", "selection"}:
        raise StageError("DoubleScanner requires explicit lags, ratios and selection")
    if scan["selection"] not in ("raw_r2", "white_partial_r2"):
        raise StageError("Unknown DoubleScanner selection rule")
    for name in ("lags_seconds", "duration_ratios"):
        values = scan[name]
        if not isinstance(values, list) or not values or not all(finite_number(v) for v in values):
            raise StageError(f"Invalid DoubleScanner {name}")
        if len(values) != len(set(values)):
            raise StageError(f"Duplicate DoubleScanner {name}")
    from bold_ab import validate_policy
    validate_policy(analysis)
    maximum = 1215 if analysis.get("double_scanner_ab") is not None else 1024
    if min(scan["duration_ratios"])<=0 or len(scan["lags_seconds"])*len(scan["duration_ratios"])>maximum:
        raise StageError(f"Duration ratios must be positive; grid limit is {maximum}")
    fir = analysis["fir_delays_tr"]
    if not isinstance(fir, list) or not fir or any(type(v) is not int or v < 0 for v in fir) or len(set(fir)) != len(fir):
        raise StageError("FIR requires distinct nonnegative integer lags in TR units")

    if config["preprocessing"]["mode"] not in ("existing","engineering_subset"):
        raise StageError("BOLD v2 registers completed preprocessing; run the reviewed preprocessing entry separately")
    if config["preprocessing"]["mode"]=="existing":
        state = workflow.workspace(config)
        data = workflow.data_root(config)
        parts = state.relative_to(data).parts if state.is_relative_to(data) else ()
        if len(parts)<4 or parts[2:4]!=("results","task_fMRI"):
            raise StageError("BOLD v2 outputs belong under the individual visit's results/task_fMRI")
    runs = {run["id"]: run for run in analysis["runs"]}
    if not runs:
        raise StageError("Configure at least one BOLD run")
    if len(runs) != len(analysis["runs"]):
        raise StageError("Run identifiers must be unique")
    if len({(workflow.derivatives(config)/run["bold"]).resolve() for run in runs.values()}) != len(runs):
        raise StageError("The same BOLD NIfTI cannot represent two runs")
    acquisitions = [re.sub(r"_(?:echo|part)-[^_]+","",stem(Path(run["bold"]))) for run in runs.values()]
    visits = {tuple(part for part in name.split("_") if part.startswith(("sub-","ses-"))) for name in acquisitions}
    if len(visits)!=1 or not all(visit and visit[0].startswith("sub-") for visit in visits):
        raise StageError("One BOLD config must identify one BIDS participant and visit")
    if len(acquisitions) != len(set(acquisitions)):
        raise StageError("Different spatial derivatives of one acquisition cannot be independent runs")
    for run_id, run in runs.items():
        require_id(run_id)
        if run.get("task_id") is not None:
            require_id(run["task_id"])
        if run.get("estimand_id") is not None:
            require_id(run["estimand_id"])
        if run.get("reference_space") is not None:
            require_id(run["reference_space"])
            tokens = {token.removeprefix("space-") for token in Path(run["bold"]).name.split("_")
                      if token.startswith("space-")}
            if tokens and tokens != {run["reference_space"]}:
                raise StageError("Declared reference space differs from the BOLD filename")
        if run.get("stc") not in (None, "on", "off"):
            raise StageError("STC is on, off, or unresolved")
        if run.get("sdc") not in (None,"off","pepolar","fieldmap","syn"):
            raise StageError("SDC route must be off, pepolar, fieldmap, syn or unresolved")
        if bool(run.get("viewer_anatomy"))!=bool(run.get("viewer_anatomy_mask")):
            raise StageError("The optional 4D viewer requires both T1 and its mask")
        contrast = contrast_for(config, run)
        if not isinstance(contrast, dict) or not contrast or not all(
            isinstance(k, str) and k and finite_number(v) for k, v in contrast.items()
        ) or not any(contrast.values()):
            raise StageError(f"Invalid static contrast for {run_id}")

    products = set(runs)
    for group in analysis.get("groups", []):
        require_id(group["id"])
        require_id(group["estimand_id"])
        if group["id"] in products:
            raise StageError("Repeated output identifier")
        products.add(group["id"])
        members = group["runs"]
        if len(members) < 2 or len(set(members)) != len(members) or not set(members) <= set(runs):
            raise StageError("A fixed-effect group requires distinct configured runs")
        if group.get("same_estimand_confirmed") is not True:
            raise StageError("Fixed effects require explicit same-estimand confirmation")
        selected = [runs[run_id] for run_id in members]
        if (not selected[0].get("task_id")
            or any(run.get("task_id") != selected[0]["task_id"] for run in selected)
            or any(run.get("estimand_id") != group["estimand_id"] for run in selected)
            or any(contrast_for(config, run) != contrast_for(config, selected[0]) for run in selected)):
            raise StageError("Fixed effects cannot pool different tasks, estimands or contrasts")
    comparisons = set()
    for pair in analysis.get("comparisons", []):
        require_id(pair["id"])
        if pair["id"] in comparisons or pair["left"] == pair["right"] or not {pair["left"], pair["right"]} <= products:
            raise StageError("Invalid or duplicate post-fit comparison")
        if pair.get("model","GLM_TD") not in ("GLM_TD","DoubleScanner","DoubleScanner_AB"):
            raise StageError("Beta comparison supports GLM_TD, DoubleScanner or DoubleScanner_AB")
        if pair.get("model") in ("DoubleScanner","DoubleScanner_AB") and not {pair["left"],pair["right"]} <= set(runs):
            raise StageError("DoubleScanner comparison requires two individual run fits")
        comparisons.add(pair["id"])


def file_for(config, registered, run, key):
    name = run.get(key)
    if not name:
        raise TimingBlocked(f"{run['id']}: missing {key}")
    base = workflow.derivatives(config).resolve()
    path = (base / name).resolve()
    if not path.is_relative_to(base) or path.relative_to(base).as_posix() not in registered:
        raise StageError(f"{run['id']}: {key} is outside the registered derivative set")
    if not path.is_file():
        raise TimingBlocked(f"{run['id']}: missing {key}: {path}")
    return path


def stem(bold):
    return re.split(r"_(?:space-|desc-)", bold.name, maxsplit=1)[0]


@dataclass(frozen=True)
class Clock:
    times: np.ndarray
    tr: float
    start: float
    stc: str


def clock_for(run, image, sidecar, source_tr):
    tr = run.get("tr_seconds")
    if not finite_number(tr) or tr <= 0 or not finite_number(source_tr) or source_tr <= 0:
        raise TimingBlocked(f"{run['id']}: verified positive TR required")
    if len(image.shape) != 4 or image.shape[3] < 2:
        raise StageError("Expected a separate four-dimensional BOLD NIfTI per run")
    header_tr = image.header.get_zooms()[3]
    time_units = image.header.get_xyzt_units()[1]
    conversion = {"sec": 1., "msec": .001, "usec": .000001, "unknown": 1.}
    if time_units not in conversion:
        raise StageError("BOLD NIfTI has an unsupported temporal unit")
    header_tr *= conversion[time_units]
    if not np.isclose(source_tr, tr, rtol=1e-5, atol=1e-6) or not np.isclose(header_tr, tr, rtol=1e-5, atol=1e-6):
        raise StageError("Nominal, source or NIfTI TR conflict")
    derivative_tr = sidecar.get("RepetitionTime", tr)
    if not finite_number(derivative_tr) or not np.isclose(derivative_tr, tr, rtol=1e-5, atol=1e-6):
        raise StageError("Derivative sidecar TR conflict")

    requested = run.get("stc")
    observed = sidecar.get("SliceTimingCorrected")
    if requested not in ("on", "off") or type(observed) is not bool:
        raise TimingBlocked(f"{run['id']}: declared and observed STC are required")
    if (requested == "on") != observed:
        raise StageError("Observed STC differs from the requested run profile")
    start = sidecar.get("StartTime", 0. if not observed else None)
    if not finite_number(start) or not 0 <= start < tr or (not observed and start != 0):
        raise TimingBlocked(f"{run['id']}: invalid or missing derivative StartTime")
    return Clock(start + np.arange(image.shape[3], dtype=float)*tr, float(tr), float(start), requested)


def inventory(config):
    """Header-only inventory. Acquisition labels do not decide scientific task equivalence."""
    base = workflow.derivatives(config)
    if not base.is_dir():
        base = workflow.source_root(config)/"bids"
    assigned = {(workflow.derivatives(config)/r["bold"]).resolve(): r for r in config["analysis"]["runs"]}
    rows = []
    for path in sorted(base.rglob("*_bold.nii.gz")):
        if "_desc-" in path.name and not path.name.endswith("_desc-preproc_bold.nii.gz"):
            continue
        sidecar = path.with_name(path.name.removesuffix(".nii.gz")+".json")
        run = assigned.get(path.resolve(), {})
        tokens = dict(token.split("-",1) for token in stem(path).split("_") if "-" in token)
        missing = [label for label,value in (("BOLD JSON",sidecar.is_file()),("scientific task_id",run.get("task_id")),("events",run.get("events"))) if not value]
        try:
            image = nib.load(path)
            frames = image.shape[3] if len(image.shape) == 4 else None
            tr = float(image.header.get_zooms()[3]) if frames is not None else None
            if frames is None:
                missing.append("4D BOLD")
            if tr is not None and not math.isfinite(tr):
                missing.append("finite TR")
                tr = None
        except (OSError, nib.filebasedimages.ImageFileError):
            frames = tr = None
            missing.append("readable NIfTI header")
        rows.append(dict(path=str(path), acquisition_entities=tokens, frames=frames, header_TR=tr,
                         sidecar=str(sidecar) if sidecar.is_file() else None, configured_run=run.get("id"),
                         scientific_task=run.get("task_id"), needs_review=missing))
    return dict(status="HEADER_INVENTORY", candidates=rows, task_equivalence_inferred=False)


def run_identity(config, pre):
    """Content-bound keys. A missing file is an explicit pending dependency."""
    rows = []
    for run in config["analysis"]["runs"]:
        timing = run.get("timing", config["analysis"]["timing"])
        row = dict(run=run, timing=timing)
        for key, value in (("events",run.get("events")), ("TR_source",run.get("tr_evidence")), ("timing_evidence",timing.get("evidence")), ("method_review",run.get("method_review"))):
            path = workflow.resolve(value) if value else None
            row[key] = workflow.sha(path) if path and path.is_file() else None
        rows.append(row)
    engine = workflow.CODE
    spec = dict(schema=SCHEMA, case_id=config["case_id"], preprocessing_fingerprint=pre["fingerprint"], runs=rows,
                settings={k:v for k,v in config["analysis"].items() if k not in ("runs","timing","comparisons")},
                engine_sha256={name:workflow.sha(engine/name) for name in
                               ("bold_runs.py","bold_models.py","bold_stage.py","bold_ab.py","bold_workflow.py")},
                thread_policy=workflow.THREADS)
    return workflow.digest(spec), spec
