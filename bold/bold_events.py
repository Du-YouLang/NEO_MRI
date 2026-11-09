
"""Version a scan-specific events override without editing any source or clearance."""
from io import BytesIO
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd

import bold_workflow as w
from bold_runs import SCHEMA


def stage_events(config_path, run_id, source, *, operator, reason, confirmed=False):
    config = w.load(config_path)
    if config["schema"] != SCHEMA:
        raise w.StageError("Event overrides require the run-aware BOLD schema")
    if confirmed is not True or not operator.strip() or not reason.strip():
        raise w.StageError("Confirm the event source and record the operator and reason")
    run = next((row for row in config["analysis"]["runs"] if row["id"]==run_id),None)
    if run is None:
        raise w.StageError("The target run is not in the configuration")
    path = Path(source).resolve(strict=True)
    suffix = path.suffix.lower()
    contents = path.read_bytes()
    if suffix==".json":
        values = json.loads(contents.decode("utf-8-sig"))
        if not isinstance(values,list) or not all(isinstance(row,dict) for row in values):
            raise w.StageError("Event JSON must contain a list of event objects")
        frame = pd.DataFrame(values)
    elif suffix in (".tsv",".csv"):
        frame = pd.read_csv(BytesIO(contents),sep="\t" if suffix==".tsv" else ",",encoding="utf-8-sig")
    else:
        raise w.StageError("Use an events TSV, CSV or JSON")
    if not {"trial_type","onset","duration"} <= set(frame) or frame.empty:
        raise w.StageError("Events require trial_type, onset and duration")
    valid_types = frame["trial_type"].map(lambda value:isinstance(value,str) and bool(value.strip()))
    timing = frame[["onset","duration"]].to_numpy(float)
    if not valid_types.all() or not np.isfinite(timing).all() or (timing[:,1]<0).any():
        raise w.StageError("Invalid event condition, onset or duration")
    tsv = frame.to_csv(sep="\t",index=False,lineterminator="\n").encode("utf-8")
    old = w.resolve(run["events"]) if run.get("events") else None
    request = dict(run=run_id,source=str(path),source_sha256=hashlib.sha256(contents).hexdigest(),
                   old_events=str(old) if old else None,old_events_sha256=w.sha(old) if old and old.is_file() else None,
                   new_events_sha256=hashlib.sha256(tsv).hexdigest(),operator=operator,reason=reason,
                   source_confirmed_by_operator=True,scan_specific_timing_clearance="MUST_REVERIFY")
    dest = w.workspace(config)/"events_overrides"/run_id/w.digest(request)
    with w.lock(config):
        if (dest/"receipt.json").is_file():
            proof = w.read(dest/"receipt.json")
            w.check_files(dest,proof["files"],True)
            return proof
        if dest.exists():
            raise w.StageError(f"Incomplete events import preserved: {dest}")
        dest.mkdir(parents=True)
        staged = dest/"events.tsv"
        staged.write_bytes(tsv)
        receipt = dict(**request,staged_events=str(staged),created_UTC=w.now(),files=[w.descriptor(staged,dest)])
        w.write(dest/"receipt.json",receipt)
        return receipt
