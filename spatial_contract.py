"""Evidence shared by the tractography and lesion-load checkout commands."""
from hashlib import file_digest
import json
import math
from pathlib import Path


MNI2009C = "MNI152NLin2009cAsym"
MNI6 = "MNI152NLin6Asym"
BRIDGE_SOURCE = "tpl-MNI152NLin2009cAsym_from-MNI152NLin6Asym_mode-image_xfm.h5"
BRIDGE_SHA256 = "2e3869a07b96aec406e0419ca2e434afc54882d37cc212b933b139d1b63a4dfe"
REVIEW_SCHEMA = "neo-spatial-input-review-1"


def digest(path):
    with Path(path).open("rb") as stream:
        return file_digest(stream, "sha256").hexdigest()


def checked(record, base):
    path = (base / record["path"]).resolve(strict=True)
    if not path.is_file() or digest(path) != record["sha256"].lower():
        raise ValueError(f"Input SHA-256 mismatch: {path}")
    return path


def visit_identity(visit):
    keys = {"dataset", "subject", "session"}
    if not isinstance(visit, dict) or set(visit) != keys or not all(
        isinstance(value, str) and value and value == value.strip() for value in visit.values()
    ):
        raise ValueError("Provide the actual dataset, subject and acquisition session in visit")
    return visit


def locked_template_bridge(record, base):
    """Read only the pinned 2009c → Gen6 LPS point map."""
    import SimpleITK as sitk

    if record["sha256"].lower() != BRIDGE_SHA256:
        raise ValueError("The template bridge must be the locked 2009c → Gen6 LPS point map")
    path = checked(record, base)
    return sitk.ReadTransform(str(path))


def require_review(record, base, *, visit, scope, inputs, directions, branch=None, landmark_errors=None):
    """Check a reviewer's decision against the actual execution contract."""
    path = checked(record, base)
    review = json.loads(path.read_text(encoding="utf-8-sig"))
    if review.get("schema") != REVIEW_SCHEMA or review.get("decision") != "APPROVED":
        raise ValueError("Spatial input review has no supported approval")
    if review.get("scope") != scope or visit_identity(review.get("visit")) != visit_identity(visit):
        raise ValueError("Spatial input review belongs to a different visit or use")
    if review.get("branch") != branch or review.get("point_directions") != directions:
        raise ValueError("Spatial input review does not match the branch or point directions")
    declared = review.get("input_sha256")
    if not isinstance(declared, dict):
        raise ValueError("Spatial input review requires an input SHA-256 table")
    for role, expected in inputs.items():
        if declared.get(role) != expected:
            raise ValueError(f"Spatial input review does not bind the actual {role}")
    for key in ("reviewer", "reviewed_utc", "basis"):
        if not isinstance(review.get(key), str) or not review[key].strip():
            raise ValueError(f"Spatial input review requires {key}")
    evidence = review.get("evidence")
    if not isinstance(evidence, dict) or not evidence:
        raise ValueError("Spatial input review requires evidence of the decision")
    for entry in evidence.values():
        checked(entry, path.parent)
    limits = review.get("landmark_max_mm")
    if landmark_errors is not None:
        if not isinstance(limits, dict) or set(limits) != set(landmark_errors):
            raise ValueError("Record the reviewer-approved landmark error limit for every segment")
        for role, error in landmark_errors.items():
            limit = limits[role]
            if type(limit) not in (float, int) or not math.isfinite(limit) or limit <= 0:
                raise ValueError(f"Invalid landmark limit for {role}")
            if error["max_mm"] > limit:
                raise ValueError(f"Landmark error exceeds the reviewed limit for {role}")
    return {"path": str(path), "sha256": record["sha256"].lower(),
            "visit": visit, "scope": scope, "decision": review["decision"], "branch": branch,
            "input_sha256": inputs, "point_directions": directions, "landmark_max_mm": limits}
