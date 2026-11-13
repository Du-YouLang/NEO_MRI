
"""Independent GLS/bootstrap reference and end-to-end portable BOLD checks."""
import base64
import json
from pathlib import Path
import sys

import nibabel as nib
import numpy as np
import pytest
from scipy.stats import t

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"bold"))
import bold_workflow as w
from bold_ab import (LAGS, RATIOS, PARAMETERS, classify, corrected_family,
                     bootstrap_choices, chosen_gls, stability_frequencies)
from bold_models import whiten
from bold_demo import make_dataset


def test_fallback_order_and_whole_family():
    p = PARAMETERS
    default, = np.flatnonzero((p==[0,1]).all(axis=1))
    interior, = np.flatnonzero((p==[1,1.1]).all(axis=1))
    winner = np.array([0,44,1170,1214,45,89,interior,interior,interior])
    partial = np.array([0.]*6+[.05,.049,.2])
    gain = np.array([0.]*6+[.01,.2,.009])
    edge,reason,final = classify(p,winner,partial,gain)
    np.testing.assert_array_equal(reason,[1,1,2,2,3,4,0,5,5])
    np.testing.assert_array_equal(final,np.where(reason==0,winner,default))
    _,q = corrected_family(np.array([0.,1.,5.,8.]),119)
    assert (q["audit_1215"]>=q["A"]).all() and (q["A"]>=q["B"]).all()
    assert q["B"][-1]<.05<q["A"][1]
    choices = np.tile([interior]*5,(1,1)).T
    choices[0,0]=0
    freq = stability_frequencies(p,edge,choices,np.array([interior]))
    np.testing.assert_array_equal(freq["near_inside"],[.8])


def test_raw_bootstrap_against_direct_gls():
    rng = np.random.default_rng(708)
    n = 40
    z = np.column_stack((rng.normal(size=(n,2)),np.ones(n)))
    g = rng.normal(size=(n,9))
    y = rng.normal(size=(n,3))+g[:,[2,3,4]]*2
    rho = np.array([.42,-.25,.42])
    anchor = np.array([2,3,4])
    draws = rng.integers(0,n,size=(11,n))
    got = bootstrap_choices(y,z,g,rho,anchor,draws,"cpu")
    got_sse,got_beta,got_t,df = chosen_gls(y,rho,z,g,anchor)
    for k,r in enumerate(rho):
        design = np.column_stack((z,g[:,anchor[k]]))
        W = np.eye(n)-np.eye(n,k=-1)*r
        wy = W@y[:,k]
        wx = W@design
        estimate = np.linalg.lstsq(wx,wy,rcond=None)[0]
        fitted = wx@estimate
        error = wy-fitted
        v = np.linalg.inv(wx.T@wx)[-1,-1]*(error@error)/df
        np.testing.assert_allclose([got_sse[k],got_beta[k],got_t[k]],
                                    [sum((y[:,k]-design@estimate)**2),estimate[-1],estimate[-1]/np.sqrt(v)],
                                    atol=1e-10,rtol=0)
        for row,draw in enumerate(draws):
            x = fitted+(error-error.mean())[draw]*np.sqrt(n/df)
            native = np.linalg.solve(W,x)
            errors = []
            for candidate in g.T:
                a = np.column_stack((z,candidate))
                beta = np.linalg.lstsq(W@a,x,rcond=None)[0]
                errors.append(np.sum((native-a@beta)**2))
            assert got[row,k]==int(np.argmin(errors))


@pytest.fixture
def dataset(tmp_path):
    return make_dataset(tmp_path)


def test_end_to_end(dataset):
    config = w.load(dataset)
    assert len(w.dispatch(config,"inventory")["candidates"])==2
    assert w.dispatch(config,"preprocess")["fmriprep_invocations"]==0
    assert w.dispatch(config,"audit")["ready_runs"]==2
    result = w.dispatch(config,"glm")
    assert result["glm_fits"]==2
    fit = Path(result["output"])
    folder = fit/"run_1/fwhm_2mm"
    assert (folder/"fir_joint_q_bh.nii.gz").is_file()
    with np.load(folder/"double_scanner_ab.npz") as f:
        assert f["parameters"].shape==(1215,2)
        assert f["draws"].shape==(1000,110)
        assert np.all(f["q_A"]>=f["q_B"])
        np.testing.assert_allclose(f["p"],2*t.sf(abs(f["t"]),f["df"]))
        assert f["B"].any()
        for name in ("fwhm_0mm","fwhm_3mm"):
            assert not (fit/"run_1"/name/"double_scanner_ab.npz").exists()
    assert w.dispatch(config,"compare",full=True)["glm_fits"]==0
    report = Path(w.dispatch(config,"report",full=True)["output"])
    assert "探索性" in report.read_text(encoding="utf-8")
    viewer = report.parent/"ab/run_1/viewer.html"
    assert viewer.is_file() and (viewer.parent/"raw_argmax.svg").is_file()
    payload = json.loads((viewer.parent/"viewer_data.js").read_text(encoding="utf-8").removeprefix("window.BOLD_AB=").removesuffix(";"))
    assert payload["task"]=="MI" and len(payload["frame"])==110
    box = np.frombuffer(base64.b64decode(payload["box"]),"<f4").reshape(110,1215)
    default, = np.flatnonzero((PARAMETERS==[0,1]).all(axis=1))
    assert box[7,default]==1. and box[12,default]==0.
    assert w.dispatch(config,"glm")["glm_fits"]==0
    assert w.dispatch(config,"report")["action"]=="REUSED_IDENTICAL_REPORT"


def test_1215_is_opt_in_and_task_labels_are_distinct(dataset):
    config = w.read(dataset)
    config["analysis"].pop("double_scanner_ab")
    w.write(dataset,config)
    with pytest.raises(w.StageError,match="1024"):
        w.load(dataset)
    config["analysis"]["double_scanner_ab"] = dict(seed=1,device="cpu")
    config["analysis"]["groups"] = [dict(id="cross_task",runs=["run_1","run_2"],estimand_id="MI",same_estimand_confirmed=True)]
    w.write(dataset,config)
    with pytest.raises(w.StageError,match="cannot pool"):
        w.load(dataset)


def test_compare_blocks_saved_beta_on_shifted_grid(dataset, monkeypatch):
    import bold_stage

    config = w.load(dataset)
    source = dataset.parent / "saved"
    pair = config["analysis"]["comparisons"][0]
    pair["model"] = "GLM_TD"
    config["analysis"]["kernels_mm"] = [2]
    w.write(source / "summary.json", {"runs": [
        {"id": side, "status": "FITTED"} for side in (pair["left"], pair["right"])
    ], "groups": []})
    for offset, key in enumerate(("left", "right")):
        folder = source / pair[key]
        w.write(folder / "summary.json", {"reference_space": "MNI152NLin2009cAsym"})
        nib.save(nib.Nifti1Image(np.ones((2, 2, 2), dtype=np.uint8), np.eye(4)),
                 folder / "mask.nii.gz")
        (folder / "fwhm_2mm").mkdir()
        affine = np.eye(4)
        affine[0, 3] = offset
        nib.save(nib.Nifti1Image(np.ones((2, 2, 2), dtype=np.float32), affine),
                 folder / "fwhm_2mm/effect.nii.gz")
    monkeypatch.setattr(bold_stage, "saved_fit", lambda *_: ("fixture", source))
    output = Path(bold_stage.compare(config)["output"])
    result, = w.read(output / "summary.json")["comparisons"]
    assert result["status"] == "BLOCKED"
    assert "physical grid" in result["reason"]
    assert not (output / pair["id"]).exists()

    members = [pair["left"], pair["right"]]
    config["analysis"]["statistics"]["fdr_q"] = 0.05
    for member in members:
        folder = source / member
        w.write(folder / "summary.json", {
            "reference_space": "MNI152NLin2009cAsym",
            "kernels": {"2": {"residual_dof": 20}},
        })
        nib.save(nib.Nifti1Image(np.ones((2, 2, 2), dtype=np.float32), np.eye(4)),
                 folder / "fwhm_2mm/variance.nii.gz")
    runs = {member: {"id": member, "status": "FITTED", "task_id": "synthetic",
                     "timing": {"status": "engineering_subset"}} for member in members}
    group = {"id": "pooled", "runs": members, "estimand_id": "synthetic"}
    with pytest.raises(w.StageError, match="physical grid"):
        bold_stage.fixed_effect(config, group, runs, source)


def test_method_review_is_pinned_and_missing_only_blocks_ab(dataset):
    c = w.read(dataset)
    c["preprocessing"]["mode"] = "existing"
    c["analysis"]["timing"] = dict(status="confirmed",evidence="timing.json",interpretation="Synthetic proof fixture")
    c["analysis"]["runs"][0]["method_review"] = "method.json"
    w.write(dataset,c)
    config = w.load(dataset)
    project = w.source_root(config)
    w.write(project/"logs/preprocessing_status.json",{"exit_code":0})
    w.write(project/"logs/fmriprep_run_configuration.json",{"arguments":["--omp-nthreads","4"]})
    source_report = w.derivatives(config)/"sub-synthetic.html"
    source_report.write_text("<p>Synthetic review fixture</p>",encoding="utf-8")
    pre = w.dispatch(config,"preprocess")
    w.write(dataset.parent/"timing.json",dict(status="CONFIRMED_SCAN_SPECIFIC_TIMING",case_id="synthetic",
            preprocessing_fingerprint=pre["fingerprint"],
            event_sha256_by_run={r["id"]:w.sha(Path(r["events"])) for r in config["analysis"]["runs"]}))
    from bold_stage import inspect_runs
    _,rows = inspect_runs(config,w.require_preprocessing(config))
    assert rows[0]["status"]=="DESIGN_READY" and "ab_unavailable" in rows[0]["timing"]
    run = config["analysis"]["runs"][0]
    w.write(dataset.parent/"method.json",dict(status="RUN_METHOD_REVIEWED",case_id="synthetic",run_id=run["id"],
            bold_sha256=w.sha(w.derivatives(config)/run["bold"]),report="sub-synthetic.html",report_sha256=w.sha(source_report),
            reviewer="Synthetic fixture",registration_review="synthetic result",observed_sdc="off",observed_stc="off"))
    _,rows = inspect_runs(config,w.require_preprocessing(config))
    assert rows[0]["timing"]["methods"]["status"]=="RUN_METHOD_REVIEWED"
    assert "ab_unavailable" not in rows[0]["timing"] and "ab_unavailable" in rows[1]["timing"]
