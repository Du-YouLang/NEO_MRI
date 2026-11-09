
"""Independent GLS/bootstrap reference and end-to-end portable BOLD checks."""
import base64
import json
from pathlib import Path
import sys

import nibabel as nib
import numpy as np
import pandas as pd
import pytest
from nilearn.glm.first_level import compute_regressor
from scipy.stats import t

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"bold"))
import bold_workflow as w
from bold_ab import (LAGS, RATIOS, PARAMETERS, classify, corrected_family,
                     bootstrap_choices, chosen_gls, stability_frequencies)
from bold_models import MOTION, whiten


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
    project = tmp_path/"source"
    derivatives = project/"derivatives/fmriprep"
    func = derivatives/"sub-synthetic/func"
    func.mkdir(parents=True)
    n = 110
    rng = np.random.default_rng(71)
    shape = (4,4,4)
    affine = np.diag([2.,2.,2.,1.])
    mask = np.zeros(shape,np.uint8)
    mask[1:3,1:3,1:3] = 1
    onsets = np.array([14.,65.,122.,171.])
    durations = np.full(len(onsets),10.)
    reg,_ = compute_regressor(np.vstack((onsets+3.,durations*1.15,np.ones(len(onsets)))),
                              "spm",np.arange(n)*2,oversampling=50)
    w.write(project/"subset_provenance.json",{"kind":"SYNTHETIC_ENGINEERING_INPUTS","RepetitionTime":2.})
    t1 = derivatives/"synthetic_T1.nii.gz"
    t1_mask = derivatives/"synthetic_T1_mask.nii.gz"
    nib.save(nib.Nifti1Image(np.arange(64,dtype=np.float32).reshape(shape)+2,affine),t1)
    nib.save(nib.Nifti1Image(np.ones(shape,np.uint8),affine),t1_mask)
    runs = []
    for i,task in enumerate(("MI","PM"),1):
        stem = f"sub-synthetic_task-test_run-{i:02d}"
        prefix = stem+"_space-MNI152NLin2009cAsym_res-2"
        bold = func/(prefix+"_desc-preproc_bold.nii.gz")
        img = nib.Nifti1Image((100+reg[:,0]*4+rng.normal(size=shape+(n,))*.35).astype(np.float32),affine)
        img.header.set_zooms((2,2,2,2))
        img.header.set_xyzt_units("mm","sec")
        nib.save(img,bold)
        mask_path = func/(prefix+"_desc-brain_mask.nii.gz")
        nib.save(nib.Nifti1Image(mask,affine),mask_path)
        sidecar = func/(prefix+"_desc-preproc_bold.json")
        w.write(sidecar,dict(RepetitionTime=2.,SliceTimingCorrected=False))
        data = {}
        for motion in MOTION:
            a = rng.normal(size=n)
            diff = np.r_[np.nan,np.diff(a)]
            data.update({motion:a,motion+"_derivative1":diff,motion+"_power2":a*a,motion+"_derivative1_power2":diff*diff})
        data.update(framewise_displacement=np.r_[np.nan,np.zeros(n-1)],std_dvars=np.r_[np.nan,np.ones(n-1)])
        confounds = func/(stem+"_desc-confounds_timeseries.tsv")
        pd.DataFrame(data).to_csv(confounds,sep="\t",index=False)
        w.write(confounds.with_suffix(".json"),{})
        events = tmp_path/f"events_{i}.tsv"
        pd.DataFrame(dict(onset=onsets,duration=durations,trial_type=task)).to_csv(events,sep="\t",index=False)
        rel = lambda path:path.relative_to(derivatives).as_posix()
        runs.append(dict(id=f"run_{i}",task_id=task,estimand_id=task,reference_space="MNI152NLin2009cAsym",
                         tr_seconds=2.,stc="off",sdc="off",bold=rel(bold),mask=rel(mask_path),bold_json=rel(sidecar),
                         confounds=rel(confounds),confounds_json=rel(confounds.with_suffix(".json")),
                         viewer_anatomy=rel(t1),viewer_anatomy_mask=rel(t1_mask),events=events.name,
                         tr_evidence="source/subset_provenance.json",dynamic_condition=task,contrast={task:1}))
    config = dict(schema="neo-bold-stages-2",case_id="synthetic",data_root="data",workspace="data/synthetic/s01/results/task_fMRI/ab",
        preprocessing=dict(project="source",mode="engineering_subset"),
        analysis=dict(runs=runs,kernels_mm=[0,2,3],fir_delays_tr=[0,1,2,3],
                      timing=dict(status="engineering_subset",evidence="source/subset_provenance.json",interpretation="Synthetic fixture"),
                      double_scanner=dict(lags_seconds=list(LAGS),duration_ratios=list(RATIOS),selection="raw_r2"),
                      double_scanner_ab=dict(seed=31,device="cpu"),
                      statistics=dict(hrf_model="spm + derivative",smoothing_fwhm_mm=2,noise_model="ar1",compcor_components=0,
                                      fd_threshold_mm=.5,dvars_threshold=1.5,contrast={"MI":1},fdr_q=.05),
                      comparisons=[dict(id="mi_pm",left="run_1",right="run_2",model="DoubleScanner_AB")]))
    path = tmp_path/"config.json"
    w.write(path,config)
    return path


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
