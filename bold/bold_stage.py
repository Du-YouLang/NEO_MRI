
"""Run-wise BOLD fitting, explicit fixed effects and saved-beta comparison."""
from pathlib import Path
import time

import nibabel as nib
import numpy as np
import pandas as pd
from nilearn.image import smooth_img
from scipy.stats import t as student_t

import bold_workflow as w
from bold_models import bh, build_designs, correlation_rows, estimate_rho, fit_models, matrix_gate, whiten
from bold_runs import RUN_FILES, clock_for, contrast_for, file_for, run_identity, stem

StageError = w.StageError
TimingBlocked = w.TimingBlocked


def method_evidence(config,run,paths,registered):
    """Pin the operator's run-specific SDC and registration review to its actual inputs."""
    source = run.get("method_review")
    if not source or not w.resolve(source).is_file():
        return dict(status="NOT_REVIEWED",sdc=run.get("sdc"))
    record = w.read(w.resolve(source))
    base = w.derivatives(config).resolve()
    bold_sha = registered[paths["bold"].relative_to(base).as_posix()]
    if (record.get("status")!="RUN_METHOD_REVIEWED" or record.get("case_id")!=config["case_id"]
        or record.get("run_id")!=run["id"] or record.get("bold_sha256")!=bold_sha
        or record.get("observed_sdc")!=run.get("sdc") or not run.get("sdc")
        or record.get("observed_stc")!=run.get("stc")
        or not record.get("registration_review") or not record.get("reviewer")):
        raise StageError("Method review does not match the declared run, source BOLD, STC, SDC or registration")
    report = file_for(config,registered,{"id":run["id"],"report":record.get("report")},"report")
    if record.get("report_sha256")!=registered[report.relative_to(base).as_posix()]:
        raise StageError("Reviewed fMRIPrep report changed")
    return dict(status="RUN_METHOD_REVIEWED",sdc=run["sdc"],registration=record["registration_review"],
                fmriprep_report=str(report),fmriprep_report_sha256=record["report_sha256"])


def prepare_run(config, pre, run):
    if not run.get("task_id"):
        raise TimingBlocked("Confirm the scientific task for this run")
    timing = run.get("timing", config["analysis"]["timing"])
    if not {"status","evidence","interpretation"} <= set(timing):
        raise TimingBlocked("Run timing clearance is incomplete")
    paths = {key:file_for(config,pre,run,key) for key in RUN_FILES}
    prefix = stem(paths["bold"])
    expected = {
        "bold_json": paths["bold"].name.removesuffix(".nii.gz")+".json",
        "mask": paths["bold"].name.replace("_desc-preproc_bold.nii.gz","_desc-brain_mask.nii.gz"),
        "confounds": prefix+"_desc-confounds_timeseries.tsv",
        "confounds_json": prefix+"_desc-confounds_timeseries.json",
    }
    if not paths["bold"].name.endswith("_desc-preproc_bold.nii.gz") or any(
        paths[key].name!=name or paths[key].parent!=paths["bold"].parent for key,name in expected.items()
    ):
        raise StageError("BOLD, mask, sidecar and confounds do not identify the same run/grid")
    for key in ("events","tr_evidence"):
        name = run.get(key)
        if not name or not w.resolve(name).is_file():
            raise TimingBlocked(f"Missing {key} for this run")
        paths[key] = w.resolve(name)
    analysis = config["analysis"]
    run_case = {**config,"analysis":{**analysis,"runs":[run],"timing":timing,
                "statistics":{**analysis["statistics"],"contrast":contrast_for(config,run)}}}
    if (timing["status"]=="engineering_subset") != (config["preprocessing"]["mode"]=="engineering_subset"):
        raise TimingBlocked("Engineering timing is restricted to engineering input receipts")
    proof = w.task_gate(run_case)
    proof["methods"] = method_evidence(config,run,paths,pre)
    if config["analysis"].get("double_scanner_ab") and config["preprocessing"]["mode"]=="existing":
        if proof["methods"]["status"]!="RUN_METHOD_REVIEWED":
            proof["ab_unavailable"] = "Run-specific fMRIPrep SDC and registration review is missing"
    image = nib.load(paths["bold"])
    clock = clock_for(run,image,w.read(paths["bold_json"]),w.read(paths["tr_evidence"]).get("RepetitionTime"))
    events = pd.read_csv(paths["events"],sep="\t")
    designs = build_designs(config,run,clock,events,pd.read_csv(paths["confounds"],sep="\t"),
                            w.read(paths["confounds_json"]))
    return paths,image,clock,proof,designs


def inspect_runs(config, receipt):
    registered = {entry["path"]:entry["sha256"] for entry in receipt["files"]}
    ready = {}
    rows = []
    for run in config["analysis"]["runs"]:
        row = dict(id=run["id"],task_id=run.get("task_id"),estimand_id=run.get("estimand_id"),
                   contrast=contrast_for(config,run))
        try:
            ready[run["id"]] = prepare_run(config,registered,run)
            _,image,clock,proof,designs = ready[run["id"]]
            row.update(status="DESIGN_READY",frames=len(clock.times),tr_seconds=clock.tr,
                       start_seconds=clock.start,stc=clock.stc,timing=proof,
                       shape=list(image.shape),audit=designs.audit,unavailable=designs.unavailable)
        except TimingBlocked as exc:
            row.update(status="BLOCKED_METADATA_OR_TIMING",reason=str(exc))
        except (StageError,KeyError,ValueError,OSError,nib.filebasedimages.ImageFileError) as exc:
            row.update(status="BLOCKED_INPUT_OR_DESIGN",reason=str(exc))
        rows.append(row)
    return ready,rows


def audit(config, full=False):
    pre = w.require_preprocessing(config,full)
    key,_ = run_identity(config,pre)
    _,rows = inspect_runs(config,pre)
    dest = w.workspace(config)/"audits"/key/"model_audit.json"
    w.write(dest,dict(runs=rows))
    return dict(stage="audit",output=str(dest),ready_runs=sum(r["status"]=="DESIGN_READY" for r in rows),
                glm_fits=0,fmriprep_invocations=0)


def same_grid(a,b):
    return a.shape[:3]==b.shape[:3] and np.allclose(a.affine,b.affine,atol=1e-5,rtol=0)


def read_mask(path, reference):
    image = nib.load(path)
    if len(image.shape)!=3 or not same_grid(image,reference):
        raise StageError(f"Mask has a different physical grid: {path}")
    values = image.get_fdata(dtype=np.float32,caching="unchanged")
    if not np.isin(values,[0.,1.]).all():
        raise StageError(f"Expected a binary mask: {path}")
    return values.astype(bool)


def map_file(folder, name, vector, mask, reference, dtype=np.float32):
    outside = 1 if name in ("p_two_sided","q_bh","fir_joint_p","fir_joint_q_bh","ab_p_nominal","ab_q_A","ab_q_B_exploratory") else 0
    data = np.full(mask.shape,outside,dtype=dtype)
    data[mask] = vector
    image = nib.Nifti1Image(data,reference.affine)
    image.header.set_xyzt_units(reference.header.get_xyzt_units()[0])
    nib.save(image,folder/f"{name}.nii.gz")


def summarize(values):
    q = np.quantile(values,[.25,.5,.75])
    return dict(n=len(values),mean=float(np.mean(values)),q25=float(q[0]),median=float(q[1]),q75=float(q[2]))


def histogram_tables(metrics):
    result = {}
    for key in next(iter(metrics.values())):
        pool = np.concatenate([m[key] for m in metrics.values()])
        lo,hi = np.quantile(pool,[.01,.99])
        bounds = (float(pool.min()),float(pool.max())),(float(lo),float(hi))
        views = {}
        for view,(left,right) in zip(("full","central"),bounds):
            if left==right:
                padding = max(abs(left)*1e-4,1e-8)
                left,right = left-padding,right+padding
            edges = np.linspace(left,right,33)
            views[view] = dict(edges=edges.tolist(),series={
                name:dict(counts=np.histogram(values[key],edges)[0].tolist(),
                          below=int((values[key]<left).sum()),above=int((values[key]>right).sum()))
                for name,values in metrics.items()})
        result[key] = views
    return result


def save_static(folder, effect, variance, statistic, p, mask, image, fdr_q):
    q = bh(p)
    for name,values in (("effect",effect),("variance",variance),("t",statistic),("p_two_sided",p),("q_bh",q)):
        map_file(folder,name,values,mask,image)
    return dict(effect_psc=summarize(effect),q_threshold=fdr_q,voxels_fdr=int((q<=fdr_q).sum()))


def fit_run(config, pre, run, prepared, dest):
    paths,image,clock,proof,designs = prepared
    analysis = config["analysis"]
    stats = analysis["statistics"]
    ab_policy = analysis.get("double_scanner_ab")
    ab_block = proof.get("ab_unavailable")
    if ab_policy and (designs.dynamic is None or not np.array_equal(designs.dynamic_nuisance,designs.nuisance.to_numpy())):
        ab_block = "A/B requires exactly one estimable task condition and the same N0 for both gains"
    mask = read_mask(paths["mask"],image)
    if analysis.get("engine_variant") == "anatomical_intersection_v1":
        anatomy = file_for(config,{r["path"] for r in pre["files"]},
                           {"id":run["id"],"anatomical_mask":analysis.get("anatomical_mask")},"anatomical_mask")
        mask &= read_mask(anatomy,image)
    source = image.get_fdata(dtype=np.float32,caching="unchanged")
    if not np.isfinite(source).all():
        raise StageError("Nonfinite BOLD; review data before smoothing or fitting")
    inside = source[mask,:]
    valid = (inside.mean(axis=1)>0) & (inside.std(axis=1)>1e-6)
    del inside
    mask[mask] = valid
    if not mask.any():
        raise StageError("No finite, nonconstant BOLD voxel in the declared mask")
    raw = source[mask,:].T.copy()
    raw = 100*(raw/raw.mean(axis=0)-1)
    rho = estimate_rho(raw,designs.static.to_numpy(),stats["noise_model"])
    if 0 not in analysis["kernels_mm"]:
        del raw
    typical_rho = float(rho[np.argmin(abs(rho-np.median(rho)))])
    white_static = whiten(designs.static.to_numpy(),typical_rho)
    white_nuisance = whiten(designs.nuisance.to_numpy(),typical_rho)
    task_count = white_static.shape[1]-white_nuisance.shape[1]
    white_task = white_static[:,:task_count]
    white_audit = dict(representative_rho=typical_rho,voxel_count_at_rho=int((rho==typical_rho).sum()),
                       **matrix_gate(white_static),
                       design_vs_confound=correlation_rows(white_task,list(designs.static.columns[:task_count]),
                                                         white_nuisance,list(designs.nuisance.columns)))
    dest.mkdir()
    map_file(dest,"mask",1,mask,image,np.uint8)
    pd.DataFrame({"seconds":clock.times}).to_csv(dest/"sample_clock.tsv",sep="\t",index=False)
    designs.static.to_csv(dest/"glm_td_design.tsv",sep="\t",index=False)
    if designs.fir is not None:
        designs.fir.to_csv(dest/"fir_design.tsv",sep="\t",index=False)
    np.savez(dest/"design_arrays.npz",nuisance=designs.nuisance.to_numpy(),
             **({"dynamic":designs.dynamic,"dynamic_nuisance":designs.dynamic_nuisance}
                if designs.dynamic is not None else {}))
    summary = dict(ab_unavailable=ab_block,voxels=int(mask.sum()),discarded_voxels=int((~valid).sum()),primary_kernel_mm=stats["smoothing_fwhm_mm"],
                   dynamic_condition=run.get("dynamic_condition"),dynamic_grid=designs.parameters,
                   dynamic_selection=analysis["double_scanner"]["selection"],rho=summarize(rho),kernels={},
                   static_contrast=contrast_for(config,run),fir_delays_tr=analysis["fir_delays_tr"],white_domain_example=white_audit,
                   reference_space=run.get("reference_space"))
    for kernel in analysis["kernels_mm"]:
        folder = dest/f"fwhm_{kernel}mm"
        folder.mkdir()
        if kernel==0:
            y = raw
        else:
            # One spatial smoothing call per kernel, starting from the same run.
            # Zero outside the declared source support is normalized away.
            support = nib.Nifti1Image(mask.astype(np.float32),image.affine)
            weight = np.asarray(smooth_img(support,kernel).dataobj)[mask]
            masked = nib.Nifti1Image(np.where(mask[...,None],source,0.),image.affine)
            y = np.asarray(smooth_img(masked,kernel).dataobj)[mask,:].T / weight
            mean = y.mean(axis=0)
            if (mean<=0).any() or (y.std(axis=0)<=1e-6).any():
                raise StageError("Kernel changes valid voxel support; review input before cross-kernel comparison")
            y = 100*(y/mean-1)
        ab_enabled = bool(ab_policy and not ab_block and kernel==stats["smoothing_fwhm_mm"])
        default = next(i for i,pair in enumerate(designs.parameters) if pair==[0.,1.]) if ab_enabled else None
        fit = fit_models(y,designs,rho,summary["dynamic_selection"],default_choice=default)
        entry = dict(residual_dof=fit["dof"],models={name:{k:summarize(v) for k,v in metrics.items()}
                    for name,metrics in fit["metrics"].items()})
        entry["GLM_TD"] = save_static(folder,fit["effect"],fit["variance"],fit["t"],fit["p_two_sided"],
                                     mask,image,stats["fdr_q"])
        np.savez(folder/"coefficients.npz",**fit["coefficients"],rho=rho,
                 **({"ds_choice":fit["best"]} if fit["best"] is not None else {}))
        if fit["best"] is not None:
            lag,ratio = np.asarray(designs.parameters)[fit["best"]].T
            for name,values in (("ds_beta",fit["coefficients"]["DoubleScanner"][0]),("ds_lag_seconds",lag),("ds_duration_ratio",ratio)):
                map_file(folder,name,values,mask,image)
        if fit["fir_p"] is not None:
            fir_q = bh(fit["fir_p"])
            for name,values in (("fir_joint_f",fit["fir_f"]),("fir_joint_p",fit["fir_p"]),("fir_joint_q_bh",fir_q)):
                map_file(folder,name,values,mask,image)
            entry["FIR"] = dict(joint_fdr_voxels=int((fir_q<=stats["fdr_q"]).sum()),
                                numerator_dof=len(designs.fir_contrasts),denominator_dof=len(y)-designs.fir.shape[1])

        # One real recorded time course per kernel. The same representative voxel
        # is used across kernels and chosen independently of fitted statistics.
        coordinates = np.argwhere(mask)
        j = int(np.argmin(((coordinates-coordinates.mean(axis=0))**2).sum(axis=1)))
        coef = fit["coefficients"]
        trace = dict(seconds=clock.times.tolist(),voxel=coordinates[j].tolist(),observed=y[:,j].tolist(),
                     N0=(designs.nuisance.to_numpy()@coef["N0"][:,j]).tolist(),
                     GLM_TD=(designs.static.to_numpy()@coef["GLM_TD"][:,j]).tolist())
        if fit["best"] is not None:
            choice = int(fit["best"][j])
            template = designs.dynamic[:,choice]
            lag,ratio = designs.parameters[choice]
            trace.update(DoubleScanner=(designs.dynamic_nuisance@coef["DoubleScanner"][1:,j]+
                                        template*coef["DoubleScanner"][0,j]).tolist(),
                         dynamic_post_hrf=(template*coef["DoubleScanner"][0,j]).tolist(),
                         dynamic_event_windows=[[float(on+lag),float(duration*ratio)]
                                                for on,duration,_ in designs.target_events.T],
                         selected_lag_seconds=lag,selected_duration_ratio=ratio)
        if "FIR" in coef:
            trace["FIR"] = (designs.fir.to_numpy()@coef["FIR"][:,j]).tolist()
        if ab_enabled:
            from bold_ab import fit_ab
            ab = fit_ab(y,designs,rho,fit,ab_policy)
            entry["DoubleScanner_AB"] = save_ab(folder,mask,image,ab,designs,clock,ab_policy)
        w.write(folder/"trace.json",trace)
        w.write(folder/"histograms.json",histogram_tables(fit["metrics"]))
        summary["kernels"][str(kernel)] = entry
    w.write(dest/"summary.json",summary)


def save_ab(folder,mask,image,result,designs,clock,policy):
    from bold_ab import ALPHA, FACTORS, REASONS, TASK_GAIN, DEFAULT_GAIN, STABILITY
    np.savez_compressed(folder/"double_scanner_ab.npz",**result,
                        target_events=designs.target_events,templates=designs.dynamic,clock=clock.times)
    scan = result["reason"]==0
    for name,key,dtype in (("ab_final_beta","beta",np.float32),("ab_final_t","t",np.float32),
                          ("ab_p_nominal","p",np.float32),("ab_q_A","q_A",np.float32),
                          ("ab_q_B_exploratory","q_B",np.float32),("ab_A","A",np.uint8),
                          ("ab_B_exploratory","B",np.uint8),("ab_A_stable","A_stable",np.uint8),
                          ("ab_B_stable_exploratory","B_stable",np.uint8)):
        map_file(folder,name,result[key],mask,image,dtype)
    map_file(folder,"ab_model_source",result["reason"]+1,mask,image,np.uint8)
    for name,axis in (("ab_lag_seconds",0),("ab_duration_ratio",1)):
        map_file(folder,name,result["parameters"][result["final"],axis],mask,image)
    counts = {name:dict(n=int(result[name].sum()),scan=int((result[name]&scan).sum()),
                        fallback=int((result[name]&~scan).sum())) for name in ("A","B","A_stable","B_stable")}
    return dict(voxels=len(scan),candidates=len(result["parameters"]),factors=FACTORS,alpha=ALPHA,df=int(result["df"]),
                counts=counts,fallback=dict(zip(REASONS,np.bincount(result["reason"],minlength=6).tolist())),
                changed_at_1215=int(np.sum(result["A"]!=(result["q_1215"]<=ALPHA))),
                thresholds=dict(raw_relative_n0=TASK_GAIN,raw_gain_over_default=DEFAULT_GAIN,stable=STABILITY),
                bootstrap_draws=len(result["draws"]),bootstrap_seed=policy["seed"],bootstrap_device=policy["device"],
                B_exploratory=True,conditional_on_estimated_ar1=True)


def compatible_masks(paths):
    images = [nib.load(path) for path in paths]
    if not all(same_grid(images[0],image) for image in images):
        raise StageError("Physical grids differ; select registered outputs on a common grid")
    mask = np.logical_and.reduce([np.asanyarray(image.dataobj)==1 for image in images])
    if not mask.any():
        raise StageError("No voxels in the intersection of the saved masks")
    return images[0],mask


def load_stat_map(path, reference):
    image = nib.load(path)
    if len(image.shape) != 3 or not same_grid(image, reference):
        raise StageError(f"Saved statistical map has a different physical grid: {path}")
    return image


def fixed_effect(config, group, runs, dest):
    members = [runs[member] for member in group["runs"]]
    if any(row["status"]!="FITTED" for row in members):
        raise StageError("All group members must have valid individual fits")
    # Conditional/exploratory runs can still be inspected and compared as saved
    # maps. Do not promote them to a joint scientific inference.
    if any(row["timing"]["status"] not in ("confirmed","engineering_subset") for row in members):
        raise StageError("Fixed effects require confirmed timing for every scientific run")
    folders = [dest/member["id"] for member in members]
    contexts = [w.read(folder/"summary.json") for folder in folders]
    if not contexts[0]["reference_space"] or any(c["reference_space"]!=contexts[0]["reference_space"] for c in contexts):
        raise StageError("Fixed effects require the same explicitly declared reference space")
    image,mask = compatible_masks([folder/"mask.nii.gz" for folder in folders])
    out = dest/group["id"]
    out.mkdir()
    map_file(out,"mask",1,mask,image,np.uint8)
    kernels = {}
    for kernel in config["analysis"]["kernels_mm"]:
        name = f"fwhm_{kernel}mm"
        variances = np.stack([load_stat_map(f/name/"variance.nii.gz", image).get_fdata()[mask] for f in folders])
        effects = np.stack([load_stat_map(f/name/"effect.nii.gz", image).get_fdata()[mask] for f in folders])
        dofs = np.asarray([s["kernels"][str(kernel)]["residual_dof"] for s in contexts],dtype=float)
        if (variances<=0).any() or not np.isfinite(variances).all():
            raise StageError("Invalid saved variance in fixed-effects domain")
        precision = 1/variances
        total = precision.sum(axis=0)
        variance = 1/total
        effect = (effects*precision).sum(axis=0)/total
        # Satterthwaite approximation for independent estimated run variances.
        dof = total**2 / ((precision**2)/dofs[:,None]).sum(axis=0)
        t = effect/np.sqrt(variance)
        folder = out/name
        folder.mkdir()
        kernels[str(kernel)] = save_static(folder,effect,variance,t,2*student_t.sf(abs(t),dof),mask,image,
                                          config["analysis"]["statistics"]["fdr_q"])
        map_file(folder,"approximate_dof",dof,mask,image)
    summary = dict(status="FITTED_FIXED_EFFECT_GLM_TD",id=group["id"],members=group["runs"],
                   task_id=members[0]["task_id"],estimand_id=group["estimand_id"],kernels=kernels,
                   voxels=int(mask.sum()),reference_space=contexts[0]["reference_space"])
    w.write(out/"summary.json",summary)
    return summary


def glm(config, full=False):
    start = time.perf_counter()
    pre = w.require_preprocessing(config,full)
    key,spec = run_identity(config,pre)
    dest = w.workspace(config)/"glm"/key
    if w.completed(dest,full):
        return dict(stage="glm",action="REUSED_IDENTICAL_STATISTICS",key=key,output=str(dest),glm_fits=0,fmriprep_invocations=0)
    if dest.exists():
        raise StageError(f"Incomplete statistical directory preserved: {dest}")
    ready,rows = inspect_runs(config,pre)
    if not ready:
        review = w.workspace(config)/"audits"/key/"model_audit.json"
        w.write(review,dict(runs=rows))
        raise TimingBlocked(f"No run cleared for fitting. See {review}")
    dest.mkdir(parents=True)
    w.write(dest/"request.json",spec)
    runs = {run["id"]:run for run in config["analysis"]["runs"]}
    groups = []
    try:
        for row in rows:
            name = row["id"]
            if name in ready:
                try:
                    fit_run(config,pre,runs[name],ready[name],dest/name)
                    row["status"] = "FITTED"
                except StageError as exc:
                    row.update(status="BLOCKED_FIT",reason=str(exc))
        fitted_count = sum(row["status"]=="FITTED" for row in rows)
        if not fitted_count:
            raise StageError(f"No run completed voxel fitting. See {dest/'failure.json'}")
        for group in config["analysis"].get("groups",[]):
            try:
                groups.append(fixed_effect(config,group,{r["id"]:r for r in rows},dest))
            except StageError as exc:
                groups.append(dict(id=group["id"],status="BLOCKED",reason=str(exc),members=group["runs"]))
        if run_identity(config,w.require_preprocessing(config,False))[0]!=key:
            raise StageError("Analysis inputs changed during fitting; this attempt cannot be registered")
        w.write(dest/"summary.json",dict(runs=rows,groups=groups,primary_kernel_mm=config["analysis"]["statistics"]["smoothing_fwhm_mm"]))
        files = [w.descriptor(path,dest) for path in sorted(dest.rglob("*")) if path.is_file()]
        w.write(dest/"receipt.json",dict(status="COMPLETED_NUMERICALLY_REVIEW_PENDING",key=key,files=files,
                preprocessing_fingerprint=pre["fingerprint"],finished_UTC=w.now(),clinical_validation=False))
    except Exception as exc:
        w.write(dest/"failure.json",dict(error=str(exc),runs=rows,UTC=w.now()))
        raise
    return dict(stage="glm",action="FITTED_NEW_STATISTICAL_VERSION",key=key,output=str(dest),glm_fits=fitted_count,
                blocked_runs=len(rows)-fitted_count,fmriprep_invocations=0,elapsed_seconds=time.perf_counter()-start)


def saved_fit(config, full=False):
    pre = w.require_preprocessing(config,full)
    key,_ = run_identity(config,pre)
    dest = w.workspace(config)/"glm"/key
    if not w.completed(dest,full):
        raise StageError("Current fit is unavailable; this stage only reads saved results")
    return key,dest


def comparison_location(config, fit_key):
    request = dict(fit_key=fit_key,pairs=config["analysis"].get("comparisons",[]),
                   kernels=config["analysis"]["kernels_mm"],engine_sha256=w.sha(Path(__file__)))
    return w.workspace(config)/"comparisons"/w.digest(request),request


def compare(config, full=False):
    fit_key,source = saved_fit(config,full)
    dest,request = comparison_location(config,fit_key)
    if not request["pairs"]:
        raise StageError("No post-fit comparisons specified")
    if w.completed(dest,full):
        return dict(stage="compare",action="REUSED_COMPARISONS",output=str(dest),glm_fits=0,fmriprep_invocations=0)
    if dest.exists():
        raise StageError(f"Incomplete comparison directory preserved: {dest}")
    summary = w.read(source/"summary.json")
    eligible = {r["id"] for r in summary["runs"]+summary["groups"] if r["status"] in ("FITTED","FITTED_FIXED_EFFECT_GLM_TD")}
    dest.mkdir(parents=True)
    rows = []
    for pair in request["pairs"]:
        row = {**pair,"status":"DESCRIPTIVE_LEFT_MINUS_RIGHT"}
        try:
            if not {pair["left"],pair["right"]} <= eligible:
                raise StageError("Both sides need completed saved fits")
            sides = [source/pair[key] for key in ("left","right")]
            contexts = [w.read(f/"summary.json") for f in sides]
            if not contexts[0]["reference_space"] or contexts[0]["reference_space"]!=contexts[1]["reference_space"]:
                raise StageError("Explicit matching reference spaces required")
            image,mask = compatible_masks([f/"mask.nii.gz" for f in sides])
            model = pair.get("model","GLM_TD")
            beta = {"GLM_TD":"effect","DoubleScanner":"ds_beta","DoubleScanner_AB":"ab_final_beta"}[model]
            if model=="DoubleScanner_AB":
                primary = config["analysis"]["statistics"]["smoothing_fwhm_mm"]
                kernels = [primary]
            else:
                kernels = request["kernels"]
            paths = [[side/f"fwhm_{kernel}mm"/f"{beta}.nii.gz" for side in sides] for kernel in kernels]
            if any(not p.is_file() for pair_paths in paths for p in pair_paths):
                raise StageError(f"Both sides need saved {model} beta maps for all kernels")
            beta_images = [[load_stat_map(path, image) for path in pair_paths] for pair_paths in paths]
            folder = dest/pair["id"]
            folder.mkdir()
            map_file(folder,"mask",1,mask,image,np.uint8)
            row.update(model=model,voxels=int(mask.sum()),kernels={},reference_space=contexts[0]["reference_space"])
            for kernel,pair_images in zip(kernels,beta_images,strict=True):
                where = folder/f"fwhm_{kernel}mm"
                where.mkdir()
                left,right = [beta_image.get_fdata()[mask] for beta_image in pair_images]
                delta = left-right
                for name,values in (("left",left),("right",right),("left_minus_right",delta)):
                    map_file(where,name,values,mask,image)
                row["kernels"][str(kernel)] = dict(left=summarize(left),right=summarize(right),difference=summarize(delta))
        except StageError as exc:
            row.update(status="BLOCKED",reason=str(exc))
        rows.append(row)
    w.write(dest/"summary.json",dict(comparisons=rows))
    files = [w.descriptor(path,dest) for path in sorted(dest.rglob("*")) if path.is_file()]
    w.write(dest/"receipt.json",dict(status="COMPLETED_NUMERICALLY_REVIEW_PENDING",files=files,request=request,
                                    finished_UTC=w.now(),clinical_validation=False))
    return dict(stage="compare",action="COMPARED_SAVED_BETA",output=str(dest),glm_fits=0,fmriprep_invocations=0)


def report(config,full=False):
    from bold_report import render
    fit_key,fit = saved_fit(config,full)
    comparison,request = comparison_location(config,fit_key)
    compared = w.completed(comparison,full) if request["pairs"] else None
    request = dict(fit_key=fit_key,comparison_key=comparison.name if compared else None,
                   pending_comparisons=request["pairs"] if not compared else [],context=config.get("report",{}),
                   renderer_sha256={name:w.sha(w.CODE/name) for name in
                                    ("bold_report.py","bold_ab_viewer.py","ab_viewer.html","ab_viewer.js","ab_viewer.css")})
    dest = w.workspace(config)/"reports"/w.digest(request)
    receipt = dest/"report_receipt.json"
    if receipt.is_file():
        w.check_files(dest,w.read(receipt)["files"],full)
        return dict(stage="report",action="REUSED_IDENTICAL_REPORT",output=str(dest/"report.html"),glm_fits=0,fmriprep_invocations=0)
    if dest.exists():
        raise StageError(f"Incomplete report directory preserved: {dest}")
    dest.mkdir(parents=True)
    render(config,fit,comparison if compared else None,request["pending_comparisons"],dest)
    w.write(receipt,dict(status="REPORT_GENERATED_REVIEW_PENDING",request=request,finished_UTC=w.now(),
                        files=[w.descriptor(path,dest) for path in sorted(dest.rglob("*")) if path.is_file()]))
    return dict(stage="report",action="RENDERED_FROM_SAVED_OUTPUTS",output=str(dest/"report.html"),glm_fits=0,fmriprep_invocations=0)
