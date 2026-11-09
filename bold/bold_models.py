"""Same-run GLM+TD, DoubleScanner and FIR with a shared validated sampling clock.

The DoubleScanner core uses the audited FWL equations, with explicit selection
and bounded voxel batches. All fitted values use percent signal change units.
"""
from dataclasses import dataclass

import numpy as np
import pandas as pd
from nilearn.glm.first_level import compute_regressor
from scipy.stats import t as student_t, f as fisher_f

MOTION=['trans_x','trans_y','trans_z','rot_x','rot_y','rot_z']
MOTION24=MOTION+[x+'_derivative1' for x in MOTION]+[x+'_power2' for x in MOTION]+[x+'_derivative1_power2' for x in MOTION]
from bold_workflow import StageError
from bold_runs import contrast_for

BATCH_VOXELS = 2048


@dataclass
class Designs:
    nuisance: pd.DataFrame
    static: pd.DataFrame
    static_contrast: np.ndarray
    fir: pd.DataFrame | None
    fir_contrasts: np.ndarray | None
    dynamic: np.ndarray | None
    dynamic_nuisance: np.ndarray | None
    parameters: list
    target_events: np.ndarray | None
    unavailable: dict
    audit: dict


def nuisance_matrix(table, meta, stats):
    compcor = [name for name in table if name.startswith("a_comp_cor_")
               and meta.get(name, {}).get("Mask") == "combined"
               and meta.get(name, {}).get("Retained", True)]
    ncomp = stats["compcor_components"]
    if len(compcor) < ncomp:
        raise StageError("Insufficient retained combined aCompCor components")
    columns = MOTION24 + compcor[:ncomp] + [name for name in table if name.startswith("cosine")]
    values = table[columns].astype(float).copy()
    for name in columns:
        missing = values[name].isna()
        if missing.any() and not (missing.sum() == 1 and missing.iloc[0] and "_derivative1" in name):
            raise StageError(f"Unexpected confound NaN: {name}")
    values = values.fillna(0)
    sd = values.std(ddof=0)
    if not np.isfinite(values.to_numpy()).all() or (sd <= 1e-12).any():
        raise StageError("Nonfinite or constant mandatory confound")
    values = (values-values.mean())/sd
    spikes = np.zeros(len(table), dtype=bool)
    for column, cutoff in (("framewise_displacement",stats["fd_threshold_mm"]),("std_dvars",stats["dvars_threshold"])):
        metric = table[column].astype(float).to_numpy()
        if np.isinf(metric).any() or np.isnan(metric[1:]).any():
            raise StageError(f"Unexpected missing/nonfinite motion metric: {column}")
        spikes |= np.nan_to_num(metric,nan=0.) > cutoff
    outliers = [name for name in table if name.startswith("non_steady_state_outlier")]
    if outliers:
        flags = table[outliers].to_numpy(dtype=float)
        if not np.isin(flags,[0.,1.]).all():
            raise StageError("Non-steady-state flags must be binary")
        spikes |= flags.any(axis=1)
    for index in np.flatnonzero(spikes):
        values[f"spike_{index:03d}"] = (np.arange(len(table)) == index).astype(float)
    values["constant"] = 1.
    return values, np.flatnonzero(spikes), compcor[:ncomp]


def matrix_gate(matrix):
    x = np.asarray(matrix,dtype=float)
    if x.ndim != 2 or not np.isfinite(x).all():
        raise StageError("Design must be a finite matrix")
    singular = np.linalg.svd(x,compute_uv=False)
    rank = int(np.linalg.matrix_rank(x))
    if rank != x.shape[1] or len(x)-rank < 10 or singular[0]/singular[-1] >= 1e8:
        raise StageError("Design is rank-deficient, ill-conditioned, or has fewer than ten residual degrees of freedom")
    return dict(rank=rank, residual_dof=len(x)-rank, condition_number=float(singular[0]/singular[-1]))


def event_array(group):
    amplitude = group["modulation"].to_numpy(float) if "modulation" in group else np.ones(len(group))
    return np.vstack((group["onset"].to_numpy(float),group["duration"].to_numpy(float),amplitude))


def regressor(events,times,model,name,delays=None):
    kwargs = {"fir_delays":delays} if delays is not None else {}
    values, labels = compute_regressor(events,model,times,con_id=name,oversampling=50,**kwargs)
    return pd.DataFrame(values,columns=labels)


def correlation_rows(matrix, labels, nuisance, nuisance_names):
    centered = matrix-matrix.mean(axis=0)
    conf = nuisance-nuisance.mean(axis=0)
    ssm = np.sum(centered**2,axis=0)
    ssc = np.sum(conf**2,axis=0)
    denom = np.sqrt(ssm[:,None]*ssc[None,:])
    corr = np.full(denom.shape,np.nan)
    np.divide(centered.T@conf,denom,out=corr,where=denom>1e-20)
    # Centered projection gives Pearson-compatible joint R². This also
    # handles the whitened intercept, whose first row has a distinct weight.
    resid = centered - conf @ (np.linalg.pinv(conf) @ centered)
    joint = np.full(len(labels),np.nan)
    np.divide(np.sum(resid**2,axis=0),ssm,out=joint,where=ssm>1e-20)
    return [dict(column=name,joint_nuisance_r2=float(np.clip(1-joint[j],0,1)) if np.isfinite(joint[j]) else None,
                 pearson={key:float(corr[j,k]) if np.isfinite(corr[j,k]) else None for k,key in enumerate(nuisance_names)})
            for j,name in enumerate(labels)]


def build_designs(config,run,clock,events,confounds,confounds_meta):
    if len(confounds) != len(clock.times):
        raise StageError("BOLD and confounds frame counts differ")
    needed = {"trial_type","onset","duration"}
    if not needed <= set(events) or events.empty or events.trial_type.isna().any():
        raise StageError("Events require nonempty trial_type, onset and duration")
    events = events.copy()
    events["trial_type"] = events["trial_type"].astype(str)
    values = events[["onset","duration"]].to_numpy(float)
    if not np.isfinite(values).all() or (values[:,1]<0).any():
        raise StageError("Events have invalid onset or duration")
    if (values[:,0]<clock.times[0]-24).any() or (values.sum(axis=1)>len(clock.times)*clock.tr+1e-4).any():
        raise StageError("Events extend outside the supported recorded time; no implicit trimming")
    nuisance, spikes, compcor = nuisance_matrix(confounds,confounds_meta,config["analysis"]["statistics"])
    matrix_gate(nuisance)
    conditions = {name:event_array(rows) for name,rows in events.groupby("trial_type",sort=True)}
    if any(not np.isfinite(e).all() for e in conditions.values()):
        raise StageError("Nonfinite event modulation")
    regressors = {name:regressor(e,clock.times,"spm + derivative",name) for name,e in conditions.items()}
    task = pd.concat(regressors.values(),axis=1)
    static = pd.concat([task,nuisance],axis=1)
    if not static.columns.is_unique:
        raise StageError("Duplicate design columns")
    diagnostic = matrix_gate(static)
    contrast = contrast_for(config,run)
    if not set(contrast) <= set(conditions):
        raise StageError("Static contrast must reference canonical event condition names")
    c = np.array([contrast.get(name,0.) for name in static],dtype=float)
    sampled = []
    for e in conditions.values():
        # Record event coverage on the saved volume interval, including impulses.
        onset,duration,_ = e
        origin = clock.times[:,None]-clock.start
        sampled.append(((origin < onset+np.maximum(duration,1e-12)) &
                        (origin+clock.tr > onset)).any(axis=1).astype(float))
    flags = np.column_stack(sampled)
    z = nuisance.to_numpy()
    x = static.to_numpy()
    vif = np.sum(np.linalg.pinv(x)**2,axis=1)*np.sum((x-x.mean(axis=0))**2,axis=0)
    audit = dict(**diagnostic,combined_compcor=compcor,spike_indices=spikes.tolist(),
                 spike_event_overlap={name:int(flags[spikes,i].sum()) for i,name in enumerate(conditions)},
                 event_vs_confound=correlation_rows(flags,list(conditions),z,list(nuisance)),
                 design_vs_confound=correlation_rows(task.to_numpy(),list(task),z,list(nuisance)),
                 VIF={name:float(vif[i]) if name!="constant" else None for i,name in enumerate(static)},
                 white_domain_diagnostic="Representative AR(1) design after fitting; the whitening transform can vary by voxel")

    fir = fir_contrasts = None
    dynamic = dynamic_nuisance = target_events = None
    params = []
    unavailable = {}
    delays = config["analysis"]["fir_delays_tr"]
    try:
        parts = {name:regressor(e,clock.times,"fir",name,delays) for name,e in conditions.items()}
        candidate = pd.concat([*parts.values(),nuisance],axis=1)
        if not candidate.columns.is_unique:
            raise StageError("Duplicate FIR columns")
        matrix_gate(candidate)
        weights = np.zeros((len(delays),candidate.shape[1]))
        for name,weight in contrast.items():
            for j,label in enumerate(parts[name].columns):
                weights[j,candidate.columns.get_loc(label)] = weight
        if np.linalg.matrix_rank(weights) != len(weights):
            raise StageError("FIR joint contrast is not full rank")
        fir,fir_contrasts = candidate,weights
    except StageError as exc:
        unavailable["FIR"] = str(exc)

    name = run.get("dynamic_condition")
    if name not in conditions:
        unavailable["DoubleScanner"] = "Select one verified dynamic event condition for this run"
    elif np.any(conditions[name][1]<=0):
        unavailable["DoubleScanner"] = "Duration-ratio search requires positive target-event durations"
    else:
        target_events = conditions[name]
        dynamic_nuisance = pd.concat([*(r for k,r in regressors.items() if k!=name),nuisance],axis=1).to_numpy()
        spec = config["analysis"]["double_scanner"]
        columns = []
        try:
            matrix_gate(dynamic_nuisance)
            for ratio in spec["duration_ratios"]:
                for lag in spec["lags_seconds"]:
                    e = target_events.copy()
                    e[0] += lag
                    e[1] *= ratio
                    if (e[0]<clock.times[0]-24).any():
                        raise StageError("DoubleScanner grid exceeds the supported pre-scan HRF window")
                    g = regressor(e,clock.times,"spm",name).to_numpy()[:,0]
                    matrix_gate(np.column_stack((dynamic_nuisance,g)))
                    columns.append(g)
                    params.append([float(lag),float(ratio)])
            dynamic = np.column_stack(columns)
        except StageError as exc:
            unavailable["DoubleScanner"] = str(exc)
            params = []
    return Designs(nuisance,static,c,fir,fir_contrasts,dynamic,dynamic_nuisance,params,target_events,unavailable,audit)


def whiten(x,rho):
    result = np.array(x,dtype=float,copy=True)
    result[1:] -= rho*np.asarray(x[:-1])
    return result


def estimate_rho(y,static,noise_model):
    if noise_model == "ols":
        return np.zeros(y.shape[1])
    pinv = np.linalg.pinv(static)
    rho = np.empty(y.shape[1])
    for start in range(0,len(rho),BATCH_VOXELS):
        where = slice(start,start+BATCH_VOXELS)
        values = y[:,where]
        error = values-static@(pinv@values)
        rho[where] = np.round(np.clip(np.sum(error[1:]*error[:-1],axis=0) /
                                     np.maximum(np.sum(error[:-1]**2,axis=0),1e-20),-.95,.95),2)
    return rho


def squared(x):
    return np.einsum("ij,ij->j",x,x)


def fit_models(y, designs, rho, selection, default_choice=None):
    y = np.asarray(y)
    rho = np.asarray(rho)
    if (y.ndim!=2 or rho.shape!=(y.shape[1],) or len(y)!=len(designs.static)
        or not np.isfinite(y).all() or not np.isfinite(rho).all() or (np.abs(rho)>=1).any()):
        raise StageError("BOLD, design and finite AR(1) coefficients must match frames and voxels")
    z = designs.nuisance.to_numpy()
    static = designs.static.to_numpy()
    fir = designs.fir.to_numpy() if designs.fir is not None else None
    g = designs.dynamic
    voxels = y.shape[1]
    models = {"N0":z.shape[1],"GLM_TD":static.shape[1]}
    if fir is not None:
        models["FIR"] = fir.shape[1]
    if g is not None:
        models["DoubleScanner"] = designs.dynamic_nuisance.shape[1]+1
    stats = {name:{key:np.full(voxels,np.nan) for key in ("r2","rmse","partial_r2_white","sse_reduction_raw","white_residual_acf1")}
             for name in models}
    coefs = {name:np.empty((n,voxels),np.float32) for name,n in models.items()}
    effect = np.empty(voxels)
    variance = np.empty(voxels)
    fir_f = np.empty(voxels) if fir is not None else None
    best = np.empty(voxels,dtype=np.int32) if g is not None else None
    if default_choice is not None:
        if g is None or selection!="raw_r2" or not 0<=default_choice<g.shape[1]:
            raise StageError("A grid default requires a valid raw-domain candidate")
        default_r2 = np.empty(voxels)
        raw_gain = np.empty(voxels)
    else:
        default_r2 = raw_gain = None
    raw = {"GLM_TD":static}
    if fir is not None:
        raw["FIR"] = fir

    for r in np.unique(rho):
        indices = np.flatnonzero(rho==r)
        wz = whiten(z,r)
        pz = np.linalg.pinv(wz)
        operators = {}
        for name,x in raw.items():
            wx = whiten(x,r)
            operators[name] = wx,np.linalg.pinv(wx)
        if g is not None:
            other = designs.dynamic_nuisance
            wo = whiten(other,r)
            po = np.linalg.pinv(wo)
            wg = whiten(g,r)
            residual_g = wg-wo@(po@wg)
            den = squared(residual_g)
            if (den<=1e-12).any():
                raise StageError("DoubleScanner template is confounded in the white domain")
            g_raw = g-other@(po@wg)
        for start in range(0,len(indices),BATCH_VOXELS):
            idx = indices[start:start+BATCH_VOXELS]
            yy = y[:,idx]
            yw = whiten(yy,r)
            sst = squared(yy-yy.mean(axis=0))
            coef0 = pz@yw
            fit0 = z@coef0
            raw0 = squared(yy-fit0)
            white0 = squared(yw-wz@coef0)

            def record(name,fit,white_fit,coef):
                residual = yy-fit
                error_w = yw-white_fit
                ss = squared(residual)
                sw = squared(error_w)
                coefs[name][:,idx] = coef
                stats[name]["r2"][idx] = 1-ss/np.maximum(sst,1e-20)
                stats[name]["rmse"][idx] = np.sqrt(ss/len(y))
                stats[name]["partial_r2_white"][idx] = 1-sw/np.maximum(white0,1e-20)
                stats[name]["sse_reduction_raw"][idx] = 1-ss/np.maximum(raw0,1e-20)
                a = error_w[:-1]-error_w[:-1].mean(axis=0)
                b = error_w[1:]-error_w[1:].mean(axis=0)
                stats[name]["white_residual_acf1"][idx] = np.sum(a*b,axis=0)/np.maximum(np.sqrt(squared(a)*squared(b)),1e-20)
                return sw

            record("N0",fit0,wz@coef0,coef0)
            stats["N0"]["partial_r2_white"][idx] = 0
            stats["N0"]["sse_reduction_raw"][idx] = 0
            for name,x in raw.items():
                wx,px = operators[name]
                coef = px@yw
                sw = record(name,x@coef,wx@coef,coef)
                if name == "GLM_TD":
                    c = designs.static_contrast
                    effect[idx] = c@coef
                    variance[idx] = sw/(len(x)-x.shape[1])*np.sum((c@px)**2)
                else:
                    c = designs.fir_contrasts
                    cv = c@coef
                    cov = (c@px)@(c@px).T
                    fir_f[idx] = np.sum(cv*np.linalg.solve(cov,cv),axis=0)/len(c)/np.maximum(sw/(len(x)-x.shape[1]),1e-20)
            if g is not None:
                yr = yw-wo@(po@yw)
                dot = residual_g.T@yr
                beta = dot/den[:,None]
                if selection == "raw_r2":
                    e0 = yy-other@(po@yw)
                    raw_rss = squared(e0)[None,:]-2*beta*(g_raw.T@e0)+beta**2*squared(g_raw)[:,None]
                    selected = np.argmin(raw_rss,axis=0)
                    if default_choice is not None:
                        denom = np.maximum(sst,1e-20)
                        default_r2[idx] = 1-raw_rss[default_choice]/denom
                        raw_gain[idx] = (raw_rss[default_choice]-raw_rss[selected,np.arange(len(idx))])/denom
                else:
                    selected = np.argmax(dot*beta,axis=0)
                b = beta[selected,np.arange(len(idx))]
                other_coef = po@(yw-wg[:,selected]*b)
                best[idx] = selected
                record("DoubleScanner",other@other_coef+g[:,selected]*b,
                       wo@other_coef+wg[:,selected]*b,np.vstack((b,other_coef)))

    dof = len(y)-static.shape[1]
    if (variance<=0).any() or not np.isfinite(variance).all():
        raise StageError("Static contrast contains invalid variance; review fit or mask")
    tv = effect/np.sqrt(variance)
    return dict(metrics=stats,coefficients=coefs,effect=effect,variance=variance,t=tv,
                p_two_sided=2*student_t.sf(abs(tv),dof),dof=dof,
                fir_f=fir_f,fir_p=(fisher_f.sf(fir_f,len(designs.fir_contrasts),len(y)-fir.shape[1])
                                   if fir is not None else None),best=best,default_r2=default_r2,raw_gain=raw_gain)


def bh(p):
    p = np.asarray(p,dtype=float)
    if p.ndim!=1 or not len(p) or not np.isfinite(p).all() or (p<0).any() or (p>1).any():
        raise StageError("Invalid p values for the declared voxel family")
    order = np.argsort(p,kind="stable")
    q = np.empty_like(p)
    q[order] = np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
    return q
