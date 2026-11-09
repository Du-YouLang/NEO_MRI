
"""Opt-in raw-R² DoubleScanner A/B protocol and conditional stability."""
import numpy as np
from scipy.stats import false_discovery_control, t as student_t

from bold_models import squared, whiten
from bold_workflow import StageError

LAGS = tuple(x/2 for x in range(-24,21))
RATIOS = tuple(x/20 for x in range(4,31))
PARAMETERS = np.asarray([(lag,ratio) for ratio in RATIOS for lag in LAGS])
FACTORS = {"A":1214,"B":1,"audit_1215":1215}
TASK_GAIN = .05
DEFAULT_GAIN = .01
ALPHA = .05
DRAWS = 1000
STABILITY = .80
REASONS = ("scan","lower_duration","upper_duration","earlier_lag","later_lag","effect_size")


def validate_policy(analysis):
    spec = analysis.get("double_scanner_ab")
    if spec is None:
        return
    if (set(spec)!={"seed","device"} or type(spec["seed"]) is not int or spec["seed"]<0
        or spec["device"] not in ("cpu","cuda")):
        raise StageError("DoubleScanner A/B requires an explicit nonnegative seed and cpu or cuda device")
    scan = analysis["double_scanner"]
    if (tuple(scan["lags_seconds"])!=LAGS or tuple(scan["duration_ratios"])!=RATIOS or scan["selection"]!="raw_r2"):
        raise StageError("A/B requires the declared 1,215-grid protocol and original-domain R2")
    if analysis["statistics"]["noise_model"]!="ar1":
        raise StageError("A/B uses the AR(1) profile with coefficients frozen from the unsmoothed static fit")


def classify(parameters, winner, partial, gain):
    lag,ratio = parameters.T
    # Duration is evaluated first so all four corners belong only to these edges.
    edges = np.select((ratio==ratio.min(),ratio==ratio.max(),lag==lag.min(),lag==lag.max()),(1,2,3,4))
    reason = edges[winner].astype(np.uint8)
    reason[(reason==0)&((partial<TASK_GAIN)|(gain<DEFAULT_GAIN))]=5
    default, = np.flatnonzero((parameters==[0.,1.]).all(axis=1))
    return edges,reason,np.where(reason==0,winner,default).astype(np.int16)


def chosen_gls(y,rho,z,g,choices):
    df = len(y)-z.shape[1]-1
    sse,beta,t = (np.empty(y.shape[1]) for _ in range(3))
    for r in np.unique(rho):
        indices = np.flatnonzero(rho==r)
        wz = whiten(z,r)
        pinv = np.linalg.pinv(wz)
        wg = whiten(g,r)
        residual_g = wg-wz@(pinv@wg)
        for start in range(0,len(indices),1024):
            ids = indices[start:start+1024]
            selected = choices[ids]
            target = residual_g[:,selected]
            yy = y[:,ids]
            yw = whiten(yy,r)
            energy = squared(target)
            if (energy<=1e-12).any():
                raise StageError("A/B chosen design is rank-deficient")
            wy = yw-wz@(pinv@yw)
            b = np.sum(wy*target,axis=0)/energy
            rest = pinv@(yw-wg[:,selected]*b)
            sw = squared(wy-target*b)
            sse[ids] = squared(yy-z@rest-g[:,selected]*b)
            beta[ids] = b
            t[ids] = b*np.sqrt(energy/np.maximum(sw/df,1e-20))
    return sse,beta,t,df


def corrected_family(t,df):
    """Submit exactly one two-sided nominal p per valid voxel, including fallback."""
    p = 2*student_t.sf(np.abs(t),df)
    if not len(p) or not np.isfinite(p).all():
        raise StageError("A/B requires a finite whole-mask p family")
    return p,{name:false_discovery_control(np.minimum(1.,factor*p),method="by")
              for name,factor in FACTORS.items()}


def bootstrap_choices(y,z,g,rho,anchor,draws,device):
    """Centered, df-scaled whitened residual bootstrap with shared row draws."""
    if len(anchor)==0:
        return np.empty((len(draws),0),np.int16)
    torch = None
    if device=="cuda":
        try:
            import torch
        except ImportError as exc:
            raise StageError("The optional CUDA PyTorch backend is not available") from exc
        if not torch.cuda.is_available():
            raise StageError("Requested CUDA Bootstrap backend is unavailable")
        torch.set_num_threads(4)

    result = np.empty((len(draws),len(anchor)),np.int16)
    n = len(y)
    df = n-z.shape[1]-1
    steps = np.arange(n)[:,None]-np.arange(n)[None,:]
    for r in np.unique(rho):
        group = np.flatnonzero(rho==r)
        wz,wg = whiten(z,r),whiten(g,r)
        p = np.linalg.pinv(wz)
        reg = p@wg
        gr = wg-wz@reg
        gr_raw = g-z@reg
        norm = np.sqrt(squared(gr))
        if (norm<=1e-6).any():
            raise StageError("A/B bootstrap contains a confounded candidate")
        inv = np.where(steps>=0,r**np.maximum(steps,0),0.)
        u = (gr.T/norm[:,None])
        u -= (u@wz)@p
        v = (gr_raw.T@inv-(gr_raw.T@z)@p)/norm[:,None]
        c = (squared(gr_raw)/norm**2)[:,None]
        if torch is not None:
            def device_array(values):
                return torch.as_tensor(values,dtype=torch.float64,device="cuda")
            cu,cv,cc = map(device_array,(u,v,c))
        for first in range(0,len(group),384):
            ix = group[first:first+384]
            yw = whiten(y[:,ix],r)
            choice = anchor[ix]
            task = wg[:,choice]*np.sum(u[choice].T*yw,axis=0)/norm[choice]
            fitted = task+wz@(p@(yw-task))
            err = yw-fitted
            err = (err-err.mean(axis=0))*np.sqrt(n/df)
            for start in range(0,len(draws),4):
                draw = draws[start:start+4]
                pseudo = (fitted[None]+err[draw]).transpose(1,0,2).reshape(n,-1)
                if torch is None:
                    a = u@pseudo
                    selection = (a*(2*(v@pseudo)-a*c)).argmax(axis=0)
                else:
                    x = device_array(pseudo)
                    a = cu@x
                    selection = (a*(2*(cv@x)-a*cc)).argmax(dim=0).cpu().numpy()
                result[start:start+len(draw),ix] = selection.reshape(len(draw),len(ix))
    return result


def stability_frequencies(parameters,edges,choices,reference):
    lag_step = np.min(np.diff(np.unique(parameters[:,0])))
    ratio_step = np.min(np.diff(np.unique(parameters[:,1])))
    integer = np.rint(parameters/[lag_step,ratio_step]).astype(np.int16)
    exact,near,interior,edge = (np.empty(len(reference)) for _ in range(4))
    for start in range(0,len(reference),512):
        where = slice(start,start+512)
        candidate = choices[:,where]
        orig = reference[where]
        close = np.abs(integer[candidate]-integer[orig]).sum(axis=-1)<=3
        inside = edges[candidate]==0
        exact[where] = (candidate==orig).mean(axis=0)
        near[where] = close.mean(axis=0)
        interior[where] = (close&inside).mean(axis=0)
        edge[where] = (~inside).mean(axis=0)
    return dict(exact=exact,near=near,near_inside=interior,edge_frequency=edge)


def fit_ab(y,designs,rho,fit,spec):
    """Consume the observed raw-R² argmax from the shared GLM and scan backend."""
    parameters = np.asarray(designs.parameters)
    if not np.array_equal(parameters,PARAMETERS):
        raise StageError("A/B candidate ordering or grid changed")
    g,z = designs.dynamic,designs.dynamic_nuisance
    if not np.array_equal(z,designs.nuisance.to_numpy()):
        raise StageError("A/B requires a single task condition so the 5% gate compares to the same N0")
    best = fit["best"]
    default_r2,gain = fit["default_r2"],fit["raw_gain"]
    if default_r2 is None or gain is None or (gain<0).any():
        raise StageError("Missing validated gain over the raw-R2 default candidate")
    observed = fit["metrics"]["DoubleScanner"]
    edge,reason,final = classify(parameters,best,observed["sse_reduction_raw"],gain)
    _,beta,t,df = chosen_gls(y,rho,z,g,final)
    p,q = corrected_family(t,df)
    a,b = q["A"]<=ALPHA,q["B"]<=ALPHA
    eligible = np.flatnonzero(b&(reason==0))
    draws = np.random.default_rng(spec["seed"]).integers(0,len(y),size=(DRAWS,len(y)),dtype=np.int64)
    choices = bootstrap_choices(y[:,eligible],z,g,rho[eligible],final[eligible],draws,spec["device"])
    freq = stability_frequencies(parameters,edge,choices,final[eligible])
    stable = np.zeros(len(a),bool)
    stable[eligible] = freq["near_inside"]>=STABILITY
    return dict(parameters=parameters,best=best,final=final,reason=reason,default_r2=default_r2,
                raw_best_r2=observed["r2"],raw_partial=observed["sse_reduction_raw"],raw_gain=gain,
                p=p,q_A=q["A"],q_B=q["B"],q_1215=q["audit_1215"],beta=beta,t=t,df=df,
                A=a,B=b,A_stable=a&stable,B_stable=b&stable,voxel=eligible,choices=choices,draws=draws,**freq)
