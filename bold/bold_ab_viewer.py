
"""Report-only DoubleScanner A/B figures and 4D payloads from a completed fit."""
import base64
from pathlib import Path
import shutil

import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
from scipy.ndimage import binary_erosion, gaussian_filter1d
from scipy.signal import find_peaks

import bold_workflow as w
from bold_runs import file_for
from bold_stage import same_grid, read_mask

MODES = ("A","A_stable","B","B_stable")


def encode(values,dtype):
    return base64.b64encode(np.ascontiguousarray(values,dtype=dtype).tobytes()).decode("ascii")


def timing_groups(f,reference,mask,folder):
    """Peak groups for common starts/ends relative to the modeled event."""
    durations = f["target_events"][1]
    if not np.allclose(durations,durations[0],rtol=0,atol=1e-8):
        return [],"可变长度事件保留逐帧模型，无法为全部事件定义同一相对终点。"
    p = f["parameters"]
    ix = f["voxel"]
    choices = f["choices"]
    final = f["final"]
    groups = []
    for mode in ("A","B"):
        selected = np.flatnonzero(f[mode+"_stable"][ix])
        if not len(selected):
            continue
        eligible = ix[selected]
        minimum = max(15,int(np.ceil(.02*len(eligible))))
        for landmark in ("onset","offset"):
            offset = durations[0]*p[:,1] if landmark=="offset" else 0
            ticks = np.rint(2*(p[:,0]+offset)).astype(int)
            t0 = ticks.min()
            observed = ticks[final[eligible]]
            hist = np.bincount(observed-t0,minlength=ticks.max()-t0+1)
            density = gaussian_filter1d(hist.astype(float),2)
            peaks,_ = find_peaks(np.r_[0.,density,0.],distance=6,prominence=max(2.,.05*density.max()))
            centers = sorted((k-1+t0 for k in peaks),key=lambda k:(-density[k-t0],k))
            assigned = np.zeros(len(eligible),bool)
            for center in centers:
                nearby = np.flatnonzero((abs(observed-center)<=2)&~assigned)
                if len(nearby)<minimum:
                    continue
                timing = ticks[choices[:,selected[nearby]]]
                frequency = (abs(timing-center)<=3).mean(axis=0)
                keep = nearby[frequency>=.8]
                if len(keep)<minimum:
                    continue
                assigned[keep] = True
                points = np.flatnonzero(mask)[eligible[keep]]
                series = sum(g["mode"]==mode and g["landmark"]==landmark for g in groups)+1
                label = f"{mode}_{landmark}_{series:02d}"
                name = label+".nii.gz"
                volume = np.zeros(mask.shape,np.uint8)
                volume.ravel()[points] = 1
                nib.save(nib.Nifti1Image(volume,reference.affine),folder/name)
                groups.append(dict(name=label,mode=mode,landmark=landmark,n=len(points),
                                   center=float(center/2),points=encode(points,"<i4"),mask_file=name))
    return groups,"共同起点和终点分开寻峰，时间相对已登记的任务参考。"


def plot_choices(f,folder):
    p,best,reason = f["parameters"],f["best"],f["reason"]
    lag,ratio = np.unique(p[:,0]),np.unique(p[:,1])
    dx,dy = (lag[1]-lag[0])/2,(ratio[1]-ratio[0])/2
    extent = (lag[0]-dx,lag[-1]+dx,ratio[0]-dy,ratio[-1]+dy)
    final = f["final"]
    scan = reason==0
    figures = [
        ("raw_argmax.svg",[("Raw argmax",best),("Edge fallback",best[(reason>=1)&(reason<=4)]),
                          ("Effect-size fallback",best[reason==5])],(1,3)),
        ("ab_heatmaps.svg",[(mode,final[f[mode]&scan]) for mode in ("A","A_stable","B","B_stable")],(2,2)),
    ]
    for filename,panels,layout in figures:
        fig,axes = plt.subplots(*layout,figsize=(10,3.6*layout[0]),layout="constrained",sharex=True,sharey=True)
        arrays = [np.bincount(index,minlength=len(p)).reshape(len(ratio),len(lag)) for _,index in panels]
        maximum = max(1,max(int(a.max()) for a in arrays))
        for ax,(name,index),array in zip(axes.flat,panels,arrays):
            im = ax.imshow(array,origin="lower",extent=extent,aspect="auto",cmap="Blues",vmin=0,vmax=maximum)
            ax.set(title=f"{name} (n={len(index)})",xlabel="Lag (s)",ylabel="Duration ratio")
        fig.colorbar(im,ax=axes.ravel().tolist(),label="Voxels",shrink=.8)
        fig.savefig(folder/filename)
        plt.close(fig)


def payload(f,mask,reference,anatomy,anatomy_mask,groups,run):
    im = nib.load(anatomy)
    if len(im.shape)!=3 or not same_grid(im,reference) or nib.aff2axcodes(im.affine)!=("R","A","S"):
        raise w.StageError("4D viewer requires T1 and BOLD on an identical RAS grid")
    brain = read_mask(anatomy_mask,reference)
    data = im.get_fdata(dtype=np.float32,caching="unchanged")
    if not brain.any() or not np.isfinite(data[brain]).all():
        raise w.StageError("4D viewer has an invalid T1 brain mask")
    lo,hi = np.quantile(data[brain],[.01,.99])
    if hi<=lo:
        raise w.StageError("4D viewer has no T1 intensity range")
    bg = np.where(brain,25+220*np.clip((data-lo)/(hi-lo),0,1),0).astype(np.uint8)
    shell = np.flatnonzero(brain&~binary_erosion(brain))
    idx = np.flatnonzero(f["B"])
    vox = np.flatnonzero(mask)[idx]
    win = f["final"][idx]
    beta,t = f["beta"][idx],f["t"][idx]
    p,g,clock = f["parameters"],f["templates"],f["clock"]
    flags = sum((f[key][idx].astype(np.uint8)<<j) for j,key in enumerate(MODES))
    box = np.zeros((len(clock),len(p)),np.float32)
    for onset,duration,amplitude in f["target_events"].T:
        start = onset+p[:,0]
        box += amplitude*((clock[:,None]>=start)&(clock[:,None]<start+duration*p[:,1]))
    post = np.abs(beta)*np.max(abs(g[:,win]),axis=0)
    scale = lambda x: max(float(np.quantile(x,.95)),1e-6) if len(x) else 1.
    return dict(run=run["id"],task=run["task_id"],shape=list(mask.shape),affine=reference.affine.tolist(),
                anatomy=encode(bg,"u1"),shell=encode(shell[::max(1,int(np.ceil(len(shell)/5200)))],"<i4"),
                voxel=encode(vox,"<i4"),choice=encode(win,"<u2"),beta=encode(beta,"<f4"),t=encode(t,"<f4"),
                reason=encode(f["reason"][idx],"u1"),bits=encode(flags,"u1"),g=encode(g,"<f4"),box=encode(box,"<f4"),
                frame=clock.tolist(),params=p.tolist(),groups=groups,counts={key:int(f[key].sum()) for key in MODES},
                pre_scale=scale(np.abs(beta)),post_scale=scale(post))


def render_ab(config,run,fit,folder,dest,registered):
    import json
    out = dest/"ab"/run["id"]
    out.mkdir(parents=True)
    reference = nib.load(fit/"mask.nii.gz")
    mask = read_mask(fit/"mask.nii.gz",reference)
    with np.load(folder/"double_scanner_ab.npz") as f:
        plot_choices(f,out)
        groups,note = timing_groups(f,reference,mask,out)
        if not run.get("viewer_anatomy"):
            return out,groups,note,"未配置与 BOLD 同网格的已登记 T1 和掩膜。"
        try:
            anatomy = file_for(config,registered,run,"viewer_anatomy")
            anatomy_mask = file_for(config,registered,run,"viewer_anatomy_mask")
            data = payload(f,mask,reference,anatomy,anatomy_mask,groups,run)
        except w.StageError as exc:
            return out,groups,note,str(exc)
        (out/"viewer_data.js").write_text("window.BOLD_AB="+json.dumps(data,allow_nan=False,separators=(",",":"))+";",
                                          encoding="utf-8")
    code = Path(__file__).parent
    for filename in ("ab_viewer.js","ab_viewer.css"):
        target = out.parent/filename
        if not target.exists():
            shutil.copyfile(code/filename,target)
    shutil.copyfile(code/"ab_viewer.html",out/"viewer.html")
    return out,groups,note,None
