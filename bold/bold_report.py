
"""Generate the standard BOLD HTML from persisted outputs alone."""
import html
import json
import os

import nibabel as nib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import bold_workflow as w


def escape(value):
    return html.escape(str(value),quote=True)


def link(path, base, label):
    from urllib.parse import quote
    href = quote(os.path.relpath(path,base).replace(os.sep,"/"),safe="/")
    return f'<a href="{href}">{escape(label)}</a>'


def details(title,data):
    return f"<details><summary>{escape(title)}</summary><pre>{escape(json.dumps(data,ensure_ascii=False,indent=2))}</pre></details>"


def metrics_table(models):
    metrics = [("r2","原始 R²"),("rmse","原始 RMSE（%）"),("partial_r2_white","白化后相对 N0 的部分 R²"),
               ("sse_reduction_raw","原始域相对 N0 的 SSE 降幅"),("white_residual_acf1","白化残差 lag-1 相关")]
    rows = []
    for key,title in metrics:
        for name,values in models.items():
            row = values[key]
            numbers = "".join(f"<td>{row[k]:.4g}</td>" for k in ("mean","q25","median","q75"))
            rows.append(f"<tr><td>{title}</td><td>{name}</td><td>{row['n']}</td>{numbers}</tr>")
    return ("<table><thead><tr><th>指标</th><th>模型</th><th>体素数</th><th>均值</th>"
            "<th>Q25</th><th>中位数</th><th>Q75</th></tr></thead><tbody>"+"".join(rows)+"</tbody></table>")


def histograms(payload):
    names = {"r2":"原始 R²","rmse":"原始 RMSE","partial_r2_white":"白化后部分 R²",
             "sse_reduction_raw":"原始 SSE 降幅","white_residual_acf1":"白化残差 lag-1"}
    sections = ["<details><summary>模型体素分布（每个指标内共享分箱和坐标轴）</summary>"]
    for metric,views in payload.items():
        for view,distribution in views.items():
            edges = distribution["edges"]
            series = distribution["series"]
            width = 405/len(next(iter(series.values()))["counts"])
            y_max = max(1,max(max(model["counts"]) for model in series.values()))
            title = names[metric] + (" · 全范围" if view=="full" else " · 共有 1–99% 区间")
            sections.append(f"<details><summary>{escape(title)}</summary>")
            for name,model in series.items():
                bars = "".join(f'<rect x="{55+i*width:.3f}" y="{55-44*value/y_max:.3f}" '
                               f'width="{width-.4:.3f}" height="{44*value/y_max:.3f}" fill="#37809a"/>'
                               for i,value in enumerate(model["counts"]))
                tails = model["below"]+model["above"]
                total = sum(model["counts"])+tails
                sections.append(f'<p>{escape(name)} · 低于范围 {model["below"]}，高于范围 {model["above"]}'
                                f'（合计 {tails/total:.1%}）</p><svg viewBox="0 0 485 74" width="485" '
                                f'aria-label="{escape(name)} 直方图">{bars}<text x="7" y="15" font-size="11">{y_max}</text>'
                                f'<text x="55" y="69" font-size="11">{edges[0]:.3g}</text>'
                                f'<text x="460" y="69" text-anchor="end" font-size="11">{edges[-1]:.3g}</text></svg>')
            sections.append("</details>")
    return "".join(sections)+"</details>"


def comparison_figure(folder, output):
    mask = nib.as_closest_canonical(nib.load(folder.parent/"mask.nii.gz")).get_fdata()>0
    names = ("left","right","left_minus_right")
    data = [nib.as_closest_canonical(nib.load(folder/(name+".nii.gz"))).get_fdata() for name in names]
    limit = max(float(np.quantile(np.abs(np.concatenate([d[mask] for d in data])),.98)),1e-8)
    center = np.argwhere(mask).sum(axis=0)//mask.sum()
    zooms = nib.as_closest_canonical(nib.load(folder.parent/"mask.nii.gz")).header.get_zooms()
    fig, axes = plt.subplots(3,3,figsize=(9,8),layout="constrained")
    for i,axisrow in enumerate(axes):
        remaining = [j for j in range(3) if j!=i]
        for j,ax in enumerate(axisrow):
            plane = np.ma.masked_where(~np.take(mask,center[i],axis=i).T,np.take(data[j],center[i],axis=i).T)
            artist = ax.imshow(plane,origin="lower",cmap="RdBu_r",vmin=-limit,vmax=limit,
                               aspect=zooms[remaining[1]]/zooms[remaining[0]],interpolation="nearest")
            ax.set_axis_off()
            if i==0:
                ax.set_title(("Left beta","Right beta","Left - right")[j])
    fig.colorbar(artist,ax=axes.ravel().tolist(),label="Beta (% signal change)",shrink=.65,extend="both")
    fig.savefig(output,dpi=115)
    plt.close(fig)


PAGE = """<!doctype html><html lang="zh-CN"><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{title}</title><style>
body{{font:16px/1.55 system-ui,sans-serif;color:#172a35;background:#f3f7f9;margin:0 auto;padding:2rem 1.2rem;max-width:1100px}}
section{{background:white;margin:1rem 0;padding:1rem 1.2rem;border:1px solid #d7e1e6;border-radius:9px}}
h1,h2,h3{{color:#114b63}} a{{color:#035e83}} table{{border-collapse:collapse;width:100%}}
td,th{{border-bottom:1px solid #d7e1e6;padding:.25rem .45rem;text-align:left}} pre{{white-space:pre-wrap;overflow-wrap:anywhere}}
details{{margin:.65rem 0}} img{{width:min(100%,780px)}} svg{{max-width:100%;background:#f6f9fa}} label{{margin-right:.75rem}}
</style><h1>{title}</h1>
<p>每条 run 独立拟合。静态主结果采用 GLM 加时间导数，报告指定条件的 canonical β 对比。DoubleScanner 为动态主结果，报告选中模板的 β、延迟与时长比例。FIR 联合 F 检验作为敏感性对照。</p>
<p>GLM 和 FIR 的体素级 BH 校正各自在对应 run 或同任务联合、平滑尺度及分析掩膜内计算。DoubleScanner 的模板是从同一批数据选择的，本版只报告选后描述量。原始域和白化域的模型指标用于描述拟合表现，均来自拟合样本。</p>
{body}
<script type="application/json" id="trace-data">{traces}</script>
<script>
const traces = JSON.parse(document.querySelector('#trace-data').textContent);
const colors = {{observed:'#111827',N0:'#64748b',GLM_TD:'#1565ad',DoubleScanner:'#a33c08',FIR:'#16835c',dynamic_post_hrf:'#a855a5'}};
for(const widget of document.querySelectorAll('[data-trace]')){{
  const data=traces[widget.dataset.trace], plot=widget.querySelector('svg'), controls=widget.querySelector('.controls');
  const keys=Object.keys(colors).filter(k=>Array.isArray(data[k]));
  for(const k of keys){{
    const label=document.createElement('label'), input=document.createElement('input');
    input.type='checkbox';input.value=k;input.checked=k!=='N0'&&k!=='dynamic_post_hrf';
    label.style.color=colors[k];label.append(input,document.createTextNode(k));controls.append(label);
    input.addEventListener('change',draw);
  }}
  const slider=widget.querySelector('input[type=range]'), info=widget.querySelector('output');
  slider.max=data.seconds.length-1;slider.addEventListener('input',draw);
  function draw(){{
    const active=Array.from(controls.querySelectorAll('input:checked'),x=>x.value);
    const values=active.flatMap(k=>data[k]); const low=Math.min(0,...values),high=Math.max(0,...values);
    const span=Math.max(high-low,1e-6), t0=data.seconds[0], t1=data.seconds.at(-1);
    const x=t=>75+665*(t-t0)/(t1-t0), y=v=>174-150*(v-low)/span;
    function node(type,attributes,text){{
      const el=document.createElementNS('http://www.w3.org/2000/svg',type);
      for(const [name,value] of Object.entries(attributes))el.setAttribute(name,value);
      if(text!==undefined)el.textContent=text;plot.append(el);
    }}
    plot.replaceChildren();
    for(const [on,length] of data.dynamic_event_windows||[]){{
      const left=x(Math.max(t0,on)),right=x(Math.min(t1,on+length));
      if(right>left)node('rect',{{x:left,y:24,width:right-left,height:150,fill:'#d9b685',opacity:'.32'}});
    }}
    for(const v of [low,(low+high)/2,high]){{
      node('line',{{x1:75,x2:740,y1:y(v),y2:y(v),stroke:'#dae2e7'}});
      node('text',{{x:8,y:y(v)+4,'font-size':11}},v.toFixed(3)+' %');
    }}
    for(const k of active)node('polyline',{{points:data[k].map((v,j)=>x(data.seconds[j])+','+y(v)).join(' '),fill:'none',stroke:colors[k],'stroke-width':1.5}});
    const i=Number(slider.value), at=data.seconds[i];
    node('line',{{x1:x(at),x2:x(at),y1:24,y2:174,stroke:'#6b7280','stroke-dasharray':'4 3'}});
    node('text',{{x:360,y:196,'font-size':12}},'Acquisition time (s)');
    info.textContent='第 '+(i+1)+' 帧，'+at.toFixed(2)+' 秒。 '+active.map(k=>k+' '+data[k][i].toFixed(3)+'%').join(' ｜ ');
  }}
  draw();
}}
</script></html>"""


def render(config,fit,comparison,pending,dest):
    summary = w.read(fit/"summary.json")
    configured = {run["id"]:run for run in config["analysis"]["runs"]}
    pre = w.read(w.workspace(config)/"preprocessing/receipt.json")
    registered = {r["path"] for r in pre["files"]}
    traces = {}
    sections = []
    for row in summary["runs"]:
        name = row["id"]
        parts = [f"<section><h2>{escape(name)} · {escape(row.get('task_id') or '任务待确认')}</h2>"]
        if row["status"]!="FITTED":
            sections.append("".join(parts)+f"<p>未拟合。{escape(row['reason'])}</p></section>")
            continue
        run = fit/name
        fitted = w.read(run/"summary.json")
        proof = row["timing"]
        parts.append(f"<p>{row['frames']} 帧，TR {row['tr_seconds']:g} 秒，首帧模型时刻 {row['start_seconds']:g} 秒，"
                     f"STC {escape(row['stc'])}。时序状态 {escape(proof['status'])}。{escape(proof['interpretation'])}</p>")
        parts.append(details("时序、模型设计与混杂审计",row))
        parts.append(details("静态对比、动态候选和 FIR 时滞",fitted))
        if fitted.get("ab_unavailable"):
            parts.append(f"<p>A/B 未运行。{escape(fitted['ab_unavailable'])}</p>")
        for kernel,entry in fitted["kernels"].items():
            folder = run/f"fwhm_{kernel}mm"
            trace_id = f"{name}-{kernel}"
            trace = w.read(folder/"trace.json")
            traces[trace_id] = trace
            mark = "（主报告尺度）" if int(kernel)==fitted["primary_kernel_mm"] else ""
            parts.extend([f"<h3>{escape(kernel)} mm{mark}</h3>",metrics_table(entry["models"]),
                          histograms(w.read(folder/"histograms.json")),
                          "<p>"+link(folder/"effect.nii.gz",dest,"静态 β")+" · "+
                          link(folder/"q_bh.nii.gz",dest,"静态 q")+" · "+
                          link(folder/"coefficients.npz",dest,"模型系数")+"</p>"])
            if "DoubleScanner" in entry["models"]:
                parts.append("<p>"+" · ".join(link(folder/(file+".nii.gz"),dest,label)
                    for file,label in (("ds_beta","动态 β"),("ds_lag_seconds","动态延迟（秒）"),("ds_duration_ratio","时长比例")))+"</p>")
            if "DoubleScanner_AB" in entry:
                from bold_ab_viewer import render_ab
                output,groups,note,issue = render_ab(config,configured[name],run,folder,dest,registered)
                ab = entry["DoubleScanner_AB"]
                parts.extend(["<h3>DoubleScanner A/B</h3>",
                              "<p>原空间 R² 选点。触边或未达到两项效应量要求的体素回到 (0 秒, 1 倍)。"
                              "同一 run 全脑有效掩膜内，扫描与回退共用一个 BY 家族。"
                              "A 使用 min(1, 1214 × p)，B 使用选后名义 p，仅供探索。"
                              "A 的 1215 倍另作核验。推断以估计后固定的 AR(1) 为条件。"
                              "上方 DoubleScanner 拟合指标与时间线显示回退前赢家，以下为回退后分支。</p>",
                              details("家族、阈值、回退与共享 Bootstrap",ab),
                              "<p>"+link(folder/"ab_final_beta.nii.gz",dest,"最终 β")+" · "+
                              link(folder/"ab_model_source.nii.gz",dest,"回退与模型来源")+" · "+
                              link(folder/"ab_q_A.nii.gz",dest,"A 的 BY q")+" · "+
                              link(folder/"ab_q_B_exploratory.nii.gz",dest,"B 的 BY q（探索性）")+"</p>",
                              f'<p>{escape(note)}</p>',
                              "<p>"+link(output/"raw_argmax.svg",dest,"回退前 argmax 热图")+" · "+
                              link(output/"ab_heatmaps.svg",dest,"A/B 筛选与稳定热图")+"</p>"])
                if issue:
                    parts.append(f"<p>三视图与 4D 未生成。{escape(issue)}</p>")
                else:
                    parts.append("<p>"+link(output/"viewer.html",dest,"打开同步三视图与 4D")+" </p>")
                if groups:
                    parts.append("<ul>")
                    for group in groups:
                        label = "共同起点" if group["landmark"]=="onset" else "共同终点"
                        parts.append("<li>"+link(output/group["mask_file"],dest,
                                     f'{group["mode"]} {label} {group["center"]:g} 秒，{group["n"]} 体素')+"</li>")
                    parts.append("</ul>")
            if "FIR" in entry:
                parts.append("<p>"+link(folder/"fir_joint_f.nii.gz",dest,"FIR 联合 F")+" · "+
                             link(folder/"fir_joint_q_bh.nii.gz",dest,"FIR 联合 q")+"</p>")
            parts.append(f'<div data-trace="{escape(trace_id)}"><p>同一代表体素 IJK {escape(trace["voxel"])}。'
                         '竖线标示真实采样帧，底色为 DoubleScanner 选中的 HRF 卷积前任务窗口。'
                         '勾选 dynamic_post_hrf 可查看卷积后的任务分量。</p><div class="controls"></div>'
                         '<svg viewBox="0 0 760 203" aria-label="实际采样与拟合时间序列"></svg><br>'
                         '<label>浏览帧 <input type="range" min="0" value="0"></label><output></output></div>')
        sections.append("".join(parts)+"</section>")
    if summary["groups"]:
        parts = ["<section><h2>同一任务重复 run 联合</h2>",
                 "<p>在共有物理网格及掩膜内按 run 方差逆数合并静态 β。动态与 FIR 保留各 run 结果。</p>"]
        for group in summary["groups"]:
            parts.append(details(group["id"],group))
            for kernel in group.get("kernels",{}):
                parts.append(f"<p>{escape(kernel)} mm "+link(fit/group["id"]/f"fwhm_{kernel}mm"/"effect.nii.gz",dest,"联合静态 β")+"</p>")
        sections.append("".join(parts)+"</section>")
    if pending:
        sections.append("<section><h2>已配置的 β 对照待生成</h2>"+details("比较请求",pending)+"</section>")
    if comparison is not None:
        parts = ["<section><h2>拟合后 β 对照</h2><p>差值为左减右。不同任务的效应大小受条件和模板定义影响。图上两侧和差值使用同一色轴。</p>"]
        for row in w.read(comparison/"summary.json")["comparisons"]:
            parts.append(details(row["id"],row))
            for kernel in row.get("kernels",{}):
                folder = comparison/row["id"]/f"fwhm_{kernel}mm"
                image = dest/f'{row["id"]}_{kernel}mm.png'
                comparison_figure(folder,image)
                parts.append(f'<p>{escape(kernel)} mm · {escape(row["model"])} · {escape(row["left"])} − {escape(row["right"])}</p>'
                             f'<img src="{escape(image.name)}" alt="两张效应图与共同掩膜内的差值切片">'
                             "<p>"+link(folder/"left_minus_right.nii.gz",dest,"差值 NIfTI")+"</p>")
        sections.append("".join(parts)+"</section>")
    title = config.get("report",{}).get("title",f'{config["case_id"]} · 任务 BOLD 分析')
    serialized = json.dumps(traces,ensure_ascii=False,allow_nan=False).replace("<","\\u003c")
    (dest/"report.html").write_text(PAGE.format(title=escape(title),body="".join(sections),traces=serialized),encoding="utf-8")
