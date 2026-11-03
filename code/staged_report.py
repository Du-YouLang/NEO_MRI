"""HTML and figures from SAVED results only; never imports a fitting engine."""
from pathlib import Path
from html.parser import HTMLParser
from urllib.parse import unquote,urlsplit
import html,json,os
from staged_workflow import read,write,resolve,derivatives,source_root

def render(c,pre,stat,reason,dest):
    e=html.escape
    def link(p,title):return '<a href="'+e(os.path.relpath(p,dest).replace('\\','/'),quote=True)+'">'+e(title)+'</a>'
    context=c.get('report',{})
    body='<h1>'+e(context.get('title',c['case_id']))+'</h1><p>分阶段流程 · '+e(c['case_id'])+'</p>'
    body+='<section><h2>这次运行做了什么</h2><p>预处理整理解剖、头动与坐标空间；GLM把已预处理BOLD与任务时序、噪声候选量一起拟合；本页只读取保存的结果来绘图和解释。生成本页没有启动fMRIPrep，也没有重新拟合GLM。</p>'
    body+='<p>原结果：<code>'+e(str(source_root(c)))+'</code>。所有新状态、统计版本与本页存放于另一目录，原始和已完成患者结果保持原字节。</p><p>新执行固定OMP/ITK/BLAS为4线程；当前引用的历史预处理实际OMP线程为'+e(str(pre['historical_actual_omp_threads']))+'。历史线程记录不被改写，也不会为了统一到4线程而重算旧患者。</p>'
    body+='<p>预处理文件指纹：<code>'+pre['fingerprint']+'</code>。'+link(resolve(c['workspace'])/'preprocessing/receipt.json','预处理登记与哈希')+'</p></section>'
    if context.get('source_report'):
        p=resolve(context['source_report'])
        if p.is_file():body+='<section><h2>已有完整实验报告</h2><p>'+link(p,'打开已有结构、病灶和BOLD质控报告')+'</p><p>本页是分阶段流程的新入口；完整的原实验解释和已审阅图像仍保留在上述报告。</p></section>'
    if stat is None:
        body+='<section class="note"><h2>当前没有可展示的新任务统计结果</h2><p>'+e(reason or '尚未拟合')+'</p><p>本页没有把BOLD均值、tSNR、病灶掩膜或已有健康人激活当作这位患者的新任务激活。</p></section>'
    else:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import nibabel as nib
        import numpy as np
        import pandas as pd
        from nilearn import plotting
        result=read(stat/'summary.json')
        if (stat/'model_audit.json').is_file():
            body+='<p>'+link(stat/'model_audit.json','设计矩阵审计：秩、名义自由度、共线性与单位噪声OLS对比方差')+'。这是设计审计，不是过拟合诊断或实际AR(1)统计功效。</p>'
        body+='<section><h2>当前统计版本</h2><p>'+e(result['scientific_scope'])+'</p><p>'+e(c['analysis']['timing']['interpretation'])+'</p><p>'+link(stat/'request.json','实际输入与统计参数指纹')+' · '+link(stat/'summary.json','统计与软件记录')+'</p><p>颜色显示模型的Z统计量；它依赖任务时序与噪声模型，不是直接测得的神经元活动强度。FDR范围是本次共同脑掩膜，工程子集测试的FDR仅限子集，不能作为全脑推断。</p>'
        s=result['combined'];body+='<p>共同掩膜体素数：'+str(s['analysis_voxels'])+'；双侧FDR q='+str(s['fdr_q'])+'；正/负方向达阈值体素：'+str(s['positive_voxels'])+' / '+str(s['negative_voxels'])+'。阈值Z：'+e(str(s['z_threshold']))+'（None表示没有发现）。</p>'
        for r in result['runs']:
            d=r['diagnostics'];body+='<p>'+e(r['id'])+'：'+str(d['volumes'])+'帧，TR='+f"{d['TR_seconds']:.3g}"+'秒；设计矩阵秩'+str(d['rank'])+'，名义残差自由度'+str(d['nominal_residual_dof'])+'；实际选用combined aCompCor '+str(len(d['combined_compcor']))+'项，异常帧回归量'+str(len(d['spikes']))+'项。自由度是模型可估计性的线索，不是数据质量评分。</p>'
        anatomy=derivatives(c)/c['analysis']['anatomy']
        mask_image=nib.load(stat/'common_mask.nii.gz')
        center=nib.affines.apply_affine(mask_image.affine,np.argwhere(mask_image.get_fdata()>.5).mean(axis=0))
        cuts=context.get('cut_coords_world_mm',center.tolist())
        for name,threshold in [('z.nii.gz',None),('z_fdr.nii.gz',1e-10)]:
            im=nib.load(stat/'combined'/name)
            d=plotting.plot_stat_map(im,bg_img=str(anatomy),display_mode='ortho',cut_coords=cuts,threshold=threshold,
                black_bg=False,title='Combined saved Z map; '+('uncorrected display' if threshold is None else 'two-sided FDR'))
            filename='combined_'+name.replace('.nii.gz','.png');d.savefig(str(dest/filename),dpi=135);d.close()
            body+='<figure><img src="'+filename+'"><figcaption>'+('未阈值化Z图，不能按颜色直接认定显著。' if threshold is None else '双侧FDR阈值后的Z图；没有着色可能是未达阈值、未覆盖或发现不在这些切面，不能直接认定组织不可用。')+' 两图使用相同的世界坐标切面（默认共同脑掩膜中心），图中坐标不一定是MNI。</figcaption></figure>'
        for run in result['runs']:
            name=run['id'];design=pd.read_csv(stat/name/'design.tsv',sep='\t',index_col=0)
            fig,ax=plt.subplots(figsize=(12,5),constrained_layout=True)
            a=design.to_numpy();scaled=a/np.maximum(np.abs(a).max(axis=0),1e-12)
            ax.imshow(scaled,aspect='auto',cmap='RdBu_r',vmin=-1,vmax=1,interpolation='nearest')
            ax.set_xticks(np.arange(len(design.columns)),design.columns,rotation=90,fontsize=6);ax.set_ylabel('Saved volume index')
            ax.set_title(name+' saved design; rank '+str(run['diagnostics']['rank'])+', nominal DOF '+str(run['diagnostics']['nominal_residual_dof']))
            f=name+'_design.png';fig.savefig(dest/f,dpi=125);plt.close(fig)
            body+='<figure><img src="'+f+'"><figcaption>每列是一项模型解释量，每行对应一个保存体积。任务列与头动等列过度重合会影响估计；本图只按列缩放以便显示。</figcaption></figure>'
        body+='</section>'
    body+='<section><h2>质量与使用边界</h2><p>数值完成与临床认可分别记录；本页自动生成不等于人工图像审阅通过。头动、覆盖不足、卒中后的血流动力学改变都可能影响BOLD解释。本流程不自动计算皮层可用率、CST损伤率或临床入组资格。</p><p>仅改统计参数请运行glm，再运行report；仅改报告样式请只运行report。报告或GLM都不会自动回退启动预处理。</p></section>'
    style='body{font:17px/1.85 "Microsoft YaHei",system-ui,sans-serif;background:#eef3f6;color:#233b4c;margin:0}main{max-width:1180px;margin:auto;padding:28px}section{background:white;border-radius:12px;padding:23px;margin:22px 0}h1,h2{color:#1c5879}a{color:#076b8d}.note{border-left:6px solid #c39330}code{overflow-wrap:anywhere}img{max-width:100%}figure{margin:24px 0}figcaption{font-size:15px}'
    page=dest/'report.html';page.write_text('<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><style>'+style+'</style><main>'+body+'</main></html>',encoding='utf-8')
    class Parser(HTMLParser):
        def __init__(self):super().__init__();self.refs=[]
        def handle_starttag(self,t,a):self.refs.extend(v for k,v in a if k in ['href','src'])
    parser=Parser();parser.feed(page.read_text(encoding='utf-8'));count=0
    for ref in parser.refs:
        u=urlsplit(ref)
        if not u.scheme and u.path:
            if not (dest/unquote(u.path)).is_file():raise RuntimeError('Report target missing: '+ref)
            count+=1
    write(dest/'static_validation.json',dict(local_links_checked=count,GLM_fits=0,fmriprep_invocations=0,browser_interaction_verified=False))
