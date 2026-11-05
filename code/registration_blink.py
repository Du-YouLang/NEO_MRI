"""Blink saved normalized T1 against its exact template, after mask approval."""
from pathlib import Path
import argparse,base64,io,json,html
import nibabel as nib
import numpy as np
from PIL import Image,ImageDraw
from anatomy_review import require_approved,read,sha,check_file,save_new,file_record

def display_tile(slice_data, pixel_spacing, box=(310, 320)):
    """Fit a slice to screen, preserving physical aspect; no volume resampling."""
    image=Image.fromarray(slice_data).convert('RGB')
    physical=np.asarray(image.size)*np.asarray(pixel_spacing)
    scale=min(box[0]/physical[0],box[1]/physical[1])
    size=tuple(max(1,int(round(v*scale))) for v in physical)
    return image.resize(size,resample=Image.Resampling.NEAREST)


def render(config):
    review=require_approved(config['anatomy_review'])
    provenance=read(config['normalization_provenance'])
    expected=dict(t1_sha256=review['t1']['sha256'],mask_sha256=review['mask']['sha256'],
                  normalized_t1_sha256=sha(config['normalized_t1']),template_sha256=sha(config['template']))
    if any(provenance.get(k)!=v for k,v in expected.items()):raise ValueError('Normalization not bound to the approved T1/mask and selected template')
    if not provenance.get('transforms'):raise ValueError('Normalization transform evidence required')
    for item in provenance['transforms']:check_file(item)
    subject=nib.load(config['normalized_t1']);template=nib.load(config['template'])
    if len(subject.shape)!=3 or subject.shape!=template.shape or not np.allclose(subject.affine,template.affine,rtol=0,atol=1e-5):
        raise ValueError('Normalized image and template must already share the exact grid; viewer never registers/rescales')
    a=nib.as_closest_canonical(subject);b=nib.as_closest_canonical(template)
    point=np.array(config['cut_coords_world_mm'],float)
    if point.shape!=(3,) or not np.isfinite(point).all():raise ValueError('Three finite world coordinates required')
    index=np.rint(nib.affines.apply_affine(np.linalg.inv(b.affine),point)).astype(int)
    if (index<0).any() or (index>=b.shape).any():raise ValueError('World coordinate outside template FOV')
    seconds=float(config.get('interval_seconds',.7))
    if not .2<=seconds<=3:raise ValueError('Blink interval must be 0.2–3 seconds')
    images=[];window=[]
    for image,name in [(b,'TEMPLATE'),(a,'NORMALIZED SUBJECT')]:
        data=image.get_fdata(dtype=np.float32)
        if not np.isfinite(data).all():raise ValueError('Nonfinite anatomy')
        positive=data[data>0]
        if not len(positive):raise ValueError('Empty anatomy')
        lo,hi=np.percentile(positive,[1,99]);window.append([float(lo),float(hi)])
        if hi<=lo:raise ValueError('Invalid display range')
        scaled=np.clip((data-lo)/(hi-lo)*255,0,255).astype(np.uint8)
        canvas=Image.new('RGB',(960,390),(10,18,25));draw=ImageDraw.Draw(canvas)
        draw.text((15,8),name+' | '+config['space']+' | world mm '+str(point.tolist()),fill='white')
        for axis in range(3):
            sl=np.flipud(np.take(scaled,index[axis],axis=axis).T)
            remaining=[i for i in range(3) if i!=axis]
            im=display_tile(sl,[image.header.get_zooms()[i] for i in remaining])
            x=axis*320+(320-im.width)//2;y=45+(330-im.height)//2;canvas.paste(im,(x,y))
            draw.text((axis*320+10,370),['Sagittal','Coronal','Axial'][axis],fill='white')
        images.append(canvas)
    output=Path(config['output'])
    if output.exists():raise FileExistsError('Fresh blink version required')
    output.mkdir(parents=True)
    images[0].save(output/'registration_blink.gif',save_all=True,append_images=images[1:],duration=round(seconds*1000),loop=0)
    frames=[]
    for image in images:
        buffer=io.BytesIO();image.save(buffer,format='PNG');frames.append('data:image/png;base64,'+base64.b64encode(buffer.getvalue()).decode())
    page='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><style>body{font:18px/1.8 system-ui;margin:25px;background:#edf3f6}img{max-width:100%}button,input{font:inherit}</style><h1>MNI配准闪烁核查</h1><p>交替显示同网格模板与已归一化T1。先看脑室、胼胝体和深部结构，再看中央沟与病灶邻近；闪烁不是配准准确性评分。</p><button id="toggle">暂停闪烁</button><button id="reference">仅模板</button><button id="moving">仅患者</button><label>间隔(ms) <input id="interval" type="number" min="200" max="3000" step="100" value="'''+str(round(seconds*1000))+'''"></label><div id="label"></div><img id="frame"><details><summary>查看/保存独立GIF</summary><img src="registration_blink.gif"></details><p>两幅影像分别按各自1–99%灰度显示，空间网格未改变。切面为指定坐标最近体素的原网格平面；未进行新的配准、分割或临床验收。空间身份与审核绑定见 blink_receipt.json。</p><script>
const frames='''+json.dumps(frames)+''';let i=0,timer=null;const frame=document.getElementById('frame'),label=document.getElementById('label');
function show(j){i=j;frame.src=frames[i];label.textContent=i?'已归一化T1':'标准模板';}
function stop(){if(timer!==null)clearInterval(timer);timer=null;document.getElementById('toggle').textContent='继续闪烁';}
function start(){stop();let ms=Number(document.getElementById('interval').value);if(!Number.isFinite(ms)||ms<200||ms>3000)return;timer=setInterval(()=>show(1-i),ms);document.getElementById('toggle').textContent='暂停闪烁';}
document.getElementById('toggle').onclick=()=>timer===null?start():stop();document.getElementById('reference').onclick=()=>{stop();show(0);};document.getElementById('moving').onclick=()=>{stop();show(1);};document.getElementById('interval').onchange=start;show(0);start();
</script></html>'''
    page=page.replace('<h1>MNI配准闪烁核查</h1>','<h1>MNI配准闪烁核查</h1><p>空间：'+html.escape(config['space'])+'；范围：'+html.escape(config.get('scope','REAL_SAVED_IMAGES_VISUAL_REVIEW_PENDING'))+'</p>')
    (output/'blink.html').write_text(page,encoding='utf-8')
    save_new(output/'blink_receipt.json',dict(schema='neo-registration-blink-1',scope=config.get('scope','REAL_SAVED_IMAGES_VISUAL_REVIEW_PENDING'),
        space=config['space'],requested_coordinates_world_mm=point.tolist(),displayed_coordinates_world_mm=nib.affines.apply_affine(b.affine,index).tolist(),
        windows=window,anatomy_review=file_record(config['anatomy_review']),normalized=file_record(config['normalized_t1']),template=file_record(config['template']),
        provenance=file_record(config['normalization_provenance']),registration_fits=0,clinical_validation=False,
        display_scaling='Screen fit only, physical pixel aspect preserved; no volume resampling'))
    return str(output/'blink.html')

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',required=True);a=p.parse_args(argv);print(render(read(a.config)));return 0
