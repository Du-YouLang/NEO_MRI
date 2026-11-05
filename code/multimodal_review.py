"""Same-grid structural comparison viewer; quantizes display only, never registers."""
from pathlib import Path
import json,html
import nibabel as nib
import numpy as np
from anatomy_review import pack

SCRIPT=r"""
async function unpack(x){return new Uint8Array(await new Response(new Blob([Uint8Array.from(atob(x),c=>c.charCodeAt(0))]).stream().pipeThrough(new DecompressionStream('gzip'))).arrayBuffer());}
const get=id=>document.getElementById(id);let spaces,space,frame,canvases=[],sliders=[],labels=[],cursor,phase=false,timer;
function idx(p){return (p[0]*space.shape[1]+p[1])*space.shape[2]+p[2];}
function point(axis,u,v,k){return axis===0?[k,u,space.shape[2]-1-v]:axis===1?[u,k,space.shape[2]-1-v]:[u,space.shape[1]-1-v,k];}
function draw(){
 frame=space.frames[+get('frame').value];const mode=get('mode').value,alpha=+get('alpha').value;
 for(let axis=0;axis<3;axis++){
  const c=canvases[axis],ctx=c.getContext('2d'),im=ctx.createImageData(c.width,c.height);
  for(let v=0;v<c.height;v++)for(let u=0;u<c.width;u++){
   const p=point(axis,u,v,cursor[axis]),i=idx(p),a=space.background[i],b=frame.image[i];
   let gray=mode==='fixed'?a:mode==='moving'?b:mode==='checker'?((Math.floor(u/12)+Math.floor(v/12))%2?a:b):mode==='blink'?(phase?b:a):(1-alpha)*a+alpha*b;
   const k=(v*c.width+u)*4;im.data[k]=im.data[k+1]=im.data[k+2]=gray;im.data[k+3]=255;
   if(get('lesion').checked&&frame.lesion&&frame.lesion[i]){im.data[k]=255;im.data[k+1]=Math.round(gray*.3);im.data[k+2]=Math.round(gray*.3);}
  }
  ctx.putImageData(im,0,0);sliders[axis].value=cursor[axis];labels[axis].textContent=['矢状','冠状','横断'][axis]+' '+cursor[axis]+' / '+(space.shape[axis]-1);
 }
 const p=[...cursor,1],world=space.affine.slice(0,3).map(row=>row.reduce((s,x,j)=>s+x*p[j],0));
 get('coord').textContent=space.name+' · RAS毫米 '+world.map(x=>x.toFixed(1)).join(', ')+' · '+frame.name;
}
function setSpace(){
 space=spaces[+get('space').value];cursor=space.shape.map(n=>Math.floor(n/2));get('frame').innerHTML=space.frames.map((x,i)=>'<option value="'+i+'">'+x.name+'</option>').join('');get('frame').value='0';get('panels').innerHTML='';canvases=[];sliders=[];labels=[];
 for(let axis=0;axis<3;axis++){
  const panel=document.createElement('div'),label=document.createElement('div'),c=document.createElement('canvas'),slider=document.createElement('input');panel.className='panel';
  const h=axis===0?1:0,v=axis===2?1:2;c.width=space.shape[h];c.height=space.shape[v];c.style.aspectRatio=String(c.width*space.spacing[h]/(c.height*space.spacing[v]));slider.type='range';slider.min=0;slider.max=space.shape[axis]-1;
  slider.oninput=()=>{cursor[axis]=+slider.value;draw();};c.onwheel=e=>{e.preventDefault();cursor[axis]=Math.max(0,Math.min(space.shape[axis]-1,cursor[axis]+(e.deltaY>0?1:-1)));draw();};
  c.onclick=e=>{const r=c.getBoundingClientRect(),u=Math.min(c.width-1,Math.max(0,Math.floor((e.clientX-r.left)/r.width*c.width))),v0=Math.min(c.height-1,Math.max(0,Math.floor((e.clientY-r.top)/r.height*c.height)));cursor=point(axis,u,v0,cursor[axis]);draw();};
  panel.append(label,c,slider);get('panels').append(panel);canvases.push(c);sliders.push(slider);labels.push(label);
 }
 draw();
}
async function initialize(){
 spaces=JSON.parse(get('payload').textContent);
 for(const p of spaces){p.background=await unpack(p.background);for(const f of p.frames){f.image=await unpack(f.image);if(f.lesion)f.lesion=await unpack(f.lesion);const n=p.shape.reduce((a,b)=>a*b,1);if(f.image.length!==n||(f.lesion&&f.lesion.length!==n))throw Error('影像长度不匹配');}}
 get('space').innerHTML=spaces.map((x,i)=>'<option value="'+i+'">'+x.name+'</option>').join('');get('space').value='0';
 get('space').onchange=setSpace;get('frame').onchange=draw;get('mode').onchange=draw;get('alpha').oninput=draw;get('lesion').onchange=draw;
 timer=setInterval(()=>{if(get('mode').value==='blink'){phase=!phase;draw();}},700);setSpace();
}
initialize().catch(e=>get('coord').textContent=e.message);
"""

def export(spaces,output,title):
    """spaces: name, background, frames[name,image,optional lesion]. All exact-grid."""
    payload=[]
    for item in spaces:
        bg=nib.as_closest_canonical(nib.load(item['background']));b=bg.get_fdata(dtype=np.float32)
        positive=b[b>0];lo,hi=np.percentile(positive,[1,99])
        def display(data,window):
            if not np.isfinite(data).all() or window[1]<=window[0]:raise ValueError('Finite image and valid window required')
            return pack(np.clip((data-window[0])/(window[1]-window[0])*255,0,255).astype('uint8'))
        frames=[];window=None
        for entry in item['frames']:
            im=nib.as_closest_canonical(nib.load(entry['image']))
            if im.shape!=bg.shape or not np.allclose(im.affine,bg.affine,atol=1e-5,rtol=0):raise ValueError('Viewer requires exact same grid, no implicit registration')
            data=im.get_fdata(dtype=np.float32)
            if window is None:window=np.percentile(data[data>0],[1,99])
            f=dict(name=entry['name'],image=display(data,window))
            if entry.get('lesion'):
                lm=nib.as_closest_canonical(nib.load(entry['lesion']))
                if lm.shape!=bg.shape or not np.allclose(lm.affine,bg.affine,atol=1e-5,rtol=0):raise ValueError('Lesion grid mismatch')
                f['lesion']=pack((lm.get_fdata()>.5).astype('uint8'))
            frames.append(f)
        payload.append(dict(name=item['name'],shape=list(bg.shape),spacing=list(map(float,bg.header.get_zooms())),affine=bg.affine.tolist(),background=display(b,[lo,hi]),frames=frames))
    page='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><style>body{font:17px/1.7 system-ui;margin:20px;background:#eef4f8;color:#183c50}select{font:inherit;max-width:95%}#panels{display:flex;gap:15px;flex-wrap:wrap}.panel{flex:1;min-width:280px}canvas{width:100%;background:black}input[type=range]{width:95%}.note{background:white;padding:18px;border-left:5px solid #d9a229}</style><h1>'''+html.escape(title)+'''</h1><p class="note">待人工审核。全空间网格保留，仅显示灰度量化为8位。T1背景与T2/FLAIR采用不同灰度窗口；同一空间的候选共用功能下拉框第一幅结构图的灰度窗口，亮度差异不等于位置误差。红色为病灶。配准内部金字塔有尺度平滑，最终输出没有额外平滑。此页不拟合或修改影像。</p><label>空间 <select id="space"></select></label><br><label>结构候选 <select id="frame"></select></label> <select id="mode"><option value="blend">透明叠加</option><option value="checker">棋盘格</option><option value="blink">交替闪烁</option><option value="fixed">只看背景</option><option value="moving">只看候选</option></select><label><input type="checkbox" id="lesion" checked>病灶</label><label>候选透明度<input id="alpha" type="range" min="0" max="1" step="0.05" value="0.5"></label><div id="coord">正在解压本地图像…</div><div id="panels"></div><p>支持逐层滑条、滚轮与点击三面联动。RAS坐标是当前空间的物理坐标；原生空间不冒充MNI。图像来自同一访视；指标只作描述，不自动批准临床使用。</p><script id="payload" type="application/json">'''+json.dumps(payload,separators=(',',':')).replace('<','\\u003c')+'</script><script>'+SCRIPT+'</script></html>'
    p=Path(output)
    if p.exists():raise FileExistsError('New viewer version required')
    p.parent.mkdir(parents=True,exist_ok=True);p.write_text(page,encoding='utf-8');return str(p)
