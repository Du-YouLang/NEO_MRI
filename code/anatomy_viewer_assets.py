"""Native T1 multi-mask viewer; no hidden processing."""
SCRIPT = r"""
// Original native grid, RAS-oriented axes. No resampling, no fitting.
async function inflate(value){
 const bytes=Uint8Array.from(atob(value),c=>c.charCodeAt(0));
 const stream=new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'));
 return new Uint8Array(await new Response(stream).arrayBuffer());
}
function flatIndex(x,y,z,shape){return (x*shape[1]+y)*shape[2]+z;}
function planePoint(axis,u,v,k,shape){
 if(axis===0)return [k,u,shape[2]-1-v];
 if(axis===1)return [u,k,shape[2]-1-v];
 return [u,shape[1]-1-v,k];
}
async function initialize(){
 const p=JSON.parse(document.getElementById('payload').textContent);
 const image=await inflate(p.image);const masks=await Promise.all(p.masks.map(m=>inflate(m.mask)));let mask=masks[0];
 const shape=p.shape;const n=shape.reduce((a,b)=>a*b,1);
 if(image.length!==n||masks.some(m=>m.length!==n))throw Error('影像长度核对失败');
 const cursor=p.crosshair.slice(),canvases=[],sliders=[],labels=[];
 for(let axis=0;axis<3;axis++){
  const panel=document.createElement('div');panel.className='panel';
  const label=document.createElement('div'),canvas=document.createElement('canvas'),slider=document.createElement('input');
  slider.type='range';slider.min=0;slider.max=shape[axis]-1;slider.value=cursor[axis];
  const horizontal=axis===0?1:0,vertical=axis===2?1:2;
  canvas.width=shape[horizontal];canvas.height=shape[vertical];
  canvas.style.aspectRatio=String(canvas.width*p.spacing[horizontal]/(canvas.height*p.spacing[vertical]));
  slider.oninput=()=>{cursor[axis]=+slider.value;render();};
  canvas.onwheel=event=>{event.preventDefault();cursor[axis]=Math.max(0,Math.min(shape[axis]-1,cursor[axis]+(event.deltaY>0?1:-1)));render();};
  canvas.onclick=event=>{
   const rect=canvas.getBoundingClientRect();
   const u=Math.min(canvas.width-1,Math.max(0,Math.floor((event.clientX-rect.left)/rect.width*canvas.width)));
   const v=Math.min(canvas.height-1,Math.max(0,Math.floor((event.clientY-rect.top)/rect.height*canvas.height)));
   const xyz=planePoint(axis,u,v,cursor[axis],shape);xyz.forEach((q,j)=>cursor[j]=q);render();
  };
  panel.append(label,canvas,slider);document.getElementById('panels').append(panel);
  canvases.push(canvas);sliders.push(slider);labels.push(label);
 }
 function render(){
  mask=masks[Number(document.getElementById('candidate').value)];const mode=document.getElementById('mode').value,outline=document.getElementById('outline').checked,shade=document.getElementById('shade').checked;
  const world=p.affine.slice(0,3).map(row=>row[3]+row.slice(0,3).reduce((sum,a,j)=>sum+a*cursor[j],0));
  document.getElementById('message').textContent='原生索引 '+cursor.join(', ')+' | RAS物理坐标(mm) '+world.map(x=>x.toFixed(2)).join(', ')+' | 非MNI。屏幕轴向：左L→右R（矢状为后P→前A），上S/前A。';
  for(let axis=0;axis<3;axis++){
   const canvas=canvases[axis],context=canvas.getContext('2d');
   const out=context.createImageData(canvas.width,canvas.height);
   for(let v=0;v<canvas.height;v++)for(let u=0;u<canvas.width;u++){
    const xyz=planePoint(axis,u,v,cursor[axis],shape),idx=flatIndex(...xyz,shape),inside=mask[idx]!==0;
    const gray=(mode==='keep'&&!inside)||(mode==='drop'&&inside)?0:image[idx];
    let red=gray,green=gray,blue=gray,boundary=false;
    if(inside&&outline){
     for(const [a,b] of [[u-1,v],[u+1,v],[u,v-1],[u,v+1]]){
      if(a<0||b<0||a>=canvas.width||b>=canvas.height){boundary=true;break;}
      const q=planePoint(axis,a,b,cursor[axis],shape);
      if(mask[flatIndex(...q,shape)]===0){boundary=true;break;}
     }
    }
    if(mode==='diff' && mask[idx]!==masks[0][idx]){red=inside?30:240;green=inside?240:30;blue=inside?70:240;}
    if(inside&&shade){const alpha=Number(document.getElementById('alpha').value);red=Math.round((1-alpha)*gray+alpha*255);green=Math.round((1-alpha)*gray);blue=green;}
    if(boundary){red=255;green=45;blue=45;}
    const offset=(v*canvas.width+u)*4;out.data.set([red,green,blue,255],offset);
   }
   context.putImageData(out,0,0);sliders[axis].value=cursor[axis];
   labels[axis].textContent=['矢状（左/右位置）','冠状（前/后位置）','轴向（上/下位置）'][axis]+'：'+cursor[axis]+' / '+(shape[axis]-1);
  }
 }
 for(const id of ['mode','outline','shade','candidate','alpha'])document.getElementById(id).onchange=render;
 document.getElementById('alpha').oninput=render;
 render();
}
if(typeof document!=='undefined')initialize().catch(e=>document.getElementById('message').textContent='查看器启动失败：'+e.message+'。请使用新版Chrome/Edge，或先查看静态总览。');
if(typeof module!=='undefined')module.exports={flatIndex,planePoint};

"""
