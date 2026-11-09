"use strict";
const $=id=>document.getElementById(id);
const S=window.BOLD_AB, [nx,ny,nz]=S.shape, nk=S.params.length;
const planes=["轴位","冠状位","矢状位"], sizes=[[nx,ny],[nx,nz],[ny,nz]];
function decode(s,Type){
 const raw=atob(s), bytes=new Uint8Array(raw.length);
 for(let i=0;i<raw.length;i++)bytes[i]=raw.charCodeAt(i);
 return new Type(bytes.buffer);
}
const anat=decode(S.anatomy,Uint8Array), vox=decode(S.voxel,Int32Array), shell=decode(S.shell,Int32Array);
const lookup=new Int32Array(anat.length).fill(-1);vox.forEach((v,i)=>lookup[v]=i);
const beta=decode(S.beta,Float32Array), tvals=decode(S.t,Float32Array);
const bits=decode(S.bits,Uint8Array), reasons=decode(S.reason,Uint8Array);
const g=decode(S.g,Float32Array), box=decode(S.box,Float32Array), winner=decode(S.choice,Uint16Array);
const groups=S.groups.map(x=>({...x,vox:decode(x.points,Int32Array)}));
const modeMask={A:1,A_stable:2,B:4,B_stable:8};
const modes=["A","A_stable","B","B_stable"];
const query=new URLSearchParams(location.search), userMode=query.get("mode");
let userGroup=query.get("group");
if(modes.includes(userMode))$("mode").value=userMode;
$("title").textContent=S.run+" · "+S.task;
$("time").max=S.frame.length-1;
const center=[Math.floor(nx/2),Math.floor(ny/2),Math.floor(nz/2)];
const dyn=[...center], staticCursor=[...center];
const shape=S.shape, clip=(n,min,max)=>Math.min(max,Math.max(min,n));
const at=(x,y,z)=>(x*ny+y)*nz+z;
function xyz(f){
 const z=f%nz, q=Math.floor(f/nz);
 return [Math.floor(q/ny),q%ny,z];
}
function world(p){
 return S.affine.slice(0,3).map(a=>a[0]*p[0]+a[1]*p[1]+a[2]*p[2]+a[3]);
}
function coords(k,h,v,cur){
 return k===0?[h,ny-1-v,cur[2]]:k===1?[h,cur[1],nz-1-v]:[cur[0],h,nz-1-v];
}
function where(k,p){
 return k===0?[p[0],ny-1-p[1]]:k===1?[p[0],nz-1-p[2]]:[p[1],nz-1-p[2]];
}
let hovered=null,group=null,groupVox=new Uint8Array(anat.length);
const az0=-.63,el0=.32;
let az=az0,el=el0,zoom=1,timer=null,hold=null;
function makePanel(root,cur,onChange,onHover){
 return sizes.map(([w,h],k)=>{
  const f=document.createElement("figure"), c=document.createElement("canvas"), caption=document.createElement("figcaption");
  c.width=w;c.height=h;c.dataset.plane=String(k);
  f.append(c,caption);root.append(f);
  const ctx=c.getContext("2d",{willReadFrequently:true}), pix=ctx.createImageData(w,h);
  const toWorld=ev=>{
   const r=c.getBoundingClientRect();
   return coords(k,clip(Math.floor((ev.clientX-r.left)*w/r.width),0,w-1),
                 clip(Math.floor((ev.clientY-r.top)*h/r.height),0,h-1),cur);
  };
  c.addEventListener("wheel",ev=>{
   ev.preventDefault();const i=2-k;
   cur[i]=clip(cur[i]+Math.sign(ev.deltaY),0,shape[i]-1);onChange();
  },{passive:false});
  c.addEventListener("click",ev=>{cur.splice(0,3,...toWorld(ev));onChange();});
  if(onHover)c.addEventListener("pointermove",ev=>onHover(toWorld(ev)));
  return {w,h,k,c,ctx,pix,caption};
 });
}
function renderSlices(panels,cur,field){
 for(const {w,h,k,c,ctx,pix,caption} of panels){
  for(let v=0;v<h;v++)for(let u=0;u<w;u++){
   const [x,y,z]=coords(k,u,v,cur), f=at(x,y,z), l=anat[f];
   const rgb=field(f,l);
   const j=4*(v*w+u);pix.data[j]=rgb[0];pix.data[j+1]=rgb[1];pix.data[j+2]=rgb[2];pix.data[j+3]=255;
  }
  ctx.putImageData(pix,0,0);
  const [u,v]=where(k,cur);
  ctx.strokeStyle="#b6d9c4";ctx.globalAlpha=.7;ctx.lineWidth=.65;
  ctx.beginPath();ctx.moveTo(u+.5,Math.max(0,v-8));ctx.lineTo(u+.5,Math.min(h,v+9));
  ctx.moveTo(Math.max(0,u-8),v+.5);ctx.lineTo(Math.min(w,u+9),v+.5);ctx.stroke();ctx.globalAlpha=1;
  caption.textContent=planes[k]+" · "+world(cur)[2-k].toFixed(1)+" mm";
  c.dataset.ready="true";
 }
}
const tint=(g,rgb,a)=>[g*(1-a)+rgb[0]*a,g*(1-a)+rgb[1]*a,g*(1-a)+rgb[2]*a];
const gray=v=>[v,v,v];
function visible(i,mode,branch){
 return i>=0 && (bits[i]&modeMask[mode])!==0 && (branch==="all"||(branch==="scan")===(reasons[i]===0));
}
function describe(){
 const p=hovered||dyn, i=lookup[at(...p)], mm=world(p).map(x=>x.toFixed(1)).join(", ");
 const mode=$("mode").value, branch=$("branch").value;
 let extra="";
 if(visible(i,mode,branch)){
  const [l,d]=S.params[winner[i]];
  extra=" · "+(reasons[i]?"回退到默认":"扫描保留")+" · lag "+l+" s · ratio "+d+
        " · β "+beta[i].toFixed(3)+" · t "+tvals[i].toFixed(2);
 }
 $("status").textContent="集合 "+S.counts[mode].toLocaleString()+" 体素 · MNI "+mm+" mm"+extra;
}
const panels=makePanel($("dynamic"),dyn,renderDynamic,p=>{hovered=p;describe();});
const staticPanels=makePanel($("group-slices"),staticCursor,renderStatic);
const signalColors=[null,[244,68,38],[38,120,246],[251,174,54],[27,192,208]];
const signalCode=new Uint8Array(vox.length), signalAlpha=new Float32Array(vox.length);
const volume=$("dynamic-orbit"), volumeCtx=volume.getContext("2d");
const volumeBackdrop=document.createElement("canvas");
volumeBackdrop.width=volume.width;volumeBackdrop.height=volume.height;
const backdropCtx=volumeBackdrop.getContext("2d");
const view3D={az:az0,el:el0,zoom:1};
const pointX=new Float32Array(vox.length), pointY=new Float32Array(vox.length), pointDepth=new Float32Array(vox.length);
const depthOrder=Uint32Array.from(vox,(_,i)=>i);
let needsProjection=true, volumeHold=null, pendingVolumeDraw=0;

function projectVolume(){
 const w=volume.width,h=volume.height;
 const ca=Math.cos(view3D.az),sa=Math.sin(view3D.az),ce=Math.cos(view3D.el),se=Math.sin(view3D.el);
 const size=Math.min(w,h)/Math.max(...shape)*1.05*view3D.zoom;
 function project(x,y,z){
  const u=ca*x-sa*y,d=sa*x+ca*y;
  return [u,ce*z-se*d,ce*d+se*z];
 }
 function projectVoxel(f){
  const [x,y,z]=xyz(f);
  return project(x-center[0],y-center[1],z-center[2]);
 }
 backdropCtx.fillStyle="#10191d";backdropCtx.fillRect(0,0,w,h);
 backdropCtx.fillStyle="rgba(190,210,212,.17)";
 for(const f of shell){
  const [x,y]=projectVoxel(f);
  backdropCtx.fillRect(w/2+size*x,h/2-size*y,1.6,1.6);
 }
 for(let i=0;i<vox.length;i++){
  const [x,y,d]=projectVoxel(vox[i]);
  pointX[i]=w/2+size*x;pointY[i]=h/2-size*y;pointDepth[i]=d;
 }
 depthOrder.sort((a,b)=>pointDepth[a]-pointDepth[b]);
 for(const [label,axis,color] of [["R",[1,0,0],"#f58f80"],["A",[0,1,0],"#91d0a2"],["S",[0,0,1],"#8ec6f0"]]){
  const [x,y]=project(...axis),ox=48,oy=h-45;
  backdropCtx.strokeStyle=color;backdropCtx.fillStyle=color;backdropCtx.font="12px sans-serif";
  backdropCtx.beginPath();backdropCtx.moveTo(ox,oy);backdropCtx.lineTo(ox+26*x,oy-26*y);backdropCtx.stroke();
  backdropCtx.fillText(label,ox+30*x-4,oy-30*y+4);
 }
 needsProjection=false;
}
function drawVolume(){
 if(needsProjection)projectVolume();
 volumeCtx.drawImage(volumeBackdrop,0,0);
 const size=clip(2.7*view3D.zoom,2,5);
 let last=0,n=0;
 for(const i of depthOrder){
  const code=signalCode[i];
  if(!code)continue;
  if(code!==last){volumeCtx.fillStyle="rgb("+signalColors[code].join(",")+")";last=code;}
  volumeCtx.globalAlpha=signalAlpha[i];
  volumeCtx.fillRect(pointX[i]-size/2,pointY[i]-size/2,size,size);
  n++;
 }
 volumeCtx.globalAlpha=1;
 const f=Number($("time").value);
 volume.dataset.ready="true";volume.dataset.frame=String(f);volume.dataset.colored=String(n);
 $("volume-label").textContent="三维 · "+S.frame[f]+" s · "+n.toLocaleString()+" 个着色体素";
}
function refreshVolume(){
 needsProjection=true;
 if(pendingVolumeDraw)return;
 pendingVolumeDraw=requestAnimationFrame(()=>{pendingVolumeDraw=0;drawVolume();});
}
volume.addEventListener("pointerdown",ev=>{
 volumeHold=[ev.clientX,ev.clientY];volume.setPointerCapture(ev.pointerId);
});
volume.addEventListener("pointermove",ev=>{
 if(!volumeHold)return;
 view3D.az+=(ev.clientX-volumeHold[0])*.013;
 view3D.el=clip(view3D.el+(ev.clientY-volumeHold[1])*.013,-1.5,1.5);
 volumeHold=[ev.clientX,ev.clientY];refreshVolume();
});
volume.addEventListener("pointerup",()=>{volumeHold=null;});
volume.addEventListener("lostpointercapture",()=>{volumeHold=null;});
volume.addEventListener("wheel",ev=>{
 ev.preventDefault();view3D.zoom=clip(view3D.zoom*Math.exp(-ev.deltaY*.001),.5,2.7);refreshVolume();
},{passive:false});
$("reset-dynamic").onclick=()=>{Object.assign(view3D,{az:az0,el:el0,zoom:1});refreshVolume();};

function renderDynamic(){
 const f=Number($("time").value), k=$("quantity").value, mode=$("mode").value, branch=$("branch").value;
 const scale=(k==="pre"||k==="beta"?S.pre_scale:k==="t"?8:S.post_scale)*Number($("contrast").value);
 $("time-label").textContent=S.frame[f]+" s · "+(f+1)+"/"+S.frame.length;
 $("scale").textContent=k==="branch"?"回退 / 扫描":"−"+scale.toPrecision(3)+" / 0 / +"+scale.toPrecision(3)+(k==="t"?" t":"");
 $("bar").style.background=k==="branch"?"linear-gradient(90deg,#f6b24b 50%,#24b4c8 50%)":"";
 signalCode.fill(0);
 for(let i=0;i<vox.length;i++){
  if(!visible(i,mode,branch))continue;
  if(k==="branch"){
   signalCode[i]=reasons[i]?3:4;signalAlpha[i]=.86;continue;
  }
  const value=k==="beta"?beta[i]:k==="t"?tvals[i]:(k==="pre"?box:g)[f*nk+winner[i]]*beta[i];
  if(Math.abs(value)<1e-8)continue;
  signalCode[i]=value>0?1:2;
  signalAlpha[i]=.84*Math.min(1,Math.abs(value)/scale);
 }
 renderSlices(panels,dyn,(v,l)=>{
  const i=lookup[v],code=i>=0?signalCode[i]:0;
  return code?tint(l,signalColors[code],signalAlpha[i]):gray(l);
 });
 drawVolume();
 hovered=null;describe();
}
function showGroup(){
 const index=Number($("group").value);
 group=groups[index]||null;groupVox.fill(0);
 const a=$("download");
 if(group){
  group.vox.forEach(x=>groupVox[x]=1);
  if(group.vox.length)staticCursor.splice(0,3,...xyz(group.vox[Math.floor(group.vox.length/2)]));
  a.href=group.mask_file;a.hidden=false;
  $("group-label").textContent=(group.landmark==="onset"?"共同起点":"共同终点")+" "+group.center+
    " s（相对任务起点） · "+group.n.toLocaleString()+" 个体素";
 }else{a.hidden=true;$("group-label").textContent="本集合没有达到预设大小与 Bootstrap 条件的时间群。";}
 renderStatic();
}
function renderStatic(){
 renderSlices(staticPanels,staticCursor,(v,l)=>groupVox[v]?tint(l,[255,160,20],.94):gray(l));
 orbit();
}
const cv=$("orbit"), cx=cv.getContext("2d");
function orbit(){
 const w=cv.width,h=cv.height;
 cx.fillStyle="#10191d";cx.fillRect(0,0,w,h);
 const ca=Math.cos(az), sa=Math.sin(az), ce=Math.cos(el), se=Math.sin(el);
 const size=Math.min(w,h)/Math.max(...shape)*.78*zoom;
 function point(f,color,r){
  const v=xyz(f), x=v[0]-center[0], y=v[1]-center[1], z=v[2]-center[2];
  const u=ca*x-sa*y, d=sa*x+ca*y, b=ce*z-se*d;
  cx.fillStyle=color;cx.fillRect(w/2+size*u,h/2-size*b,r,r);
 }
 for(const f of shell)point(f,"rgba(190,210,212,.19)",1.7);
 if(group)for(const f of group.vox)point(f,"rgba(82,233,157,.84)",2.6);
 cv.dataset.ready="true";
}
cv.addEventListener("pointerdown",e=>{hold=[e.clientX,e.clientY];cv.setPointerCapture(e.pointerId);});
cv.addEventListener("pointermove",e=>{
 if(!hold)return;az+=(e.clientX-hold[0])*.013;el=clip(el+(e.clientY-hold[1])*.013,-1.5,1.5);
 hold=[e.clientX,e.clientY];orbit();
});
cv.addEventListener("lostpointercapture",()=>{hold=null;});
cv.addEventListener("pointerup",()=>{hold=null;});
cv.addEventListener("wheel",e=>{e.preventDefault();zoom=clip(zoom*Math.exp(-e.deltaY*.001),.5,2.7);orbit();},{passive:false});
$("reset").onclick=()=>{az=az0;el=el0;zoom=1;orbit();};
for(const id of ["time","quantity","branch","contrast"])$(id).addEventListener("input",renderDynamic);
function fillGroups(){
 const key=$("mode").value.startsWith("A")?"A":"B", select=$("group");
 select.replaceChildren();
 groups.forEach((v,i)=>{
  if(v.mode===key)select.add(new Option((v.landmark==="onset"?"起点":"终点")+" "+v.center+" s · "+v.n+" 体素",i));
 });
 if(!select.length)select.add(new Option("没有满足条件的群", "-1"));
 const chosen=groups.findIndex(g=>g.mode===key&&g.name===userGroup);
 if(chosen>=0)select.value=chosen;
 userGroup=null;
 showGroup();
}
$("mode").addEventListener("change",()=>{renderDynamic();fillGroups();});
$("group").addEventListener("change",showGroup);
$("play").onclick=()=>{
 if(timer!==null){clearInterval(timer);timer=null;$("play").textContent="播放";return;}
 $("play").textContent="暂停";
 timer=setInterval(()=>{$("time").value=(Number($("time").value)+1)%S.frame.length;renderDynamic();},145);
};
renderDynamic();fillGroups();
