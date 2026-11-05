/* DOM/canvas simulation only. Does not claim a real-browser or anatomical review. */
const fs=require('node:fs'),vm=require('node:vm'),assert=require('node:assert/strict'),crypto=require('node:crypto');
class E{
 constructor(){this.children=[];this.style={};this.value='0';this.checked=true;this.textContent='';}
 append(...xs){this.children.push(...xs);}
 set innerHTML(x){this.html=x;this.children=[];}get innerHTML(){return this.html;}
 getBoundingClientRect(){return {left:0,top:0,width:this.width,height:this.height};}
 getContext(){return {createImageData:(w,h)=>({data:new Uint8ClampedArray(w*h*4)}),putImageData:x=>{this.pixels=x.data;}};}
}
async function check(path){
 const page=fs.readFileSync(path,'utf8'),scripts=[...page.matchAll(/<script(?:\s[^>]*)?>([\s\S]*?)<\/script>/g)].map(x=>x[1]);
 const es={},get=id=>es[id]||(es[id]=new E());get('payload').textContent=scripts[0];get('alpha').value='.5';get('mode').value='blend';
 const ctx=vm.createContext({document:{getElementById:get,createElement:()=>new E()},Blob,Response,DecompressionStream,atob,Uint8Array,Uint8ClampedArray,console,setInterval:cb=>{ctx.tick=cb;return 1;}});
 vm.runInContext(scripts.at(-1).replace("initialize().catch(e=>get('coord').textContent=e.message);",'globalThis.ready=initialize();'),ctx);await ctx.ready;
 const payload=JSON.parse(scripts[0]);let count=0;
 for(let i=0;i<payload.length;i++){
  get('space').value=String(i);get('space').onchange();assert.equal(get('panels').children.length,3);
  const canvases=get('panels').children.map(p=>p.children[1]);
  for(const c of canvases)assert.equal(c.pixels.length,c.width*c.height*4);
  assert.match(get('coord').textContent,/RAS/);
  get('mode').value='fixed';get('mode').onchange();const digest=()=>crypto.createHash('sha256').update(canvases[2].pixels).digest('hex');const first=digest();
  get('mode').value='moving';get('mode').onchange();
  for(let j=0;j<payload[i].frames.length;j++){get('frame').value=String(j);get('frame').onchange();count++;}
  get('mode').value='checker';get('mode').onchange();get('alpha').value='.8';get('alpha').oninput();get('lesion').checked=false;get('lesion').onchange();
  get('mode').value='blink';get('mode').onchange();ctx.tick();assert.ok(digest().length);
  const p=get('panels').children[2];p.children[2].value='1';p.children[2].oninput();assert.match(p.children[0].textContent,/1 \/ /);
  canvases[2].onwheel({deltaY:1,preventDefault(){}});assert.equal(+p.children[2].value,2);
  canvases[2].onclick({clientX:canvases[2].width/2,clientY:canvases[2].height/2});
 }
 return {file:path,spaces:payload.length,frames_checked:count,scope:'DOM_SIMULATION_ONLY',passed:true};
}
Promise.all(process.argv.slice(2).map(check)).then(x=>console.log(JSON.stringify(x,null,2))).catch(e=>{console.error(e);process.exitCode=1;});
