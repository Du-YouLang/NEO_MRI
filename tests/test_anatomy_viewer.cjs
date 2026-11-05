/* Local DOM/canvas simulation, not a real-browser or anatomical acceptance. */
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const crypto = require('node:crypto');

class Element {
  constructor(tag='div') { this.tag=tag;this.children=[];this.style={};this.value='0';this.checked=false;this.textContent=''; }
  append(...children) { this.children.push(...children); }
  getBoundingClientRect() { return {left:0,top:0,width:this.width,height:this.height}; }
  getContext() { return {createImageData:(w,h)=>({data:new Uint8ClampedArray(w*h*4)}),putImageData:im=>{this.pixels=im.data;}}; }
}
function digest(data) { return crypto.createHash('sha256').update(data).digest('hex'); }

async function check(file) {
  const html=fs.readFileSync(file,'utf8');
  const scripts=[...html.matchAll(/<script(?:\s[^>]*)?>([\s\S]*?)<\/script>/g)].map(m=>m[1]);
  const elements={};
  const get=id=>elements[id]||(elements[id]=new Element());
  const document={getElementById:get,createElement:tag=>new Element(tag)};
  const context=vm.createContext({Uint8Array,Uint8ClampedArray,Blob,Response,DecompressionStream,atob,console,
    setInterval:callback=>{context.tick=callback;return 1;},clearInterval:()=>{context.tick=null;},module:{exports:{}}});
  let checks=[];
  if(html.includes('id="payload"')) {
    get('payload').textContent=scripts[0];get('mode').value='full';get('outline').checked=true;get('alpha').value='.35';
    const selected=html.match(/<option value="(\d+)" selected>/);
    get('candidate').value=selected?selected[1]:'0';
    vm.runInContext(scripts.at(-1),context);context.document=document;
    await vm.runInContext('initialize()',context);
    const panels=get('panels').children;
    assert.equal(panels.length,3);assert.match(get('message').textContent,/RAS/);
    const canvases=panels.map(p=>p.children[1]);
    const snap=()=>canvases.map(c=>digest(c.pixels)).join('');
    const first=snap();
    for(const c of canvases)assert.equal(c.pixels.length,c.width*c.height*4);
    get('mode').value='keep';get('mode').onchange();const keep=snap();
    get('mode').value='drop';get('mode').onchange();assert.notEqual(snap(),keep);
    get('mode').value='full';get('mode').onchange();
    get('outline').checked=false;get('outline').onchange();get('shade').checked=true;
    get('alpha').value='0';get('alpha').oninput();const alpha0=snap();
    get('alpha').value='1';get('alpha').oninput();assert.notEqual(snap(),alpha0);
    get('shade').checked=false;get('outline').checked=true;get('outline').onchange();
    const payload=JSON.parse(scripts[0]);const choices=[];
    for(let i=0;i<payload.masks.length;i++) {get('candidate').value=String(i);get('candidate').onchange();choices.push(snap());}
    for(const p of panels){const slider=p.children[2];slider.value='1';slider.oninput();assert.match(p.children[0].textContent,/1 \/ /);}
    canvases[2].onwheel({deltaY:1,preventDefault(){}});assert.equal(Number(panels[2].children[2].value),2);
    canvases[2].onclick({clientX:canvases[2].width/2,clientY:canvases[2].height/2});
    assert.equal(Number(panels[0].children[2].value),Math.floor(canvases[2].width/2));
    get('mode').value='diff';get('mode').onchange();assert.ok(canvases.every(c=>c.pixels.length>0));
    checks=['gzip payload decoded','3 orthogonal canvases','native RAS coordinate text','keep/drop modes',
      'live alpha input','all candidate switches','slice sliders','wheel','linked click','difference rendering'];
    return {file,scope:'DOM_SIMULATION_NOT_REAL_BROWSER',checks,candidates:choices.length,distinct_center_views:new Set(choices).size};
  }
  context.document=document;get('interval').value='700';vm.runInContext(scripts.at(-1),context);
  const first=get('frame').src;assert.ok(first.startsWith('data:image/png;base64,'));assert.ok(context.tick);
  context.tick();assert.notEqual(get('frame').src,first);
  get('toggle').onclick();assert.equal(context.tick,null);
  get('reference').onclick();assert.equal(get('frame').src,first);
  get('moving').onclick();assert.notEqual(get('frame').src,first);
  get('interval').value='200';get('interval').onchange();assert.ok(context.tick);
  get('interval').value='199';get('interval').onchange();assert.equal(context.tick,null);
  checks=['two embedded PNG frames','timer alternates','pause','template only','subject only','interval bounds'];
  return {file,scope:'DOM_SIMULATION_NOT_REAL_BROWSER',checks};
}

(async()=>{
  if(process.argv.length<3)throw Error('Pass generated local review.html or blink.html paths');
  const results=[];for(const file of process.argv.slice(2))results.push(await check(file));
  console.log(JSON.stringify({passed:true,results},null,2));
})().catch(error=>{console.error(error);process.exitCode=1;});
