// Exercise the offline report's controls against a synthetic HTML result.
const assert = require("node:assert/strict");
const fs = require("node:fs");
const vm = require("node:vm");

class Element {
  constructor() {
    this.listeners = {};
    this.dataset = {};
    this.style = {};
    this.children = [];
    this.checked = true;
    this.value = "";
  }
  addEventListener(name, callback) { this.listeners[name] = callback; }
  getBoundingClientRect() { return {left: 0, top: 0, width: 640, height: 440}; }
  getContext() {
    return new Proxy({}, {get: (target, key) => target[key] ??
      (key === "measureText" ? () => ({width: 30}) : () => {})});
  }
  append(...children) { this.children.push(...children); }
  replaceChildren(...children) { this.children = children; }
  setPointerCapture() {}
}

const source = fs.readFileSync(process.argv[2], "utf8");
const scripts = [...source.matchAll(/<script(?:\s[^>]*)?>([\s\S]*?)<\/script>/g)].map(match => match[1]);
const data = JSON.parse(scripts[0]);
const expectedOrigin = {SMATT: "M1", HCPA: "CST", Boyne2025: "Whole"}[data.method];
assert.deepEqual([...data.origins].sort(), ["Empty", expectedOrigin].sort());
const elements = new Map();
const by = id => {
  if (!elements.has(id)) elements.set(id, new Element());
  return elements.get(id);
};
by("smatt-data").textContent = scripts[0];
by("side").value = "Left";
const pending = [];
const document = {
  documentElement: new Element(),
  getElementById: by,
  createElement: () => new Element(),
  createTextNode: value => value,
  querySelectorAll: () => [],
};
const globals = {
  document, window: {devicePixelRatio: 1}, atob, Uint16Array,
  Image: class { constructor() { this.naturalWidth = 24; } },
  ResizeObserver: class { observe() {} },
  requestAnimationFrame: callback => pending.push(callback),
};
const flush = () => { while (pending.length) pending.shift()(); };
const assertBars = (panel, value) => {
  const bars = JSON.parse(by(panel).dataset.values);
  assert.deepEqual(bars, data.origins.map(origin => origin === expectedOrigin ? value : null));
};
vm.runInNewContext(scripts[1], globals);
flush();
assert.equal(document.documentElement.dataset.ready, "true");
assert.equal(by("C").dataset.z, 2);
assertBars("C", 37.5);
assertBars("D", .625);
by("z").value = "1";
by("z").listeners.input();
flush();
assert.equal(by("C").dataset.z, 1);
assertBars("C", 20);
assert.match(by("axial-title").textContent, /\+1 mm/);
by("F").listeners.click({clientX: Number(by("F").dataset.right)});
flush();
assert.equal(by("D").dataset.z, 4);
by("side").value = "Right";
by("side").listeners.change();
flush();
assertBars("E", 0);
by("scene").listeners.keydown({key: "ArrowRight", preventDefault() {}});
flush();
console.log("Offline viewer controls passed DOM/canvas simulation");
