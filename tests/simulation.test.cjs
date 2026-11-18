const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const { ODESolvers, MethodsManifest } = require('../assets/engine.js');

test('Heun acoplado reproduz o exemplo numérico documentado no arquivo Octave', () => {
  const data = ODESolvers.solveRK2Sistema();
  assert.equal(data.steps.length, 12);
  const first = data.steps[0];
  assert.ok(Math.abs(first.yNext - 0.458) < 1e-12);
  assert.ok(Math.abs(first.zNext - 1.10725) < 1e-12);
  assert.ok(Math.abs(data.steps.at(-1).xNext - 1.2) < 1e-12);
});

test('Heun acoplado converge com ordem dois para uma solução analítica conhecida', () => {
  const solve = h => ODESolvers.solveRK2Sistema({h,xf:1,y0:1,z0:1,
    f:(_x,y)=>-y,g:(_x,_y,z)=>-2*z}).steps.at(-1);
  const coarse = solve(0.1), fine = solve(0.05);
  for (const [field,exact] of [['yNext',Math.exp(-1)],['zNext',Math.exp(-2)]]) {
    const ratio = Math.abs(coarse[field]-exact)/Math.abs(fine[field]-exact);
    assert.ok(ratio > 3.5 && ratio < 5, `${field}: razão dos erros ${ratio}`);
  }
});

test('Todos os métodos animados alimentam valores e colunas coerentes no painel', () => {
  const context = {window:{}};
  vm.runInNewContext(fs.readFileSync(require.resolve('../assets/simulation-panel.js'),'utf8'), context);
  for (const [key,manifest] of Object.entries(MethodsManifest)) {
    const data = manifest.solver();
    for (const index of [0,data.steps.length-1]) {
      const view = context.window.simulationPresentation({methodKey:key,manifest,methodData:data,
        stepData:data.steps[index],currentStep:index,totalSteps:data.steps.length});
      assert.equal(view.fields.length,3,key);
      assert.equal(view.rows.length,data.steps.length,key);
      for (const row of view.rows) assert.equal(row.length,view.headers.length,key);
      assert.ok(view.status,key);
    }
  }
});
