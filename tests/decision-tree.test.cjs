const test = require('node:test');
const assert = require('node:assert/strict');
const { DecisionTreeComponent, ALGORITHMS_DATABASE, CATEGORIES } = require('../assets/decision_tree.js');
const { getMethodKeyByFile } = require('../assets/engine.js');

function model() {
  const element = attrs => {
    const states = new Set();
    return { getAttribute:key=>attrs[key],setAttribute:(key,value)=>{attrs[key]=value;},
      classList:{toggle:(key,on)=>on?states.add(key):states.delete(key),contains:key=>states.has(key)} };
  };
  const cards = ALGORITHMS_DATABASE.map(a=>element({'data-slug':a.slug,'data-category':a.category}));
  const headers = CATEGORIES.map(c=>element({'data-category':c.id}));
  const branches = CATEGORIES.map(c=>element({'data-category':c.id}));
  const columns = CATEGORIES.map((c,i)=>Object.assign(element({'data-category':c.id}),{querySelector:()=>headers[i]}));
  const selectors = {'.dt-algorithm-card':cards,'.cat-header-group':headers,'.bus-branch':branches,'.dt-category-column':columns};
  const tree = Object.create(DecisionTreeComponent.prototype);
  Object.assign(tree,{currentCategory:'all',hoveredCategory:null,currentFilter:'',container:{querySelectorAll:s=>selectors[s]}});
  const visibleCategories = () => [...new Set(cards.filter(c=>!c.classList.contains('dimmed')).map(c=>c.getAttribute('data-category')))];
  return {tree,visibleCategories,cards,headers};
}

test('prévia por coluna é temporária, clique fixa, troca, alterna e núcleo limpa', () => {
  const {tree,visibleCategories} = model();
  tree.previewCategory('raizes'); assert.deepEqual(visibleCategories(),['raizes']);
  tree.previewCategory(null); assert.equal(visibleCategories().length,6);
  tree.toggleCategory('raizes'); tree.previewCategory(null); assert.deepEqual(visibleCategories(),['raizes']);
  tree.previewCategory('edos'); assert.deepEqual(visibleCategories(),['raizes']);
  tree.previewCategory(null); assert.deepEqual(visibleCategories(),['raizes']);
  tree.toggleCategory('calculo'); assert.deepEqual(visibleCategories(),['calculo']);
  tree.previewCategory('sistemas'); assert.deepEqual(visibleCategories(),['calculo']);
  tree.toggleCategory('calculo'); assert.equal(visibleCategories().length,6);
  tree.toggleCategory('sistemas'); tree.clearCategory(); assert.equal(visibleCategories().length,6);
});

test('busca e destaque por coluna funcionam juntos', () => {
  const {tree,cards} = model();
  tree.setFilter('bissecao'); tree.previewCategory('raizes');
  assert.deepEqual(cards.filter(c=>!c.classList.contains('dimmed')).map(c=>c.getAttribute('data-slug')),['bissecao']);
  tree.previewCategory('edos');
  assert.equal(cards.filter(c=>!c.classList.contains('dimmed')).length,0);
  tree.clearCategory();
  assert.equal(cards.filter(c=>!c.classList.contains('dimmed')).length,1);
});

test('indicador de animação acompanha apenas arquivos com resolvedor registrado', () => {
  const tree = Object.create(DecisionTreeComponent.prototype);
  tree.options = {hasAnimation:slug=>Boolean(getMethodKeyByFile(ALGORITHMS_DATABASE.find(a=>a.slug===slug).file))};
  tree.selectedSlug = null;
  for (const method of ALGORITHMS_DATABASE) {
    const svg=tree.buildCardSVG(method);
    assert.equal(svg.includes('class="card-animation-icon"'),Boolean(getMethodKeyByFile(method.file)),method.slug);
  }
  const svg=tree.buildSVGString();
  assert.equal((svg.match(/data-animation="true"/g)||[]).length,13);
  assert.equal((svg.match(/class="dt-category-column"/g)||[]).length,6);
});
