// Gera a versão clássica do catálogo para uso direto em file://.
// O original mantém suas exportações como módulo JavaScript.
const fs = require('node:fs');
const path = require('node:path');
const assets = path.join(__dirname, '..', 'assets');
const source = fs.readFileSync(path.join(assets, 'catalog.js'), 'utf8');
const browser = source.replace(/^export default ALGORITHMS;$/m, '').replace(/^export /gm, '');
fs.writeFileSync(path.join(assets, 'catalog.browser.js'), '// Gerado de catalog.js para funcionamento também em file://.\n' + browser);
