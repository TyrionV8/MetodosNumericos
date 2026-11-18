// Traduz os dados do motor para os mesmos painéis de gráfico, valores e tabela.
(function () {
  const format = value => Array.isArray(value) ? '[' + value.map(format).join(', ') + ']'
    : typeof value === 'number' && Number.isFinite(value)
      ? (Math.abs(value) >= 10000 || (value !== 0 && Math.abs(value) < 0.0001)
        ? value.toExponential(4) : Number(value.toFixed(6)).toString())
      : typeof value === 'string' ? value : '—';

  window.simulationPresentation = function (state) {
    const step = state.stepData;
    const family = state.manifest.family;
    let fields, headers, row;
    if (family === 'linear_system') {
      fields = [['LINHA DO PIVÔ', step.pivotRow >= 0 ? step.pivotRow+1 : null],
        ['COLUNA DO PIVÔ', step.pivotCol >= 0 ? step.pivotCol+1 : null], ['MULTIPLICADOR', step.multiplier]];
      headers = ['Passo', 'Operação', 'Linha do pivô', 'Coluna do pivô', 'Multiplicador'];
      row = s => [s.description, s.pivotRow >= 0 ? s.pivotRow+1 : null, s.pivotCol >= 0 ? s.pivotCol+1 : null, s.multiplier];
    } else if (family === 'quadrature') {
      fields = [['INÍCIO DO TRECHO', step.x0], ['FIM DO TRECHO', step.x2 ?? step.x1], ['ÁREA ACUMULADA', step.cumulativeArea]];
      headers = ['Passo', 'Início', 'Fim', 'Área do trecho', 'Área acumulada'];
      row = s => [s.x0, s.x2 ?? s.x1, s.sliceArea, s.cumulativeArea];
    } else if (family === 'ode_system') {
      fields = [['VALOR DE x', step.x], ['VALOR DE y', step.y], ['VALOR DE z', step.z]];
      headers = ['Passo', 'x', 'y', 'z', 'Próximo y', 'Próximo z'];
      row = s => [s.x, s.y, s.z, s.yNext, s.zNext];
    } else if (family === 'ode') {
      fields = [['VALOR DE x', step.x], ['VALOR DE y', step.y], ['ERRO ABSOLUTO DO PRÓXIMO PASSO', step.errAbs]];
      headers = ['Passo', 'x', 'y', 'Próximo y', 'Erro absoluto'];
      row = s => [s.x, s.y, s.yNext, s.errAbs];
    } else {
      const interval = state.methodKey === 'bissecao' || state.methodKey === 'falsa_posicao';
      fields = [['ESTIMATIVA ATUAL (x)', step.x], ['RESÍDUO |f(x)|', Math.abs(step.fx)],
        [interval ? 'LARGURA DO INTERVALO' : 'VARIAÇÃO |Δx|', Math.abs(step.delta)]];
      headers = ['Passo', 'x', 'f(x)', 'Derivada', 'Próximo x', interval ? 'Largura do intervalo' : '|Δx|'];
      row = s => [s.x, s.fx, s.dfx, s.xNext, typeof s.delta === 'number' ? Math.abs(s.delta) : null];
    }
    const status = step.description || (family === 'scalar_root'
      ? ({CONVERGED:'CONVERGIU',STEPPING:'EM ANDAMENTO',SINGULAR_DERIVATIVE:'DERIVADA PRÓXIMA DE ZERO',DIV_BY_ZERO:'DIVISÃO POR ZERO'}[step.status] || step.status)
      : state.currentStep === state.totalSteps-1 ? 'ÚLTIMO PASSO DO EXEMPLO' : 'EM ANDAMENTO');
    return {format, fields, headers, rows:state.methodData.steps.map((s,i)=>[i+1,...row(s)]), status};
  };
})();
