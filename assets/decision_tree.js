/**
 * decision_tree.js — Componente de Árvore de Decisão Interativa & Navegação Arcana
 * Projeto: Métodos Numéricos & Computação Científica (GNU Octave / MATLAB)
 * Estética: Arcane-Tech Cósmico /make-it-glow (#060b14, trilhas de silício, neon, JetBrains Mono)
 * Módulos: 10 módulos oficiais refatorados, 25 algoritmos canônicos
 */

(function (global) {
  'use strict';

  // Banco de Dados Canônico de Algoritmos (25 rotinas oficiais)
  const ALGORITHMS_DATABASE = [
    // ------------------------------------------------------------------------
    // Coluna 0: Raízes de Equações Não Lineares f(x) = 0
    // ------------------------------------------------------------------------
    {
      id: '11-metodo-da-bissecao',
      slug: '11-metodo-da-bissecao',
      catalogId: 'bissecao',
      engineKey: 'bissecao',
      shortTitle: 'Bisseção',
      title: 'Método da Bisseção Intervalar',
      category: 'raizes',
      colIndex: 0,
      rowIndex: 0,
      file: 'raizes_equacoes/bissecao.m',
      criterion: '• Troca de sinal f(a)·f(b) < 0',
      badge: 'O(log₂ 1/ε)',
      badgeType: 'order',
      subnote: 'Bolzano • f(x) contínua',
      convergence: 'Linear: O(log₂((b-a)/ε))',
      description: 'Garantia incondicional de convergência para funções contínuas pelo Teorema do Valor Intermediário de Bolzano através da bisseção sucessiva de intervalos.',
      formula: 'x_{k+1} = \\frac{a_k + b_k}{2}, \\quad k_{max} \\ge \\left\\lceil \\frac{\\ln(b-a) - \\ln(\\varepsilon)}{\\ln 2} \\right\\rceil',
      inputs: 'f(x) contínua, intervalo [a, b] com f(a)·f(b) < 0, tolerância tol, max_iter',
      outputs: 'Raiz estimada x*, resíduo |f(x*)|, número de iterações, tabela passo a passo',
      terminalCmd: 'octave --no-gui "raizes_equacoes/bissecao.m"'
    },
    {
      id: '12-metodo-da-falsa-posicao-regula-falsi',
      slug: '12-metodo-da-falsa-posicao-regula-falsi',
      catalogId: 'falsa-posicao',
      engineKey: 'falsa_posicao',
      shortTitle: 'Falsa Posição',
      title: 'Método da Falsa Posição (Regula Falsi)',
      category: 'raizes',
      colIndex: 0,
      rowIndex: 1,
      file: 'raizes_equacoes/falsa_posicao.m',
      criterion: '• Intervalo com reta secante',
      badge: 'Linear Acel.',
      badgeType: 'order',
      subnote: 'Média ponderada f(a), f(b)',
      convergence: 'Linear assintótica acelerada',
      description: 'Interpolação linear que preserva o enquadramento estrito da raiz, conectando (a, f(a)) e (b, f(b)) e calculando a interseção com o eixo das abscissas.',
      formula: 'x = \\frac{a \\cdot f(b) - b \\cdot f(a)}{f(b) - f(a)}',
      inputs: 'f(x) contínua, intervalo [a, b] com f(a)·f(b) < 0, tolerância tol, max_iter',
      outputs: 'Raiz x*, erro relativo percentual, histórico de convergência',
      terminalCmd: 'octave --no-gui "raizes_equacoes/falsa_posicao.m"'
    },
    {
      id: '13-metodo-de-newton-raphson-metodo-das-tangentes',
      slug: '13-metodo-de-newton-raphson-metodo-das-tangentes',
      catalogId: 'newton-raphson',
      engineKey: 'newton',
      shortTitle: 'Newton-Raphson',
      title: 'Método de Newton-Raphson (Tangentes)',
      category: 'raizes',
      colIndex: 0,
      rowIndex: 2,
      file: 'raizes_equacoes/newton_raphson.m',
      criterion: '• Com derivada analítica f\'(x)',
      badge: 'Ordem p = 2',
      badgeType: 'order',
      subnote: 'Duplica dígitos por iteração',
      convergence: 'Quadrática (Ordem 2 local)',
      description: 'Aproximação local pela reta tangente à curva f(x) no ponto x_k. Apresenta convergência assintótica quadrática de ordem p = 2 na vizinhança de raízes simples onde f\'(x) ≠ 0.',
      formula: 'x_{k+1} = x_k - \\frac{f(x_k)}{f\'(x_k)}',
      inputs: 'f(x), derivada analítica f\'(x), chute inicial x₀, tolerância tol, max_iter',
      outputs: 'Raiz refinada x*, resíduo ||f(x*)||, taxa quadrática demonstrada',
      terminalCmd: 'octave --no-gui "raizes_equacoes/newton_raphson.m"'
    },
    {
      id: '14-metodo-da-secante',
      slug: '14-metodo-da-secante',
      catalogId: 'secante',
      engineKey: 'secante',
      shortTitle: 'Secante',
      title: 'Método da Secante (Diferenças Finitas)',
      category: 'raizes',
      colIndex: 0,
      rowIndex: 3,
      file: 'raizes_equacoes/secante.m',
      criterion: '• Sem f\'(x) analítica (2 chutes)',
      badge: 'p ≈ 1.618',
      badgeType: 'order',
      subnote: 'Razão Áurea • Superlinear',
      convergence: 'Superlinear: p = (1 + √5)/2 ≈ 1.618',
      description: 'Substitui a derivada analítica de Newton por uma aproximação de diferença finita entre duas estimativas consecutivas x_{k-1} e x_k.',
      formula: 'x_{k+1} = x_k - f(x_k) \\frac{x_k - x_{k-1}}{f(x_k) - f(x_{k-1})}',
      inputs: 'f(x), duas aproximações iniciais x₀ e x₁, tolerância tol, max_iter',
      outputs: 'Raiz estimada x*, histórico iterativo, tabela de diferenças',
      terminalCmd: 'octave --no-gui "raizes_equacoes/secante.m"'
    },
    {
      id: '15-metodo-do-ponto-fixo-iteracao-linear-de-picard-banach',
      slug: '15-metodo-do-ponto-fixo-iteracao-linear-de-picard-banach',
      catalogId: 'ponto-fixo',
      engineKey: 'ponto_fixo',
      shortTitle: 'Ponto Fixo',
      title: 'Método do Ponto Fixo (Picard-Banach)',
      category: 'raizes',
      colIndex: 0,
      rowIndex: 4,
      file: 'raizes_equacoes/ponto_fixo.m',
      criterion: '• Forma x = g(x) com |g\'(x)| < 1',
      badge: 'Contração',
      badgeType: 'order',
      subnote: 'Picard • Teorema Contração',
      convergence: 'Linear: depende do fator de contração M = max |g\'(x)|',
      description: 'Reescreve a equação f(x) = 0 na forma equivalente x = g(x). Convergência garantida se a função g for uma contração (|g\'(x)| < 1) no domínio de busca.',
      formula: 'x_{k+1} = g(x_k), \\quad |x_{k+1} - x^*| \\le M |x_k - x^*| \\quad (M < 1)',
      inputs: 'Função de iteração g(x), chute x₀, tolerância tol, max_iter',
      outputs: 'Ponto fixo x*, verificação da condição de contração, tabela iterativa',
      terminalCmd: 'octave --no-gui "raizes_equacoes/ponto_fixo.m"'
    },

    // ------------------------------------------------------------------------
    // Coluna 1: Sistemas de Equações Lineares Ax = b
    // ------------------------------------------------------------------------
    {
      id: '21-eliminacao-gaussiana-com-pivoteamento-parcial-e-retrossubstituicao',
      slug: '21-eliminacao-gaussiana-com-pivoteamento-parcial-e-retrossubstituicao',
      catalogId: 'gauss-pivoteamento',
      engineKey: 'gauss_pivoteamento',
      shortTitle: 'Gauss Pivoteamento',
      title: 'Eliminação Gaussiana com Pivoteamento Parcial',
      category: 'sistemas',
      colIndex: 1,
      rowIndex: 0,
      file: 'sistemas_lineares/gauss_pivoteamento.m',
      criterion: '• Sistema determinado Ax = b',
      badge: 'O(n³/3)',
      badgeType: 'complexity',
      subnote: 'Triangular U + Retrossubst.',
      convergence: 'Método Direto Exato (salvo erros de arredondamento IEEE 754)',
      description: 'Triangularização da matriz de coeficientes com seleção de pivô de módulo máximo por coluna, evitando divisão por valores nulos e amplificação de resíduos.',
      formula: 'm_{ik} = \\frac{A_{ik}}{A_{kk}}, \\quad L_i \\leftarrow L_i - m_{ik} L_k, \\quad x_i = \\frac{b_i - \\sum_{j=i+1}^n U_{ij}x_j}{U_{ii}}',
      inputs: 'Matriz quadrada de coeficientes A (n×n), vetor coluna b (n×1)',
      outputs: 'Matriz triangular superior U, vetor solução x, resíduo ||Ax - b||₂',
      terminalCmd: 'octave --no-gui "sistemas_lineares/gauss_pivoteamento.m"'
    },
    {
      id: '22-metodo-de-gauss-jordan-com-normalizacao-unitaria',
      slug: '22-metodo-de-gauss-jordan-com-normalizacao-unitaria',
      catalogId: 'gauss-jordan',
      engineKey: 'gauss_jordan',
      shortTitle: 'Gauss-Jordan',
      title: 'Método de Gauss-Jordan com Normalização Unitária',
      category: 'sistemas',
      colIndex: 1,
      rowIndex: 1,
      file: 'sistemas_lineares/gauss_jordan.m',
      criterion: '• Diagonalização direta [I | x]',
      badge: 'O(n³/2)',
      badgeType: 'complexity',
      subnote: 'Sem retrossubstituição',
      convergence: 'Método Direto Exato',
      description: 'Transformação simultânea acima e abaixo da diagonal principal, convertendo a matriz aumentada [A|b] diretamente na forma escalonada reduzida [I|x].',
      formula: '[A \\mid b] \\xrightarrow{\\text{operações elementares}} [I_n \\mid x]',
      inputs: 'Matriz de coeficientes A (n×n), vetor de termos independentes b (n×1)',
      outputs: 'Matriz escalonada reduzida [I|x], vetor solução desacoplado x, resíduo',
      terminalCmd: 'octave --no-gui "sistemas_lineares/gauss_jordan.m"'
    },
    {
      id: '10-utilitarios-conversoes-e-visualizacao-grafica',
      slug: '10-utilitarios-conversoes-e-visualizacao-grafica',
      catalogId: 'operacoes-matrizes',
      engineKey: null,
      shortTitle: 'Álgebra & Matrizes',
      title: 'Diagnóstico Matricial e Condicionamento κ(A)',
      category: 'sistemas',
      colIndex: 1,
      rowIndex: 2,
      file: 'utilitarios/operacoes_matrizes.m',
      criterion: '• Diagnóstico, det(A), inv e κ(A)',
      badge: 'Diagnóstico',
      badgeType: 'metric',
      subnote: 'Condicionamento & Posto',
      convergence: 'Análise Espectral e Algébrica',
      description: 'Avaliação de integridade algébrica, cálculo de determinante, inversão analítica e número de condicionamento κ(A) para aferir sensibilidade a perturbações.',
      formula: '\\kappa(A) = \\|A\\| \\cdot \\|A^{-1}\\|, \\quad \\frac{\\|\\delta x\\|}{\\|x\\|} \\le \\kappa(A) \\frac{\\|\\delta b\\|}{\\|b\\|}',
      inputs: 'Matriz quadrada arbitrária A ou expressão matemática escalar',
      outputs: 'Determinante det(A), matriz inversa A⁻¹, número de condicionamento κ(A)',
      terminalCmd: 'octave --no-gui "utilitarios/operacoes_matrizes.m"'
    },

    // ------------------------------------------------------------------------
    // Coluna 2: Cálculo Diferencial e Integral Numérico
    // ------------------------------------------------------------------------
    {
      id: '31-regra-do-trapezio-composta',
      slug: '31-regra-do-trapezio-composta',
      catalogId: 'trapezio-composto',
      engineKey: 'trapezio_composto',
      shortTitle: 'Trapézio Composto',
      title: 'Regra do Trapézio Composta',
      category: 'calculo',
      colIndex: 2,
      rowIndex: 0,
      file: 'integracao_numerica/trapezio_composto.m',
      criterion: '• Integral definida genérica',
      badge: 'O(h²)',
      badgeType: 'order',
      subnote: 'Newton-Cotes Grau 1',
      convergence: 'Erro de Truncamento O(h²)',
      description: 'Quadratura numérica que aproxima a função f(x) por segmentos lineares em n subintervalos de comprimento uniforme h = (b - a)/n.',
      formula: 'I \\approx \\frac{h}{2} \\left[ f(a) + 2 \\sum_{i=1}^{n-1} f(x_i) + f(b) \\right], \\quad E_T = -\\frac{(b-a)h^2}{12} f\'\'(\\xi)',
      inputs: 'Integrando f(x), limites de integração [a, b], número de intervalos n',
      outputs: 'Valor aproximado da integral I, passo de malha h, estimativa do erro',
      terminalCmd: 'octave --no-gui "integracao_numerica/trapezio_composto.m"'
    },
    {
      id: '32-regra-de-simpson-13-composta',
      slug: '32-regra-de-simpson-13-composta',
      catalogId: 'simpson-13-composto',
      engineKey: 'simpson_13_composto',
      shortTitle: 'Simpson 1/3 Composto',
      title: 'Regra de Simpson 1/3 Composta',
      category: 'calculo',
      colIndex: 2,
      rowIndex: 1,
      file: 'integracao_numerica/simpson_13_composto.m',
      criterion: '• Alta exatidão (n par obrigatório)',
      badge: 'O(h⁴)',
      badgeType: 'order',
      subnote: 'Parábolas • Newton-Cotes 2',
      convergence: 'Erro de Truncamento O(h⁴)',
      description: 'Aproximação de f(x) por arcos de parábola de 2º grau aos pares de subintervalos. Apresenta erro de truncamento de quarta ordem O(h⁴) para funções com derivadas contínuas até quarta ordem.',
      formula: 'I \\approx \\frac{h}{3} \\left[ f(a) + 4\\sum_{\\text{ímpares}} f(x_i) + 2\\sum_{\\text{pares}} f(x_i) + f(b) \\right], \\; E_S = -\\frac{(b-a)h^4}{180} f^{(4)}(\\xi)',
      inputs: 'Integrando f(x), limites [a, b], número de intervalos n (estritamente par)',
      outputs: 'Valor numérico refinado da quadratura, estimativa de resíduo de alta ordem',
      terminalCmd: 'octave --no-gui "integracao_numerica/simpson_13_composto.m"'
    },
    {
      id: '4-diferenciacao-numerica-e-diferencas-finitas',
      slug: '4-diferenciacao-numerica-e-diferencas-finitas',
      catalogId: 'diferenca-progressiva-3pontos',
      engineKey: null,
      shortTitle: 'Dif. Progressiva 3P',
      title: 'Diferença Finita Progressiva de 3 Pontos',
      category: 'calculo',
      colIndex: 2,
      rowIndex: 2,
      file: 'diferenciacao_numerica/diferenca_progressiva_3pontos.m',
      criterion: '• Derivada f\'(x₀) no início da malha',
      badge: 'O(h²)',
      badgeType: 'order',
      subnote: 'Sem nós anteriores à esquerda',
      convergence: 'Ordem de Truncamento O(h²)',
      description: 'Fórmula assimétrica progressiva derivada da expansão em Série de Taylor para calcular a derivada primeira na fronteira inicial sem perder precisão quadrática.',
      formula: 'f\'(x_0) \\approx \\frac{-3f(x_0) + 4f(x_0+h) - f(x_0+2h)}{2h} + O(h^2)',
      inputs: 'Vetor discreto y ou função contínua f(x), passo de discretização h',
      outputs: 'Vetor de derivadas numéricas f\'(x), gráfico comparativo de declividade',
      terminalCmd: 'octave --no-gui "diferenciacao_numerica/diferenca_progressiva_3pontos.m"'
    },

    // ------------------------------------------------------------------------
    // Coluna 3: Equações Diferenciais Ordinárias (PVI)
    // ------------------------------------------------------------------------
    {
      id: '51-euler-explicito-primeira-ordem',
      slug: '51-euler-explicito-primeira-ordem',
      catalogId: 'euler-explicito',
      engineKey: 'euler_explicito',
      shortTitle: 'Euler Explícito',
      title: 'Método de Euler Explícito (Primeira Ordem)',
      category: 'edos',
      colIndex: 3,
      rowIndex: 0,
      file: 'equacoes_diferenciais/euler_explicito.m',
      criterion: '• PVI escalar simples y\' = f(x, y)',
      badge: 'O(h)',
      badgeType: 'order',
      subnote: 'Unistep • Avanço Tangente',
      convergence: 'Erro Global O(h), Erro Local O(h²)',
      description: 'Método unistep explícito fundamental para solução de Problemas de Valor Inicial (PVI). Avança o estado ao longo da inclinação da reta tangente local.',
      formula: 'y_{n+1} = y_n + h \\cdot f(x_n, y_n)',
      inputs: 'EDO f(x, y), condição inicial y(x₀) = y₀, passo h, intervalo [x₀, x_fim]',
      outputs: 'Vetores temporais x e y_num, tabela comparativa com solução exata',
      terminalCmd: 'octave --no-gui "equacoes_diferenciais/euler_explicito.m"'
    },
    {
      id: '52-euler-implicito-com-solucao-nao-linear-via-newton-raphson',
      slug: '52-euler-implicito-com-solucao-nao-linear-via-newton-raphson',
      catalogId: 'euler-implicito-newton',
      engineKey: 'euler_implicito',
      shortTitle: 'Euler Implícito + Newton',
      title: 'Euler Implícito com Solução Não Linear via Newton',
      category: 'edos',
      colIndex: 3,
      rowIndex: 1,
      file: 'equacoes_diferenciais/euler_implicito_newton.m',
      criterion: '• EDO rígida (stiff) ou não linear',
      badge: 'A-Estável',
      badgeType: 'stability',
      subnote: 'Incondicionalmente Estável',
      convergence: 'Ordem O(h) com estabilidade assintótica incondicional (A-estabilidade)',
      description: 'Método implícito regressivo onde o estado futuro y_{n+1} é raiz da equação algébrica F(y_{n+1}) = y_{n+1} - y_n - h f(x_{n+1}, y_{n+1}) = 0 resolvida por Newton.',
      formula: 'y_{n+1} = y_n + h \\cdot f(x_{n+1}, y_{n+1}), \\quad y_{n+1}^{(k+1)} = y_{n+1}^{(k)} - \\frac{F(y_{n+1}^{(k)})}{F\'(y_{n+1}^{(k)})}',
      inputs: 'EDO f(x, y), jacobiana ∂f/∂y, x₀, y₀, passo h, tolerância do Newton',
      outputs: 'Trajetória sem explosão numérica para problemas rígidos, resíduos iterativos',
      terminalCmd: 'octave --no-gui "equacoes_diferenciais/euler_implicito_newton.m"'
    },
    {
      id: '53-runge-kutta-de-2a-ordem-metodo-de-heun-escalar',
      slug: '53-runge-kutta-de-2a-ordem-metodo-de-heun-escalar',
      catalogId: 'rk2-individual',
      engineKey: 'rk2_heun',
      shortTitle: 'Heun RK2 Escalar',
      title: 'Runge-Kutta 2ª Ordem (Método de Heun Escalar)',
      category: 'edos',
      colIndex: 3,
      rowIndex: 2,
      file: 'equacoes_diferenciais/rk2_individual.m',
      criterion: '• PVI escalar com precisão O(h²)',
      badge: 'O(h²)',
      badgeType: 'order',
      subnote: 'Preditor-Corretor Trapezoidal',
      convergence: 'Erro Global O(h²), Erro Local O(h³)',
      description: 'Algoritmo de Heun de 2 estágios: calcula uma predição provisória de Euler no início do passo e a corrige usando a média das derivadas nas duas extremidades.',
      formula: 'k_1 = f(x_n, y_n), \\; k_2 = f(x_n + h, y_n + h k_1), \\; y_{n+1} = y_n + \\frac{h}{2}(k_1 + k_2)',
      inputs: 'EDO f(x, y), condição inicial y(x₀) = y₀, passo h, número de passos',
      outputs: 'Trajetória numérica discretizada, tabela comparativa com solução teórica',
      terminalCmd: 'octave --no-gui "equacoes_diferenciais/rk2_individual.m"'
    },
    {
      id: '54-runge-kutta-2-heun-para-sistema-acoplado-de-2-edos',
      slug: '54-runge-kutta-2-heun-para-sistema-acoplado-de-2-edos',
      catalogId: 'rk2-heun-sistema',
      engineKey: null,
      shortTitle: 'Heun RK2 Sistemas',
      title: 'Runge-Kutta 2 (Heun) para Sistema Acoplado de 2 EDOs',
      category: 'edos',
      colIndex: 3,
      rowIndex: 3,
      file: 'equacoes_diferenciais/rk2_heun_sistema.m',
      criterion: '• Sistema acoplado y\' = f e z\' = g',
      badge: 'Sistemas O(h²)',
      badgeType: 'order',
      subnote: '2 EDOs Simultâneas',
      convergence: 'Segunda Ordem Global O(h²)',
      description: 'Integração simultânea de duas EDOs de 1ª ordem acopladas y\' = f(x, y, z) e z\' = g(x, y, z), típica de osciladores harmônicos e dinâmica de corpos.',
      formula: '\\begin{cases} y_{n+1} = y_n + \\frac{h}{2}(k_{1y} + k_{2y}) \\\\ z_{n+1} = z_n + \\frac{h}{2}(k_{1z} + k_{2z}) \\end{cases}',
      inputs: 'Funções f(x,y,z) e g(x,y,z), condições iniciais y₀ e z₀, passo h',
      outputs: 'Curvas de estado temporais y(x) e z(x), gráfico de plano de fase',
      terminalCmd: 'octave --no-gui "equacoes_diferenciais/rk2_heun_sistema.m"'
    },

    // ------------------------------------------------------------------------
    // Coluna 4: Séries, Sinais & Utilitários
    // ------------------------------------------------------------------------
    {
      id: '6-series-de-funcoes-e-aproximacoes-polinomiais',
      slug: '6-series-de-funcoes-e-aproximacoes-polinomiais',
      shortTitle: 'Taylor cos(x)',
      title: 'Série de Taylor / Maclaurin para Cosseno',
      category: 'series_sinais',
      colIndex: 4,
      rowIndex: 0,
      file: 'series_aproximacoes/taylor_cosseno.m',
      criterion: '• Aproximação polinomial de cos(x)',
      badge: 'Maclaurin',
      badgeType: 'metric',
      subnote: 'Termos pares alternados',
      convergence: 'Convergência Absoluta em todo ℝ: Erro R_n(x) ~ O(x^{2N+2}/(2N+2)!)',
      description: 'Expansão em série de potências pares de Maclaurin para avaliação precisa de cos(x) sem tabelas trigonométricas, demonstrando o efeito da ordem de truncamento.',
      formula: '\\cos(x) = \\sum_{n=0}^{N} \\frac{(-1)^n x^{2n}}{(2n)!} + R_N(x)',
      inputs: 'Ângulo x em radianos, número de termos da expansão N',
      outputs: 'Valor aproximado cos_aprox, valor exato de biblioteca, erro relativo',
      terminalCmd: 'octave --no-gui "series_aproximacoes/taylor_cosseno.m"'
    },
    {
      id: '7-processamento-de-sinais-e-metricas-estatisticas',
      slug: '7-processamento-de-sinais-e-metricas-estatisticas',
      shortTitle: 'Síntese & RMS',
      title: 'Gerador de Formas de Onda e Cálculo de RMS',
      category: 'series_sinais',
      colIndex: 4,
      rowIndex: 1,
      file: 'processamento_sinais/gerador_sinais_rms.m',
      criterion: '• Sinais periódicos e valor eficaz',
      badge: 'RMS / V_ef',
      badgeType: 'metric',
      subnote: 'Senoidal, Quadrada, Triangular',
      convergence: 'Amostragem Discreta O(N)',
      description: 'Síntese digital de sinais senoidais, de onda quadrada e triangular, com cálculo numérico exato do valor médio, valor eficaz RMS e fator de crista.',
      formula: 'V_{\\text{RMS}} = \\sqrt{\\frac{1}{N} \\sum_{i=1}^N v_i^2}, \\quad V_{\\text{médio}} = \\frac{1}{N} \\sum_{i=1}^N v_i',
      inputs: 'Frequência fundamental f, amplitude A, offset contínuo DC, amostragem Fs',
      outputs: 'Formas de onda sintetizadas, valor médio, valor RMS eficaz, gráficos',
      terminalCmd: 'octave --no-gui "processamento_sinais/gerador_sinais_rms.m"'
    },
    {
      id: '10-utilitarios-conversoes-e-visualizacao-grafica-si',
      slug: '10-utilitarios-conversoes-e-visualizacao-grafica',
      shortTitle: 'Conversor SI / Imperial',
      title: 'Conversão Precisa de Unidades Físicas SI ↔ Imperial',
      category: 'series_sinais',
      colIndex: 4,
      rowIndex: 2,
      file: 'utilitarios/conversao_si_imperial.m',
      criterion: '• Conversão dimensional de engenharia',
      badge: 'SI / Imperial',
      badgeType: 'util',
      subnote: 'NIST • Massa & Comprimento',
      convergence: 'Fatores de Conversão Exatos por Definição',
      description: 'Conversor bidirecional preciso entre unidades do Sistema Internacional (SI) e unidades imperiais/anglo-saxãs (quilogramas, libras, onças, metros, polegadas, pés, jardas).',
      formula: '1 \\text{ lb} = 0{,}45359237 \\text{ kg}, \\quad 1 \\text{ in} = 0{,}0254 \\text{ m}',
      inputs: 'Massa m em quilogramas [kg] e comprimento L em metros [m]',
      outputs: 'Estrutura com valores equivalentes em lb, oz, in, ft e yd',
      terminalCmd: 'octave --no-gui "utilitarios/conversao_si_imperial.m"'
    },
    {
      id: '10-utilitarios-conversoes-e-visualizacao-grafica-plot',
      slug: '10-utilitarios-conversoes-e-visualizacao-grafica',
      shortTitle: 'Visualizador Gráfico',
      title: 'Visualizador Gráfico de Catenárias e Polinômios',
      category: 'series_sinais',
      colIndex: 4,
      rowIndex: 3,
      file: 'utilitarios/plotar_funcao.m',
      criterion: '• Visualização 2D de curvas & raízes',
      badge: 'Plot 2D',
      badgeType: 'util',
      subnote: 'Catenária • Tangentes',
      convergence: 'Renderização Gráfica Vetorial',
      description: 'Script utilitário para traçado de equações transcendentais da catenária física, localização geométrica de raízes e curvas de aproximação polinomial.',
      formula: 'y(x) = a \\cosh\\left(\\frac{x}{a}\\right) = \\frac{a}{2}\\left(e^{x/a} + e^{-x/a}\\right)',
      inputs: 'Equação analítica ou vetor de coeficientes, limites de plotagem',
      outputs: 'Janela gráfica anotada com indicação de pontos críticos e tangentes',
      terminalCmd: 'octave --no-gui "utilitarios/plotar_funcao.m"'
    },

    // ------------------------------------------------------------------------
    // Coluna 5: Físico-Química & Farmacocinética Aplicada
    // ------------------------------------------------------------------------
    {
      id: '81-equacao-de-born-mayer-para-energia-reticular',
      slug: '81-equacao-de-born-mayer-para-energia-reticular',
      shortTitle: 'Born-Mayer Reticular',
      title: 'Equação de Born-Mayer para Energia Reticular Iônica',
      category: 'aplicacoes',
      colIndex: 5,
      rowIndex: 0,
      file: 'fisico_quimica/born_mayer_energia_reticular.m',
      criterion: '• Energia de rede de cristais iônicos',
      badge: 'kJ/mol',
      badgeType: 'metric',
      subnote: 'Madelung • Repulsão Pauli',
      convergence: 'Modelo Eletrostático Mecano-Quântico',
      description: 'Cálculo da energia de rede cristalina U combinando a atração eletrostática de longo alcance de Madelung e o potencial repulsivo exponencial de Pauli.',
      formula: 'U = -\\frac{N_A A z_+ z_- e^2}{4\\pi\\varepsilon_0 r_0} \\left( 1 - \\frac{\\rho}{r_0} \\right)',
      inputs: 'Cargas iônicas z+, z-, constante de Madelung A, raios iônicos r+, r-, parâmetro ρ',
      outputs: 'Energia reticular teórica U (kJ/mol), comparação com ciclo de Born-Haber',
      terminalCmd: 'octave --no-gui "fisico_quimica/born_mayer_energia_reticular.m"'
    },
    {
      id: '82-momento-magnetico-e-susceptibilidade-molar',
      slug: '82-momento-magnetico-e-susceptibilidade-molar',
      shortTitle: 'Curie / Magnetismo',
      title: 'Momento Magnético Spin-Only e Susceptibilidade Molar',
      category: 'aplicacoes',
      colIndex: 5,
      rowIndex: 1,
      file: 'fisico_quimica/susceptibilidade_magnetica.m',
      criterion: '• Momento magnético e Lei de Curie',
      badge: 'Spin Puro',
      badgeType: 'metric',
      subnote: 'Complexos de Transição 3d',
      convergence: 'Fórmula de Spin Puro de Curie: χ_M = C / T',
      description: 'Determinação do momento magnético teórico spin-only μ_so para complexos de metais de transição e cálculo da susceptibilidade magnética molar paramagnética χ_M(T).',
      formula: '\\mu_{so} = \\sqrt{n(n+2)}\\,\\mu_B, \\quad \\chi_M = \\frac{N_A \\mu_0 \\mu_{eff}^2}{3 k_B T} = \\frac{C}{T}',
      inputs: 'Número de elétrons desemparelhados n, temperatura absoluta T [K]',
      outputs: 'Momento de spin μ_so em magnetons de Bohr (μ_B), susceptibilidade χ_M',
      terminalCmd: 'octave --no-gui "fisico_quimica/susceptibilidade_magnetica.m"'
    },
    {
      id: '83-determinacao-do-numero-de-eletrons-desemparelhados',
      slug: '83-determinacao-do-numero-de-eletrons-desemparelhados',
      shortTitle: 'Elétrons Desemparelhados',
      title: 'Resolução Inversa: Determinação de Elétrons via χ_M',
      category: 'aplicacoes',
      colIndex: 5,
      rowIndex: 2,
      file: 'fisico_quimica/eletrons_desemparelhados.m',
      criterion: '• Problema inverso: deduzir n de χ_M',
      badge: 'Inversão χ_M',
      badgeType: 'order',
      subnote: 'Raiz física n² + 2n - C = 0',
      convergence: 'Resolução Polinomial Quadrática',
      description: 'Inversão experimental da lei magnetoquímica para determinar a quantidade física inteira de elétrons d desemparelhados a partir da susceptibilidade medida em laboratório.',
      formula: 'n^2 + 2n - \\frac{3 k_B T \\chi_M}{N_A \\mu_0 \\mu_B^2} = 0 \\implies n = -1 + \\sqrt{1 + \\frac{3 k_B T \\chi_M}{N_A \\mu_0 \\mu_B^2}}',
      inputs: 'Susceptibilidade molar experimental χ_M, temperatura T [K]',
      outputs: 'Valor contínuo e número inteiro arredondado de elétrons desemparelhados n',
      terminalCmd: 'octave --no-gui "fisico_quimica/eletrons_desemparelhados.m"'
    },
    {
      id: '91-doses-repetidas-em-bolus-instantaneo',
      slug: '91-doses-repetidas-em-bolus-instantaneo',
      shortTitle: 'Farma: Bolus IV',
      title: 'Cinética Monocompartimental em Doses de Bolus IV',
      category: 'aplicacoes',
      colIndex: 5,
      rowIndex: 3,
      file: 'farmacocinetica/concentracao_dose_instantanea.m',
      criterion: '• Injeções em bolus sem absorção',
      badge: 'Bolus IV',
      badgeType: 'metric',
      subnote: 'Decaimento Exponencial 1ª Ordem',
      convergence: 'Modelo Monocompartimental Aberto',
      description: 'Modelagem do perfil de concentração plasmática de fármacos administrados por injeções endovenosas repetidas em bolus, com eliminação de primeira ordem.',
      formula: 'C(t) = \\sum_{i=0}^{m-1} \\frac{D}{V_d} e^{-k_e (t - i\\tau)}, \\quad C_{\\text{ss,max}} = \\frac{D/V_d}{1 - e^{-k_e \\tau}}',
      inputs: 'Dose padrão D, intervalo de dose τ, meia-vida de eliminação t₁/₂, volume V_d',
      outputs: 'Curva temporal C(t), concentração máxima, mínima no vale e no platô estacionário',
      terminalCmd: 'octave --no-gui "farmacocinetica/concentracao_dose_instantanea.m"'
    },
    {
      id: '92-doses-multiplas-com-cinetica-de-absorcao',
      slug: '92-doses-multiplas-com-cinetica-de-absorcao',
      shortTitle: 'Farma: Com Absorção',
      title: 'Cinética Monocompartimental com Absorção Contínua',
      category: 'aplicacoes',
      colIndex: 5,
      rowIndex: 4,
      file: 'farmacocinetica/concentracao_com_absorcao.m',
      criterion: '• Administração oral com taxa ka e ke',
      badge: 'Via Oral',
      badgeType: 'metric',
      subnote: 'Curva Suave • Pico t_max',
      convergence: 'Modelo Biexponencial de Bateman',
      description: 'Simulação farmacocinética com fase de absorção gastrointestinal contínua de taxa k_a seguida por eliminação sistêmica k_e, prevendo o pico plasmático C_max.',
      formula: 'C(t) = \\frac{F \\cdot D \\cdot k_a}{V_d (k_a - k_e)} \\left( e^{-k_e t} - e^{-k_a t} \\right), \\quad t_{\\max} = \\frac{\\ln(k_a / k_e)}{k_a - k_e}',
      inputs: 'Dose D, constante ka, meia-vida de eliminação, intervalo posológico τ',
      outputs: 'Curva suave de concentração plasmática, tempo de pico t_max e platô',
      terminalCmd: 'octave --no-gui "farmacocinetica/concentracao_com_absorcao.m"'
    },
    {
      id: '93-estudo-parametrico-intervalo-meia-vida',
      slug: '93-estudo-parametrico-intervalo-meia-vida',
      shortTitle: 'Farma: Razão τ / t₁/₂',
      title: 'Estudo Paramétrico: Razão Intervalo Posológico / Meia-Vida',
      category: 'aplicacoes',
      colIndex: 5,
      rowIndex: 5,
      file: 'farmacocinetica/estudo_intervalo_meia_vida.m',
      criterion: '• Otimização posológica τ / t₁/₂',
      badge: 'Acúmulo R',
      badgeType: 'metric',
      subnote: 'Fator R = 1 / (1 - e^(-kτ))',
      convergence: 'Varredura Paramétrica de Estabilidade',
      description: 'Análise da dependência do fator de acúmulo no estado estacionário R em função da razão adimensional entre o intervalo posológico e a meia-vida do fármaco.',
      formula: 'R = \\frac{1}{1 - e^{-k_e \\tau}} = \\frac{1}{1 - 2^{-\\tau / t_{1/2}}}, \\quad \\frac{\\tau}{t_{1/2}} \\in [0{,}25;\\, 4{,}0]',
      inputs: 'Vetor de razões de intervalo posológico τ, meia-vida de eliminação t₁/₂, dose padrão',
      outputs: 'Curva paramétrica do fator de acúmulo R, flutuação terapêutica percentual',
      terminalCmd: 'octave --no-gui "farmacocinetica/estudo_intervalo_meia_vida.m"'
    }
  ];

  // Metadados das 6 Colunas Principais da Árvore
  const CATEGORIES = [
    {
      id: 'raizes',
      colIndex: 0,
      title: 'Raízes de Equações: f(x) = 0',
      subtitle: '5 Métodos • Escalares',
      countLabel: '5 Algoritmos',
      accentColor: '#00f0ff',
      glowFilter: 'glowCyan',
      busColor: '#00f0ff',
      gradientId: 'gradCyan'
    },
    {
      id: 'sistemas',
      colIndex: 1,
      title: 'Sistemas Lineares: Ax = b',
      subtitle: '2 Diretos + Matrizes',
      countLabel: '3 Algoritmos',
      accentColor: '#818cf8',
      glowFilter: 'glowIndigo',
      busColor: '#818cf8',
      gradientId: 'gradIndigo'
    },
    {
      id: 'calculo',
      colIndex: 2,
      title: 'Cálculo Diferencial e Integral',
      subtitle: 'Quadratura & Diferenças',
      countLabel: '3 Algoritmos',
      accentColor: '#10b981',
      glowFilter: 'glowEmerald',
      busColor: '#10b981',
      gradientId: 'gradEmerald'
    },
    {
      id: 'edos',
      colIndex: 3,
      title: 'Equações Diferenciais (PVI)',
      subtitle: 'Euler Explícito, Impl. & RK2',
      countLabel: '4 Algoritmos',
      accentColor: '#f59e0b',
      glowFilter: 'glowAmber',
      busColor: '#f59e0b',
      gradientId: 'gradAmber'
    },
    {
      id: 'series_sinais',
      colIndex: 4,
      title: 'Séries e Sinais',
      subtitle: 'Taylor, RMS & Utilitários',
      countLabel: '4 Algoritmos',
      accentColor: '#38bdf8',
      glowFilter: 'glowSky',
      busColor: '#38bdf8',
      gradientId: 'gradSky'
    },
    {
      id: 'aplicacoes',
      colIndex: 5,
      title: 'Físico-Química e Farmacocinética',
      subtitle: 'Rede, Curie & Farmacocinética',
      countLabel: '6 Algoritmos',
      accentColor: '#ec4899',
      glowFilter: 'glowPink',
      busColor: '#ec4899',
      gradientId: 'gradPink'
    }
  ];

  // Geometria e Layout do SVG
  const GEOMETRY = {
    viewBoxWidth: 1600,
    viewBoxHeight: 800,
    marginX: 35,
    colWidth: 230,
    colGap: 30, // 35 + 6*230 + 5*30 = 35 + 1380 + 150 = 1565 (+ 35 right margin = 1600)
    rootX: 580,
    rootY: 18,
    rootWidth: 440,
    rootHeight: 52,
    busY: 104,
    categoryY: 138,
    categoryHeight: 48,
    trunkBottomY: 204,
    cardsStartY: 208,
    cardHeight: 84,
    cardGapY: 10
  };

  /**
   * Classe Principal do Componente da Árvore de Decisão
   */
  class DecisionTreeComponent {
    constructor(containerElement, options = {}) {
      this.container = typeof containerElement === 'string' 
        ? document.querySelector(containerElement) 
        : containerElement;

      if (!this.container) {
        throw new Error('DecisionTreeComponent: Elemento de contêiner não encontrado.');
      }

      this.options = Object.assign({
        onSelectAlgorithm: null, // Callback fn(slug, algorithm)
        isCompact: false,
        theme: 'dark-cosmic',
        initialFilter: '',
        initialCategory: 'all'
      }, options);

      this.currentFilter = this.options.initialFilter;
      this.currentCategory = this.options.initialCategory;
      this.selectedSlug = null;
      this.isCompact = this.options.isCompact;

      this.init();
    }

    init() {
      this.render();
      this.attachEventListeners();
    }

    /**
     * Renderiza o container completo (HUD superior de controle + SVG central + Gaveta de Inspeção)
     */
    render() {
      this.container.innerHTML = `
        <div class="decision-tree-wrapper cosmic-theme ${this.isCompact ? 'mode-compact' : 'mode-expanded'}">
          <!-- Barra de Controle e Filtros Holográfica -->
          <div class="dt-control-bar">
            <div class="dt-brand">
              <span class="dt-rune">⚡</span>
              <div class="dt-brand-text">
                <span class="dt-title">Árvore de Decisão Interativa</span>
                <span class="dt-subtitle">Qual arquivo e algoritmo usar? Selecione para navegar instantaneamente</span>
              </div>
            </div>

            <!-- Campo de Busca em Tempo Real -->
            <div class="dt-search-box">
              <span class="dt-search-icon">🔍</span>
              <input type="text" class="dt-search-input" placeholder="Buscar algoritmo, critério, ordem (ex.: stiff, simpson, O(h⁴))..." value="${this.escapeHtml(this.currentFilter)}">
              <button class="dt-clear-search ${this.currentFilter ? '' : 'hidden'}" title="Limpar busca">&times;</button>
            </div>

            <!-- Seletor de Modo de Visualização -->
            <div class="dt-view-actions">
              <button class="dt-btn dt-btn-compact ${this.isCompact ? 'active' : ''}" title="Alternar Modo Compacto / Expandido">
                <span class="dt-btn-icon">${this.isCompact ? '⤢ Expandir' : '⊟ Recolher'}</span>
              </button>
            </div>
          </div>

          <!-- Pílulas de Categoria -->
          <div class="dt-category-pills">
            <button class="dt-pill ${this.currentCategory === 'all' ? 'active' : ''}" data-category="all">
              Todas as 6 Áreas (25 Algoritmos)
            </button>
            ${CATEGORIES.map(cat => `
              <button class="dt-pill ${this.currentCategory === cat.id ? 'active' : ''}" data-category="${cat.id}" style="--pill-color: ${cat.accentColor};">
                ${cat.title.split(' ')[0]} ${cat.title.split(' ').slice(1).join(' ')} (${cat.countLabel.split(' ')[0]})
              </button>
            `).join('')}
          </div>

          <!-- Área do Diagrama SVG -->
          <div class="dt-canvas-area">
            ${this.buildSVGString()}
          </div>

          <!-- Drawer de Inspeção do Algoritmo Selecionado (Holográfico) -->
          <div class="dt-inspector-drawer" id="dtInspectorDrawer">
            <div class="dt-inspector-empty">
              <span class="dt-inspector-icon">⬡</span>
              <p>Passe o mouse para inspecionar ou clique em qualquer nó para abrir os parâmetros detalhados e navegar na Wiki.</p>
            </div>
            <div class="dt-inspector-content hidden" id="dtInspectorContent"></div>
          </div>
        </div>
      `;
    }

    /**
     * Constrói a string SVG completa com trilhas de silício, gradientes neon e cartões interativos
     */
    buildSVGString() {
      const g = GEOMETRY;
      const svgWidth = g.viewBoxWidth;
      const svgHeight = g.viewBoxHeight;

      // Calcular centros das colunas
      const colCenters = CATEGORIES.map((_, i) => g.marginX + i * (g.colWidth + g.colGap) + g.colWidth / 2);

      return `
        <svg class="dt-svg-element" viewBox="0 0 ${svgWidth} ${svgHeight}" preserveAspectRatio="xMidYMid meet" xmlns="http://www.w3.org/2000/svg">
          <defs>
            <!-- Padrão de Fundo de Grade Tecnológica -->
            <pattern id="siliconGrid" width="40" height="40" patternUnits="userSpaceOnUse">
              <path d="M 40 0 L 0 0 0 40" fill="none" stroke="#0e1a2f" stroke-width="0.75" stroke-dasharray="2 4"/>
              <circle cx="40" cy="40" r="1" fill="#1e3a5f" opacity="0.6"/>
            </pattern>

            <!-- Filtros de Brilho Neon (Make It Glow) -->
            <filter id="glowCyan" x="-20%" y="-20%" width="140%" height="140%">
              <feDropShadow dx="0" dy="0" stdDeviation="4" flood-color="#00f0ff" flood-opacity="0.6"/>
            </filter>
            <filter id="glowIndigo" x="-20%" y="-20%" width="140%" height="140%">
              <feDropShadow dx="0" dy="0" stdDeviation="4" flood-color="#818cf8" flood-opacity="0.6"/>
            </filter>
            <filter id="glowEmerald" x="-20%" y="-20%" width="140%" height="140%">
              <feDropShadow dx="0" dy="0" stdDeviation="4" flood-color="#10b981" flood-opacity="0.6"/>
            </filter>
            <filter id="glowAmber" x="-20%" y="-20%" width="140%" height="140%">
              <feDropShadow dx="0" dy="0" stdDeviation="4" flood-color="#f59e0b" flood-opacity="0.6"/>
            </filter>
            <filter id="glowSky" x="-20%" y="-20%" width="140%" height="140%">
              <feDropShadow dx="0" dy="0" stdDeviation="4" flood-color="#38bdf8" flood-opacity="0.6"/>
            </filter>
            <filter id="glowPink" x="-20%" y="-20%" width="140%" height="140%">
              <feDropShadow dx="0" dy="0" stdDeviation="4" flood-color="#ec4899" flood-opacity="0.6"/>
            </filter>
            <filter id="glowRoot" x="-15%" y="-15%" width="130%" height="130%">
              <feDropShadow dx="0" dy="2" stdDeviation="6" flood-color="#38bdf8" flood-opacity="0.45"/>
            </filter>
            <filter id="cardShadow" x="-10%" y="-10%" width="120%" height="120%">
              <feDropShadow dx="0" dy="4" stdDeviation="5" flood-color="#000000" flood-opacity="0.7"/>
            </filter>

            <!-- Gradientes de Cores Neon -->
            <linearGradient id="gradRoot" x1="0%" y1="0%" x2="100%" y2="100%">
              <stop offset="0%" stop-color="#0f2b48"/>
              <stop offset="50%" stop-color="#163c66"/>
              <stop offset="100%" stop-color="#091b30"/>
            </linearGradient>
            <linearGradient id="gradCyan" x1="0%" y1="0%" x2="100%" y2="0%">
              <stop offset="0%" stop-color="#00f0ff"/>
              <stop offset="100%" stop-color="#0284c7"/>
            </linearGradient>
            <linearGradient id="gradIndigo" x1="0%" y1="0%" x2="100%" y2="0%">
              <stop offset="0%" stop-color="#818cf8"/>
              <stop offset="100%" stop-color="#4f46e5"/>
            </linearGradient>
            <linearGradient id="gradEmerald" x1="0%" y1="0%" x2="100%" y2="0%">
              <stop offset="0%" stop-color="#10b981"/>
              <stop offset="100%" stop-color="#059669"/>
            </linearGradient>
            <linearGradient id="gradAmber" x1="0%" y1="0%" x2="100%" y2="0%">
              <stop offset="0%" stop-color="#f59e0b"/>
              <stop offset="100%" stop-color="#d97706"/>
            </linearGradient>
            <linearGradient id="gradSky" x1="0%" y1="0%" x2="100%" y2="0%">
              <stop offset="0%" stop-color="#38bdf8"/>
              <stop offset="100%" stop-color="#2563eb"/>
            </linearGradient>
            <linearGradient id="gradPink" x1="0%" y1="0%" x2="100%" y2="0%">
              <stop offset="0%" stop-color="#ec4899"/>
              <stop offset="100%" stop-color="#be185d"/>
            </linearGradient>
          </defs>

          <!-- 1. Fundo Cósmico com Grade de Silício -->
          <rect width="100%" height="100%" fill="#060b14"/>
          <rect width="100%" height="100%" fill="url(#siliconGrid)"/>

          <!-- 2. Linhas de Silício e Barramento Elétrico Principal (Silicon Bus) -->
          <g class="silicon-tracks" opacity="0.85">
            <!-- Linha Tronco Vertical que sai do Nó Raiz -->
            <line x1="800" y1="${g.rootY + g.rootHeight}" x2="800" y2="${g.busY}" stroke="#38bdf8" stroke-width="2.5"/>
            <circle cx="800" cy="${g.busY}" r="4" fill="#38bdf8" filter="url(#glowSky)"/>

            <!-- Barramento Horizontal Superior Interligando as 6 Categorias -->
            <path d="M ${colCenters[0]} ${g.busY} L ${colCenters[5]} ${g.busY}" stroke="#1e3a5f" stroke-width="3"/>
            <path d="M ${colCenters[0]} ${g.busY} L ${colCenters[5]} ${g.busY}" stroke="#38bdf8" stroke-width="1.5" stroke-dasharray="12 6" opacity="0.7" class="silicon-bus-glow"/>

            <!-- Derivações e Vias para Cada Coluna -->
            ${colCenters.map((cx, i) => `
              <g class="bus-branch" data-category="${CATEGORIES[i].id}">
                <!-- Linha Vertical até o Cabeçalho da Categoria -->
                <line x1="${cx}" y1="${g.busY}" x2="${cx}" y2="${g.categoryY}" stroke="${CATEGORIES[i].accentColor}" stroke-width="2"/>
                <!-- Nó de Conexão (Via de Silício) -->
                <circle cx="${cx}" cy="${g.busY}" r="4.5" fill="#060b14" stroke="${CATEGORIES[i].accentColor}" stroke-width="2"/>
                <circle cx="${cx}" cy="${g.busY}" r="2" fill="${CATEGORIES[i].accentColor}"/>
                <!-- Linha Tronco que sai da Categoria e vai para os Cards -->
                <line x1="${cx}" y1="${g.categoryY + g.categoryHeight}" x2="${cx}" y2="${g.trunkBottomY}" stroke="${CATEGORIES[i].accentColor}" stroke-width="1.8" stroke-dasharray="4 3"/>
              </g>
            `).join('')}
          </g>

          <!-- 3. Nó Raiz Central: Núcleo de Decisão Taxonômica -->
          <g class="dt-root-node" filter="url(#glowRoot)">
            <!-- Borda Externa Tecnológica com Chanfros -->
            <rect x="${g.rootX}" y="${g.rootY}" width="${g.rootWidth}" height="${g.rootHeight}" rx="10" fill="url(#gradRoot)" stroke="#38bdf8" stroke-width="1.8"/>
            <!-- Runa e Decorações nos Cantos -->
            <path d="M ${g.rootX + 6} ${g.rootY + 16} L ${g.rootX + 6} ${g.rootY + 6} L ${g.rootX + 16} ${g.rootY + 6}" fill="none" stroke="#38bdf8" stroke-width="2"/>
            <path d="M ${g.rootX + g.rootWidth - 6} ${g.rootY + 16} L ${g.rootX + g.rootWidth - 6} ${g.rootY + 6} L ${g.rootX + g.rootWidth - 16} ${g.rootY + 6}" fill="none" stroke="#38bdf8" stroke-width="2"/>
            <path d="M ${g.rootX + 6} ${g.rootY + g.rootHeight - 16} L ${g.rootX + 6} ${g.rootY + g.rootHeight - 6} L ${g.rootX + 16} ${g.rootY + g.rootHeight - 6}" fill="none" stroke="#38bdf8" stroke-width="2"/>
            <path d="M ${g.rootX + g.rootWidth - 6} ${g.rootY + g.rootHeight - 16} L ${g.rootX + g.rootWidth - 6} ${g.rootY + g.rootHeight - 6} L ${g.rootX + g.rootWidth - 16} ${g.rootY + g.rootHeight - 6}" fill="none" stroke="#38bdf8" stroke-width="2"/>
            <!-- Título do Nó Raiz -->
            <text x="800" y="${g.rootY + 23}" fill="#ffffff" font-family="'JetBrains Mono', monospace, system-ui" font-weight="700" font-size="13" text-anchor="middle" letter-spacing="0.5">
              ⚡ NÚCLEO DE DECISÃO COMPUTACIONAL
            </text>
            <text x="800" y="${g.rootY + 41}" fill="#94a3b8" font-family="system-ui, -apple-system, sans-serif" font-weight="500" font-size="11" text-anchor="middle">
              Qual é a classe do seu problema? Selecione a rota ótima
            </text>
          </g>

          <!-- 4. Cabeçalhos de Cada Categoria (Nós Mestre de Coluna) -->
          <g class="category-headers">
            ${CATEGORIES.map((cat, i) => {
              const x = g.marginX + i * (g.colWidth + g.colGap);
              const cx = x + g.colWidth / 2;
              return `
                <g class="cat-header-group" data-category="${cat.id}">
                  <rect x="${x}" y="${g.categoryY}" width="${g.colWidth}" height="${g.categoryHeight}" rx="8" fill="#0d182b" stroke="${cat.accentColor}" stroke-width="1.6" filter="url(#cardShadow)"/>
                  <!-- Faixa Superior Neon -->
                  <rect x="${x + 8}" y="${g.categoryY}" width="${g.colWidth - 16}" height="2.5" fill="url(#${cat.gradientId})"/>
                  <!-- Título da Categoria -->
                  <text x="${cx}" y="${g.categoryY + 22}" fill="${cat.accentColor}" font-family="'JetBrains Mono', monospace, system-ui" font-weight="700" font-size="12" text-anchor="middle">
                    ${cat.title}
                  </text>
                  <!-- Subtítulo e Contagem -->
                  <text x="${cx}" y="${g.categoryY + 38}" fill="#94a3b8" font-family="system-ui, sans-serif" font-size="10" text-anchor="middle">
                    ${cat.subtitle}
                  </text>
                </g>
              `;
            }).join('')}
          </g>

          <!-- 5. Cartões Interativos dos Algoritmos (Nós Folha Clicáveis) -->
          <g class="algorithm-leaf-nodes">
            ${ALGORITHMS_DATABASE.map(algo => this.buildCardSVG(algo)).join('')}
          </g>
        </svg>
      `;
    }

    /**
     * Gera o código SVG de um nó cartão individual interativo
     */
    buildCardSVG(algo) {
      const g = GEOMETRY;
      const cat = CATEGORIES[algo.colIndex];
      const cardX = g.marginX + algo.colIndex * (g.colWidth + g.colGap);
      const cardY = g.cardsStartY + algo.rowIndex * (g.cardHeight + g.cardGapY);
      const isSelected = this.selectedSlug === algo.slug;

      // Trilha de conexão que liga o nó anterior a este card
      const connectionLine = `
        <line x1="${cardX + g.colWidth / 2}" y1="${cardY - g.cardGapY}" x2="${cardX + g.colWidth / 2}" y2="${cardY}" stroke="${cat.accentColor}" stroke-width="1" stroke-dasharray="2 3" opacity="0.4"/>
      `;

      return `
        <g class="dt-algorithm-card ${isSelected ? 'selected' : ''}" 
           data-slug="${algo.slug}" 
           data-id="${algo.id}"
           data-category="${algo.category}"
           tabindex="0"
           role="button"
           aria-label="${algo.title}">
          ${connectionLine}

          <!-- Retângulo de Fundo do Card -->
          <rect class="card-bg" x="${cardX}" y="${cardY}" width="${g.colWidth}" height="${g.cardHeight}" rx="7" 
                fill="#0a1220" stroke="#1b2a47" stroke-width="1.3" filter="url(#cardShadow)"/>

          <!-- Indicador Lateral Neon da Área -->
          <rect class="card-edge" x="${cardX}" y="${cardY + 7}" width="3.5" height="${g.cardHeight - 14}" rx="1.5" fill="${cat.accentColor}"/>

          <!-- 1. Critério de Escolha / Pré-condição -->
          <text class="card-criterion" x="${cardX + 12}" y="${cardY + 17}" fill="#a5b4fc" font-family="system-ui, sans-serif" font-size="9.5" font-weight="600">
            ${this.truncateText(algo.criterion, 23)}
          </text>

          <!-- 2. Badge de Complexidade / Ordem de Convergência -->
          <g class="card-badge" transform="translate(${cardX + g.colWidth - 76}, ${cardY + 7})">
            <rect width="68" height="17" rx="3.5" fill="#0f1b2e" stroke="${cat.accentColor}" stroke-width="1"/>
            <text x="34" y="12" fill="${cat.accentColor}" font-family="'JetBrains Mono', monospace" font-size="8.5" font-weight="700" text-anchor="middle">
              ${algo.badge}
            </text>
          </g>

          <!-- 3. Nome do Algoritmo -->
          <text class="card-title" x="${cardX + 12}" y="${cardY + 36}" fill="#f8fafc" font-family="'JetBrains Mono', monospace" font-size="11.5" font-weight="700">
            ${algo.shortTitle}
          </text>

          <!-- 4. Caminho Monospaçado do Arquivo Oficial -->
          <text class="card-file" x="${cardX + 12}" y="${cardY + 53}" fill="#38bdf8" font-family="'JetBrains Mono', monospace" font-size="9">
            ${this.truncateText(algo.file, 28)}
          </text>

          <!-- 5. Nota Técnica / Propriedade Chave -->
          <text class="card-subnote" x="${cardX + 12}" y="${cardY + 70}" fill="#64748b" font-family="system-ui, sans-serif" font-size="8.8">
            ${algo.subnote}
          </text>

          <!-- Marcador de Ação / Ícone de Navegação -->
          <text class="card-nav-icon" x="${cardX + g.colWidth - 16}" y="${cardY + 71}" fill="${cat.accentColor}" font-family="'JetBrains Mono', monospace" font-size="10" font-weight="bold">
            ➔
          </text>
        </g>
      `;
    }

    /**
     * Vincula todos os eventos de clique, hover, busca e filtragem
     */
    attachEventListeners() {
      // 1. Clique nos cartões de algoritmo
      const cards = this.container.querySelectorAll('.dt-algorithm-card');
      cards.forEach(card => {
        const slug = card.getAttribute('data-slug');
        const algo = ALGORITHMS_DATABASE.find(a => a.slug === slug);

        card.addEventListener('click', (e) => {
          e.preventDefault();
          this.selectAlgorithm(slug);
        });

        card.addEventListener('mouseenter', () => {
          this.previewAlgorithm(algo);
        });

        card.addEventListener('keydown', (e) => {
          if (e.key === 'Enter' || e.key === ' ') {
            e.preventDefault();
            this.selectAlgorithm(slug);
          }
        });
      });

      // 2. Campo de busca em tempo real
      const searchInput = this.container.querySelector('.dt-search-input');
      const clearBtn = this.container.querySelector('.dt-clear-search');

      if (searchInput) {
        searchInput.addEventListener('input', (e) => {
          this.setFilter(e.target.value);
          if (clearBtn) {
            clearBtn.classList.toggle('hidden', !e.target.value);
          }
        });
      }

      if (clearBtn) {
        clearBtn.addEventListener('click', () => {
          searchInput.value = '';
          clearBtn.classList.add('hidden');
          this.setFilter('');
        });
      }

      // 3. Botões de filtro por categoria
      const pills = this.container.querySelectorAll('.dt-pill');
      pills.forEach(pill => {
        pill.addEventListener('click', () => {
          pills.forEach(p => p.classList.remove('active'));
          pill.classList.add('active');
          this.setCategory(pill.getAttribute('data-category'));
        });
      });

      // 4. Botão de alternância modo compacto / expandido
      const toggleBtn = this.container.querySelector('.dt-btn-compact');
      if (toggleBtn) {
        toggleBtn.addEventListener('click', () => {
          this.toggleCompactMode();
        });
      }
    }

    /**
     * Seleciona um algoritmo, atualiza o visual, dispara evento e abre detalhes
     */
    selectAlgorithm(slug) {
      const algo = ALGORITHMS_DATABASE.find(a => a.slug === slug);
      if (!algo) return;

      this.selectedSlug = slug;

      // Destacar o card no SVG
      const cards = this.container.querySelectorAll('.dt-algorithm-card');
      cards.forEach(c => {
        if (c.getAttribute('data-slug') === slug) {
          c.classList.add('selected');
          c.scrollIntoView({ behavior: 'smooth', block: 'nearest', inline: 'nearest' });
        } else {
          c.classList.remove('selected');
        }
      });

      // Renderizar o Inspector Drawer com o algoritmo selecionado
      this.displayInspectorDetails(algo);

      // Disparar o Callback de Opção
      if (typeof this.options.onSelectAlgorithm === 'function') {
        this.options.onSelectAlgorithm(slug, algo);
      }

      // Disparar Evento Customizado no Window para integração com a Wiki
      const event = new CustomEvent('select-algorithm', {
        detail: {
          slug: slug,
          algorithm: algo,
          timestamp: Date.now()
        },
        bubbles: true
      });
      window.dispatchEvent(event);
    }

    /**
     * Pré-visualiza o algoritmo no Inspector ao passar o mouse
     */
    previewAlgorithm(algo) {
      if (!algo) return;
      this.displayInspectorDetails(algo, true);
    }

    /**
     * Exibe os detalhes no painel/drawer holográfico
     */
    displayInspectorDetails(algo, isHover = false) {
      const contentEl = this.container.querySelector('#dtInspectorContent');
      const emptyEl = this.container.querySelector('.dt-inspector-empty');
      if (!contentEl) return;

      if (emptyEl) emptyEl.classList.add('hidden');
      contentEl.classList.remove('hidden');

      const cat = CATEGORIES[algo.colIndex];

      contentEl.innerHTML = `
        <div class="dt-inspector-card" style="--accent: ${cat.accentColor};">
          <div class="dt-insp-header">
            <span class="dt-insp-cat">${cat.title}</span>
            <span class="dt-insp-badge">${algo.badge} &bull; ${algo.convergence.split(':')[0]}</span>
          </div>

          <h3 class="dt-insp-title">${algo.title}</h3>
          
          <div class="dt-insp-file-bar">
            <span class="dt-insp-path"><code>${algo.file}</code></span>
            <button class="dt-btn-copy" data-copy="${algo.file}" title="Copiar caminho">📋 Copiar Caminho</button>
          </div>

          <p class="dt-insp-desc">${algo.description}</p>

          <div class="dt-insp-formula-box">
            <div class="dt-formula-label">📐 Modelo Matemático:</div>
            <div class="dt-formula-math">$$${algo.formula}$$</div>
          </div>

          <div class="dt-insp-grid">
            <div class="dt-insp-col">
              <span class="dt-label">📥 Entradas Necessárias:</span>
              <p class="dt-val">${algo.inputs}</p>
            </div>
            <div class="dt-insp-col">
              <span class="dt-label">📤 Saídas Entregues:</span>
              <p class="dt-val">${algo.outputs}</p>
            </div>
          </div>

          <div class="dt-insp-cmd-bar">
            <div class="dt-cmd-title">⚡ Executar no Terminal (Octave CLI):</div>
            <div class="dt-cmd-content">
              <code>${algo.terminalCmd}</code>
              <button class="dt-btn-copy" data-copy="${algo.terminalCmd}" title="Copiar comando">Copiar</button>
            </div>
          </div>

          <div class="dt-insp-actions">
            <a href="#${algo.slug}" class="dt-btn-navigate" data-slug="${algo.slug}">
              🚀 Abrir Seção Detalhada na Wiki &bull; [${algo.slug}]
            </a>
          </div>
        </div>
      `;

      // Renderizar fórmulas KaTeX se disponível na página global
      if (typeof renderMathInElement === 'function') {
        renderMathInElement(contentEl, {
          delimiters: [
            { left: '$$', right: '$$', display: true },
            { left: '$', right: '$', display: false }
          ]
        });
      }

      // Adicionar listeners para botões de cópia
      const copyBtns = contentEl.querySelectorAll('.dt-btn-copy');
      copyBtns.forEach(btn => {
        btn.addEventListener('click', (e) => {
          e.preventDefault();
          const text = btn.getAttribute('data-copy');
          navigator.clipboard.writeText(text).then(() => {
            const originalText = btn.textContent;
            btn.textContent = '✓ Copiado!';
            btn.classList.add('copied');
            setTimeout(() => {
              btn.textContent = originalText;
              btn.classList.remove('copied');
            }, 1800);
          });
        });
      });

      // Listener para o botão de navegação
      const navBtn = contentEl.querySelector('.dt-btn-navigate');
      if (navBtn) {
        navBtn.addEventListener('click', (e) => {
          const targetSlug = navBtn.getAttribute('data-slug');
          const targetElement = document.getElementById(targetSlug);
          if (targetElement) {
            e.preventDefault();
            targetElement.scrollIntoView({ behavior: 'smooth', block: 'start' });
            targetElement.classList.add('highlight-section');
            setTimeout(() => targetElement.classList.remove('highlight-section'), 2500);
          }
        });
      }
    }

    /**
     * Aplica filtro de busca textual nos nós da árvore
     */
    setFilter(queryText) {
      this.currentFilter = (queryText || '').toLowerCase().trim();
      this.applyVisibilityFilters();
    }

    /**
     * Filtra por categoria
     */
    setCategory(categoryId) {
      this.currentCategory = categoryId;
      this.applyVisibilityFilters();
    }

    /**
     * Aplica os estados visuais (ativo/esmaecido) nos cartões de acordo com filtros ativos
     */
    applyVisibilityFilters() {
      const cards = this.container.querySelectorAll('.dt-algorithm-card');
      const headers = this.container.querySelectorAll('.cat-header-group');
      const branches = this.container.querySelectorAll('.bus-branch');

      cards.forEach(card => {
        const slug = card.getAttribute('data-slug');
        const cat = card.getAttribute('data-category');
        const algo = ALGORITHMS_DATABASE.find(a => a.slug === slug);

        let matchesCategory = (this.currentCategory === 'all' || cat === this.currentCategory);
        let matchesQuery = true;

        if (this.currentFilter && algo) {
          const fullText = [
            algo.title,
            algo.shortTitle,
            algo.criterion,
            algo.badge,
            algo.file,
            algo.description,
            algo.subnote,
            algo.convergence
          ].join(' ').toLowerCase();

          matchesQuery = fullText.includes(this.currentFilter);
        }

        const isVisible = matchesCategory && matchesQuery;
        card.classList.toggle('dimmed', !isVisible);
        card.classList.toggle('match-highlight', Boolean(isVisible && this.currentFilter));
      });

      // Ajustar opacidade dos cabeçalhos das colunas
      headers.forEach(h => {
        const cat = h.getAttribute('data-category');
        const isCatVisible = (this.currentCategory === 'all' || cat === this.currentCategory);
        h.classList.toggle('dimmed', !isCatVisible);
      });

      branches.forEach(b => {
        const cat = b.getAttribute('data-category');
        const isCatVisible = (this.currentCategory === 'all' || cat === this.currentCategory);
        b.classList.toggle('dimmed', !isCatVisible);
      });
    }

    /**
     * Alterna entre modo compacto (cabeçalho resumido) e modo expandido
     */
    toggleCompactMode(forceState) {
      this.isCompact = forceState !== undefined ? forceState : !this.isCompact;
      const wrapper = this.container.querySelector('.decision-tree-wrapper');
      const btn = this.container.querySelector('.dt-btn-compact');

      if (wrapper) {
        wrapper.classList.toggle('mode-compact', this.isCompact);
        wrapper.classList.toggle('mode-expanded', !this.isCompact);
      }

      if (btn) {
        btn.classList.toggle('active', this.isCompact);
        btn.querySelector('.dt-btn-icon').textContent = this.isCompact ? '⤢ Expandir' : '⊟ Recolher';
      }
    }

    /**
     * Utilitário para truncar textos longos nos nós SVG
     */
    truncateText(text, maxChars) {
      if (!text) return '';
      return text.length > maxChars ? text.substring(0, maxChars - 1) + '…' : text;
    }

    escapeHtml(str) {
      if (!str) return '';
      return str.replace(/[&<>"']/g, m => ({
        '&': '&amp;',
        '<': '&lt;',
        '>': '&gt;',
        '"': '&quot;',
        "'": '&#39;'
      })[m]);
    }
  }

  // Injeção de Estilos CSS Arcane-Tech Cósmico
  const CSS_STYLES = `
    .decision-tree-wrapper.cosmic-theme {
      background-color: #060b14;
      color: #e2e8f0;
      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif;
      border: 1px solid #1a2a44;
      border-radius: 12px;
      overflow: hidden;
      box-shadow: 0 10px 30px rgba(0, 0, 0, 0.7);
      margin: 20px 0;
      position: relative;
    }

    /* Barra Superior de Controle */
    .dt-control-bar {
      display: flex;
      flex-wrap: wrap;
      align-items: center;
      justify-content: space-between;
      gap: 16px;
      padding: 16px 20px;
      background: #091222;
      border-bottom: 1px solid #16243d;
    }

    .dt-brand {
      display: flex;
      align-items: center;
      gap: 12px;
    }

    .dt-rune {
      font-size: 24px;
      color: #00f0ff;
      filter: drop-shadow(0 0 6px rgba(0, 240, 255, 0.7));
      animation: dtPulse 2.5s infinite ease-in-out;
    }

    @keyframes dtPulse {
      0%, 100% { transform: scale(1); opacity: 0.9; }
      50% { transform: scale(1.15); opacity: 1; filter: drop-shadow(0 0 12px rgba(0, 240, 255, 1)); }
    }

    .dt-brand-text {
      display: flex;
      flex-direction: column;
    }

    .dt-title {
      font-family: 'JetBrains Mono', monospace;
      font-weight: 700;
      font-size: 15px;
      color: #f8fafc;
      letter-spacing: 0.5px;
    }

    .dt-subtitle {
      font-size: 11.5px;
      color: #94a3b8;
    }

    /* Caixa de Busca */
    .dt-search-box {
      display: flex;
      align-items: center;
      background: #0d1a30;
      border: 1px solid #20355a;
      border-radius: 8px;
      padding: 6px 12px;
      min-width: 340px;
      flex-grow: 1;
      max-width: 520px;
      transition: all 0.25s ease;
    }

    .dt-search-box:focus-within {
      border-color: #00f0ff;
      box-shadow: 0 0 10px rgba(0, 240, 255, 0.25);
    }

    .dt-search-icon {
      font-size: 14px;
      margin-right: 8px;
      opacity: 0.7;
    }

    .dt-search-input {
      background: transparent;
      border: none;
      outline: none;
      color: #f1f5f9;
      font-family: 'JetBrains Mono', monospace;
      font-size: 12px;
      width: 100%;
    }

    .dt-search-input::placeholder {
      color: #64748b;
    }

    .dt-clear-search {
      background: transparent;
      border: none;
      color: #94a3b8;
      font-size: 18px;
      cursor: pointer;
      padding: 0 4px;
      line-height: 1;
    }

    .dt-clear-search:hover {
      color: #f43f5e;
    }

    .dt-view-actions .dt-btn {
      background: #112038;
      border: 1px solid #233b66;
      color: #e2e8f0;
      padding: 8px 14px;
      border-radius: 6px;
      font-family: 'JetBrains Mono', monospace;
      font-size: 12px;
      font-weight: 600;
      cursor: pointer;
      transition: all 0.2s ease;
    }

    .dt-view-actions .dt-btn:hover {
      background: #193158;
      border-color: #38bdf8;
      color: #38bdf8;
    }

    /* Pílulas de Categoria */
    .dt-category-pills {
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
      padding: 10px 20px;
      background: #080e1b;
      border-bottom: 1px solid #132036;
    }

    .dt-pill {
      background: #0f1c31;
      border: 1px solid #1e3355;
      color: #94a3b8;
      border-radius: 20px;
      padding: 5px 12px;
      font-size: 11px;
      font-weight: 600;
      cursor: pointer;
      transition: all 0.2s ease;
    }

    .dt-pill:hover {
      background: #162846;
      color: #f8fafc;
      border-color: var(--pill-color, #38bdf8);
    }

    .dt-pill.active {
      background: #162e52;
      color: #ffffff;
      border-color: var(--pill-color, #00f0ff);
      box-shadow: 0 0 8px rgba(0, 240, 255, 0.25);
    }

    /* Área do SVG */
    .dt-canvas-area {
      position: relative;
      width: 100%;
      overflow-x: auto;
      padding: 10px 0;
      background: #060b14;
    }

    .dt-svg-element {
      display: block;
      width: 100%;
      height: auto;
      min-width: 1100px;
    }

    /* Efeitos de Hover e Seleção nos Nós SVG */
    .dt-algorithm-card {
      cursor: pointer;
      outline: none;
      transition: transform 0.2s cubic-bezier(0.2, 0, 0, 1);
    }

    .dt-algorithm-card:hover rect.card-bg {
      fill: #101c33;
      stroke: #38bdf8;
      stroke-width: 1.8;
      filter: drop-shadow(0 0 8px rgba(56, 189, 248, 0.4));
    }

    .dt-algorithm-card:hover .card-title {
      fill: #00f0ff;
    }

    .dt-algorithm-card:hover .card-nav-icon {
      transform: translateX(3px);
    }

    .dt-algorithm-card.selected rect.card-bg {
      fill: #122442;
      stroke: #00f0ff;
      stroke-width: 2.2;
      filter: drop-shadow(0 0 12px rgba(0, 240, 255, 0.6));
    }

    .dt-algorithm-card.match-highlight rect.card-bg {
      stroke: #34d399;
      stroke-width: 2;
      filter: drop-shadow(0 0 8px rgba(52, 211, 153, 0.5));
    }

    .dt-algorithm-card.dimmed,
    .cat-header-group.dimmed,
    .bus-branch.dimmed {
      opacity: 0.22;
      transition: opacity 0.25s ease;
    }

    /* Drawer de Inspeção */
    .dt-inspector-drawer {
      background: #08101e;
      border-top: 1px solid #162642;
      padding: 18px 24px;
      min-height: 120px;
    }

    .dt-inspector-empty {
      display: flex;
      align-items: center;
      gap: 12px;
      color: #64748b;
      font-size: 13px;
    }

    .dt-inspector-empty .dt-inspector-icon {
      font-size: 22px;
      color: #38bdf8;
      opacity: 0.6;
    }

    .dt-inspector-card {
      display: flex;
      flex-direction: column;
      gap: 12px;
      animation: fadeIn 0.25s ease-out;
    }

    @keyframes fadeIn {
      from { opacity: 0; transform: translateY(6px); }
      to { opacity: 1; transform: translateY(0); }
    }

    .dt-insp-header {
      display: flex;
      justify-content: space-between;
      align-items: center;
    }

    .dt-insp-cat {
      font-family: 'JetBrains Mono', monospace;
      font-size: 12px;
      font-weight: 700;
      color: var(--accent, #00f0ff);
      text-transform: uppercase;
      letter-spacing: 0.5px;
    }

    .dt-insp-badge {
      background: #11223b;
      border: 1px solid var(--accent, #00f0ff);
      color: var(--accent, #00f0ff);
      padding: 3px 10px;
      border-radius: 4px;
      font-family: 'JetBrains Mono', monospace;
      font-size: 11px;
      font-weight: 700;
    }

    .dt-insp-title {
      font-family: 'JetBrains Mono', monospace;
      font-size: 18px;
      font-weight: 700;
      color: #f8fafc;
      margin: 0;
    }

    .dt-insp-file-bar {
      display: flex;
      align-items: center;
      justify-content: space-between;
      background: #0d192e;
      border: 1px solid #1c3258;
      padding: 6px 12px;
      border-radius: 6px;
    }

    .dt-insp-path code {
      color: #38bdf8;
      font-family: 'JetBrains Mono', monospace;
      font-size: 12px;
    }

    .dt-btn-copy {
      background: #172d50;
      border: 1px solid #2b4b80;
      color: #cbd5e1;
      padding: 4px 10px;
      border-radius: 4px;
      font-size: 11px;
      cursor: pointer;
      font-family: 'JetBrains Mono', monospace;
      transition: all 0.2s ease;
    }

    .dt-btn-copy:hover {
      background: #203f70;
      color: #ffffff;
      border-color: #38bdf8;
    }

    .dt-btn-copy.copied {
      background: #065f46;
      border-color: #10b981;
      color: #a7f3d0;
    }

    .dt-insp-desc {
      font-size: 13.5px;
      line-height: 1.5;
      color: #cbd5e1;
      margin: 0;
    }

    .dt-insp-formula-box {
      background: #060e1c;
      border: 1px solid #182c4e;
      border-radius: 6px;
      padding: 10px 14px;
    }

    .dt-formula-label {
      font-size: 11px;
      font-weight: 600;
      color: #94a3b8;
      margin-bottom: 4px;
    }

    .dt-formula-math {
      font-size: 14px;
      color: #f1f5f9;
      overflow-x: auto;
    }

    .dt-insp-grid {
      display: grid;
      grid-template-columns: 1fr 1fr;
      gap: 16px;
    }

    .dt-insp-col {
      background: #0c182c;
      border: 1px solid #192e52;
      border-radius: 6px;
      padding: 10px 12px;
    }

    .dt-label {
      font-size: 11px;
      font-weight: 700;
      color: #94a3b8;
      display: block;
      margin-bottom: 4px;
    }

    .dt-val {
      font-size: 12.5px;
      color: #e2e8f0;
      margin: 0;
      line-height: 1.4;
    }

    .dt-insp-cmd-bar {
      background: #091426;
      border: 1px solid #1c3258;
      border-radius: 6px;
      padding: 8px 12px;
    }

    .dt-cmd-title {
      font-size: 11px;
      font-weight: 600;
      color: #94a3b8;
      margin-bottom: 6px;
    }

    .dt-cmd-content {
      display: flex;
      justify-content: space-between;
      align-items: center;
      gap: 10px;
    }

    .dt-cmd-content code {
      font-family: 'JetBrains Mono', monospace;
      font-size: 12px;
      color: #34d399;
      word-break: break-all;
    }

    .dt-insp-actions {
      display: flex;
      margin-top: 4px;
    }

    .dt-btn-navigate {
      display: inline-flex;
      align-items: center;
      justify-content: center;
      width: 100%;
      background: linear-gradient(90deg, #0284c7 0%, #0369a1 100%);
      color: #ffffff;
      text-decoration: none;
      padding: 10px 18px;
      border-radius: 6px;
      font-family: 'JetBrains Mono', monospace;
      font-size: 13px;
      font-weight: 700;
      letter-spacing: 0.5px;
      transition: all 0.2s ease;
      box-shadow: 0 4px 12px rgba(2, 132, 199, 0.35);
    }

    .dt-btn-navigate:hover {
      background: linear-gradient(90deg, #0ea5e9 0%, #0284c7 100%);
      box-shadow: 0 4px 18px rgba(14, 165, 233, 0.55);
      transform: translateY(-1px);
    }

    /* Modo Compacto */
    .decision-tree-wrapper.mode-compact .dt-canvas-area {
      max-height: 280px;
      overflow-y: hidden;
    }

    .decision-tree-wrapper.mode-compact .dt-inspector-drawer {
      padding: 12px 18px;
    }

    .hidden {
      display: none !important;
    }

    /* Destaque ao rolar até a seção */
    .highlight-section {
      animation: pulseHighlight 2.5s ease-out;
    }

    @keyframes pulseHighlight {
      0% { outline: 3px solid #00f0ff; background-color: rgba(0, 240, 255, 0.15); }
      100% { outline: 3px solid transparent; background-color: transparent; }
    }
  `;

  // Injetar estilos no <head> se estiver no navegador
  if (typeof document !== 'undefined') {
    const styleTag = document.createElement('style');
    styleTag.id = 'decision-tree-cosmic-styles';
    styleTag.textContent = CSS_STYLES;
    document.head.appendChild(styleTag);
  }

  // Exportar para escopo global ou módulos
  const exportsObj = {
    DecisionTreeComponent,
    ALGORITHMS_DATABASE,
    CATEGORIES,
    GEOMETRY
  };

  if (typeof module !== 'undefined' && module.exports) {
    module.exports = exportsObj;
  } else {
    global.DecisionTreeComponent = DecisionTreeComponent;
    global.DECISION_TREE_DATA = exportsObj;
  }

})(typeof window !== 'undefined' ? window : this);
