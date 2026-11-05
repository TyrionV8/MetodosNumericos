/**
 * =========================================================================
 * WIKI DE MATEMÁTICA COMPUTACIONAL — MOTOR MATEMÁTICO & CANVAS CRT
 * REPOSITÓRIO: /home/davi/repos/MetodosNumericos
 * ARQUIVO: engine.js
 * 
 * Descrição:
 * Motor matemático de cálculo numérico em tempo real e renderizador vetorial
 * de osciloscópio CRT (estilo terminal militar / workstation curses) para
 * algoritmos científicos de Métodos Numéricos em IEEE-754 FP64.
 * 
 * Famílias implementadas:
 *  1. Raízes de Equações Escalares:
 *     - Newton-Raphson (convergência quadrática O(h²))
 *     - Bisseção (convergência linear global via Teorema de Bolzano)
 *     - Secante (superlinear p ≈ 1.618)
 *     - Falsa Posição / Regula Falsi (interpolação linear de bracketing)
 *     - Ponto Fixo de Picard-Banach (iteração de contração x = g(x))
 *  2. Sistemas de Equações Lineares:
 *     - Eliminação Gaussiana com Pivoteamento Parcial e Retrossubstituição
 *     - Método de Gauss-Jordan com Normalização Unitária Diagonal
 *  3. Quadratura Numérica (Newton-Cotes):
 *     - Regra do Trapézio Composta (convergência O(h²))
 *     - Regra de Simpson 1/3 Composta (convergência O(h⁴))
 *  4. Equações Diferenciais Ordinárias (PVI):
 *     - Método de Euler Explícito (1ª ordem O(h))
 *     - Método de Euler Implícito com Solução Não Linear via Newton-Raphson (A-estável)
 *     - Runge-Kutta de 2ª Ordem / Método de Heun (preditor-corretor O(h²))
 *  5. Controles Interativos & Teclado:
 *     - play(), pause(), togglePlay(), stepForward(), stepBackward(), reset(), setSpeed(), loadMethod(), seek()
 *     - Atalhos de teclado: [SPACE], [, ], R, 1-9
 * =========================================================================
 */

(function (global) {
  'use strict';

  // =========================================================================
  // 1. TEMA VISUAL E PALETA CRT OSCILLOSCOPE
  // =========================================================================
  const CRTTheme = {
    bg: '#030508',
    bgPanel: '#060b12',
    bgCard: '#09121d',
    border: '#1a324d',
    borderBright: '#00ff66',
    
    // Fósforos CRT
    green: '#00ff66',
    greenGlow: 'rgba(0, 255, 102, 0.55)',
    greenDim: '#008f3b',
    
    cyan: '#00f0ff',
    cyanGlow: 'rgba(0, 240, 255, 0.55)',
    cyanDim: '#0a4b5c',
    
    amber: '#ffb700',
    amberGlow: 'rgba(255, 183, 0, 0.5)',
    amberDim: '#5c4300',
    
    red: '#ff3366',
    redGlow: 'rgba(255, 51, 102, 0.6)',
    
    text: '#e2f0fc',
    textDim: '#7da5c9',
    textBright: '#ffffff',
    
    reticleLine: '#0e2033',
    axisLine: '#294f77',
    
    font: '15px "JetBrains Mono", "Courier New", Courier, monospace',
    fontBold: 'bold 15px "JetBrains Mono", "Courier New", Courier, monospace',
    fontLarge: 'bold 20px "JetBrains Mono", "Courier New", Courier, monospace',
    fontSmall: '12px "JetBrains Mono", "Courier New", Courier, monospace'
  };

  // =========================================================================
  // 2. UTILITÁRIOS MATEMÁTICOS E DE FORMATAÇÃO
  // =========================================================================
  const MathUtils = {
    clamp(val, min, max) {
      return Math.max(min, Math.min(max, val));
    },

    formatFloat(num, decimals = 6) {
      if (num === undefined || num === null || Number.isNaN(num)) return 'NaN';
      if (!Number.isFinite(num)) return num > 0 ? '+Inf' : '-Inf';
      if (Math.abs(num) < 1e-12) return '0.000000';
      if (Math.abs(num) >= 1e6 || Math.abs(num) < 1e-4) {
        return num.toExponential(decimals);
      }
      return num.toFixed(decimals);
    },

    signPrefix(num, decimals = 6) {
      if (num > 0) return '+' + this.formatFloat(num, decimals);
      return this.formatFloat(num, decimals);
    },

    // Mapeamento linear para coordenadas de tela
    mapRange(x, inMin, inMax, outMin, outMax) {
      if (Math.abs(inMax - inMin) < 1e-12) return (outMin + outMax) / 2;
      return outMin + ((x - inMin) / (inMax - inMin)) * (outMax - outMin);
    }
  };

  // =========================================================================
  // 3. FAMÍLIA 1: RAÍZES DE EQUAÇÕES ESCALARES (f(x) = 0)
  // =========================================================================
  const ScalarRootSolvers = {
    /**
     * Newton-Raphson: x_{k+1} = x_k - f(x_k)/f'(x_k)
     * Convergência quadrática O(h²)
     */
    solveNewton(params = {}) {
      const f = params.f || ((x) => x * x - 2);
      const df = params.df || ((x) => 2 * x);
      const x0 = params.x0 !== undefined ? params.x0 : 3.0;
      const tol = params.tol || 1e-6;
      const maxIter = params.maxIter || 20;
      const target = params.target !== undefined ? params.target : Math.SQRT2;

      const steps = [];
      let xCurr = x0;

      for (let k = 0; k < maxIter; k++) {
        const fx = f(xCurr);
        const dfx = df(xCurr);

        if (Math.abs(dfx) < 1e-14) {
          steps.push({
            step: k,
            k,
            x: xCurr,
            fx,
            dfx,
            xNext: xCurr,
            delta: 0,
            err: Math.abs(fx),
            status: 'SINGULAR_DERIVATIVE',
            line: 'ERROR: f\'(x) ≈ 0 (derivada nula / singularidade)'
          });
          break;
        }

        const delta = -(fx / dfx);
        const xNext = xCurr + delta;
        const err = Math.abs(delta);
        const converged = Math.abs(fx) <= tol || err <= tol;

        steps.push({
          step: k,
          k,
          x: xCurr,
          fx,
          dfx,
          xNext,
          delta,
          err,
          status: converged ? 'CONVERGED' : 'STEPPING',
          line: `xi = x0 - (F(x0)/dF(x0)); % Δx = ${delta.toFixed(6)}`
        });

        if (converged) break;
        xCurr = xNext;
      }

      return {
        key: 'newton',
        title: 'NEWTON-RAPHSON',
        category: 'Raízes de Equações Escalares',
        file: 'raizes_equacoes/newton_raphson.m',
        order: 'O(h²) QUADRÁTICA',
        stability: 'LOCAL (f\'(x*) ≠ 0)',
        formula: 'x_{k+1} = x_k - [ f(x_k) / f\'(x_k) ]',
        formulaSub: 'Truncamento de Taylor O(Δx²)',
        funcLabel: 'f(x) = x² - 2 // f\'(x) = 2x',
        target,
        xRange: [0.0, 3.6],
        yRange: [-2.5, 8.5],
        f,
        df,
        steps
      };
    },

    /**
     * Bisseção: Ponto médio x_m = (a + b) / 2
     * Convergência linear garantida pelo Teorema de Bolzano
     */
    solveBissecao(params = {}) {
      const f = params.f || ((x) => x * x - 2);
      let a = params.a !== undefined ? params.a : 0.5;
      let b = params.b !== undefined ? params.b : 2.5;
      const tol = params.tol || 1e-4;
      const maxIter = params.maxIter || 25;
      const target = params.target !== undefined ? params.target : Math.SQRT2;

      const steps = [];

      for (let k = 0; k < maxIter; k++) {
        const fa = f(a);
        const fb = f(b);
        const xm = (a + b) / 2.0;
        const fxm = f(xm);
        const err = (b - a) / 2.0;
        const converged = err <= tol || Math.abs(fxm) <= tol;

        const nextA = fa * fxm < 0 ? a : xm;
        const nextB = fa * fxm < 0 ? xm : b;

        steps.push({
          step: k,
          k,
          a,
          b,
          xm,
          x: xm,
          fa,
          fb,
          fx: fxm,
          dfx: null,
          xNext: xm,
          nextA,
          nextB,
          err,
          delta: b - a,
          status: converged ? 'CONVERGED' : 'STEPPING',
          line: fa * fxm < 0 ? 'b = PontoMedio; % troca de sinal em [a, xm]' : 'a = PontoMedio; % troca de sinal em [xm, b]'
        });

        if (converged) break;
        a = nextA;
        b = nextB;
      }

      return {
        key: 'bissecao',
        title: 'MÉTODO DA BISSEÇÃO',
        category: 'Raízes de Equações Escalares',
        file: 'raizes_equacoes/bissecao.m',
        order: 'O(h) LINEAR (fator 1/2)',
        stability: 'GLOBAL (f(a)·f(b) < 0)',
        formula: 'x_m = (a + b) / 2 // Erro: (b - a)/2',
        formulaSub: 'Teorema de Bolzano / Dicotomia',
        funcLabel: 'f(x) = x² - 2 // Intervalo [a, b]',
        target,
        xRange: [0.0, 3.2],
        yRange: [-2.5, 6.0],
        f,
        steps
      };
    },

    /**
     * Secante: aproximação da derivada por diferenças finitas entre x_k e x_{k-1}
     * Ordem superlinear: p = (1 + √5)/2 ≈ 1.618
     */
    solveSecante(params = {}) {
      const f = params.f || ((x) => x * x - 2);
      let xPrev = params.x0 !== undefined ? params.x0 : 3.0;
      let xCurr = params.x1 !== undefined ? params.x1 : 2.5;
      const tol = params.tol || 1e-6;
      const maxIter = params.maxIter || 20;
      const target = params.target !== undefined ? params.target : Math.SQRT2;

      const steps = [];

      for (let k = 0; k < maxIter; k++) {
        const fPrev = f(xPrev);
        const fCurr = f(xCurr);
        const denom = fCurr - fPrev;

        if (Math.abs(denom) < 1e-14) {
          steps.push({
            step: k,
            k,
            xPrev,
            x: xCurr,
            fx: fCurr,
            dfx: 0,
            xNext: xCurr,
            delta: 0,
            err: Math.abs(fCurr),
            status: 'DIV_BY_ZERO',
            line: 'ERROR: f(x_k) - f(x_{k-1}) ≈ 0'
          });
          break;
        }

        const xNext = xCurr - fCurr * ((xCurr - xPrev) / denom);
        const delta = xNext - xCurr;
        const err = Math.abs(delta);
        const approxDfx = denom / (xCurr - xPrev);
        const converged = Math.abs(fCurr) <= tol || err <= tol;

        steps.push({
          step: k,
          k,
          xPrev,
          x: xCurr,
          fx: fCurr,
          dfx: approxDfx,
          fPrev,
          xNext,
          delta,
          err,
          status: converged ? 'CONVERGED' : 'STEPPING',
          line: 'x_prox = x - f*(x - x_ant)/(f - f_ant);'
        });

        if (converged) break;
        xPrev = xCurr;
        xCurr = xNext;
      }

      return {
        key: 'secante',
        title: 'MÉTODO DA SECANTE',
        category: 'Raízes de Equações Escalares',
        file: 'raizes_equacoes/secante.m',
        order: 'SUPERLINEAR (p ≈ 1.618)',
        stability: 'LOCAL (sem derivada analítica)',
        formula: 'x_{k+1} = x_k - f(x_k) \\cdot \\frac{x_k - x_{k-1}}{f(x_k) - f(x_{k-1})}',
        formulaSub: 'Diferenças finitas sucessivas',
        funcLabel: 'f(x) = x² - 2 // Chutes [3.0, 2.5]',
        target,
        xRange: [0.5, 3.5],
        yRange: [-2.5, 8.0],
        f,
        steps
      };
    },

    /**
     * Falsa Posição (Regula Falsi): interseção da reta secante com y = 0 com manutenção de bracketing
     */
    solveFalsaPosicao(params = {}) {
      const f = params.f || ((x) => x * x - 2);
      let a = params.a !== undefined ? params.a : 0.5;
      let b = params.b !== undefined ? params.b : 2.5;
      const tol = params.tol || 1e-4;
      const maxIter = params.maxIter || 25;
      const target = params.target !== undefined ? params.target : Math.SQRT2;

      const steps = [];
      let xAnt = a;

      for (let k = 0; k < maxIter; k++) {
        const fa = f(a);
        const fb = f(b);
        const denom = fb - fa;

        if (Math.abs(denom) < 1e-14) break;

        const xr = (a * fb - b * fa) / denom;
        const fxr = f(xr);
        const err = Math.abs(xr - xAnt);
        const converged = Math.abs(fxr) <= tol || (k > 0 && err <= tol);

        const nextA = fa * fxr < 0 ? a : xr;
        const nextB = fa * fxr < 0 ? xr : b;

        steps.push({
          step: k,
          k,
          a,
          b,
          fa,
          fb,
          xr,
          x: xr,
          fx: fxr,
          dfx: denom / (b - a),
          xNext: xr,
          nextA,
          nextB,
          err,
          delta: b - a,
          status: converged ? 'CONVERGED' : 'STEPPING',
          line: 'xr = (a*f(b) - b*f(a)) / (f(b) - f(a));'
        });

        if (converged) break;
        xAnt = xr;
        a = nextA;
        b = nextB;
      }

      return {
        key: 'falsa_posicao',
        title: 'FALSA POSIÇÃO (REGULA FALSI)',
        category: 'Raízes de Equações Escalares',
        file: 'raizes_equacoes/falsa_posicao.m',
        order: 'LINEAR / SUPERLINEAR',
        stability: 'GLOBAL (f(a)·f(b) < 0)',
        formula: 'x_r = \\frac{a \\cdot f(b) - b \\cdot f(a)}{f(b) - f(a)}',
        formulaSub: 'Interpolação linear sob bracketing',
        funcLabel: 'f(x) = x² - 2 // [a, b]',
        target,
        xRange: [0.0, 3.2],
        yRange: [-2.5, 6.0],
        f,
        steps
      };
    },

    /**
     * Ponto Fixo (MPF / Picard-Banach): x_{k+1} = g(x_k)
     * Teorema da Contração: |g'(x)| <= L < 1
     */
    solvePontoFixo(params = {}) {
      // Para f(x) = x² - 2 = 0, forma de contração babilônica: g(x) = 0.5 * (x + 2/x)
      const f = params.f || ((x) => x * x - 2);
      const g = params.g || ((x) => 0.5 * (x + 2.0 / x));
      let xCurr = params.x0 !== undefined ? params.x0 : 2.5;
      const tol = params.tol || 1e-6;
      const maxIter = params.maxIter || 20;
      const target = params.target !== undefined ? params.target : Math.SQRT2;

      const steps = [];

      for (let k = 0; k < maxIter; k++) {
        const gx = g(xCurr);
        const fx = f(xCurr);
        const delta = gx - xCurr;
        const err = Math.abs(delta);
        const converged = Math.abs(fx) <= tol || err <= tol;

        steps.push({
          step: k,
          k,
          x: xCurr,
          gx,
          fx,
          dfx: null,
          xNext: gx,
          delta,
          err,
          status: converged ? 'CONVERGED' : 'STEPPING',
          line: 'x_prox = g(x); % Iteração de Banach'
        });

        if (converged) break;
        xCurr = gx;
      }

      return {
        key: 'ponto_fixo',
        title: 'MÉTODO DO PONTO FIXO (PICARD-BANACH)',
        category: 'Raízes de Equações Escalares',
        file: 'raizes_equacoes/ponto_fixo.m',
        order: 'O(h) LINEAR (contrativo |g\'| < 1)',
        stability: 'LOCAL (depende da escolha de g(x))',
        formula: 'x_{k+1} = g(x_k) \\iff f(x) = 0',
        formulaSub: 'Teorema do Ponto Fixo de Banach',
        funcLabel: 'g(x) = 0.5·(x + 2/x) // Alvo √2',
        target,
        xRange: [0.5, 3.2],
        yRange: [0.0, 3.5],
        f,
        g,
        steps
      };
    }
  };

  // =========================================================================
  // 4. FAMÍLIA 2: SISTEMAS LINEARES (A x = b)
  // =========================================================================
  const LinearSystemSolvers = {
    /**
     * Eliminação de Gauss com Pivoteamento Parcial e Retrossubstituição
     */
    solveGaussPivoteamento(A_in, b_in) {
      // Sistema padrão validado do repositório
      const A = A_in ? JSON.parse(JSON.stringify(A_in)) : [
        [4, 3, 2],
        [1, 3, 1],
        [2, 1, 3]
      ];
      const b = b_in ? [...b_in] : [960, 510, 610];
      const n = A.length;

      // Monta matriz aumentada [A|b]
      let M = A.map((row, i) => [...row, b[i]]);
      const steps = [];

      steps.push({
        step: 0,
        type: 'initial',
        stepIndex: 0,
        matrix: JSON.parse(JSON.stringify(M)),
        pivotRow: -1,
        pivotCol: -1,
        activeRow: -1,
        activeCol: -1,
        multiplier: 0,
        description: 'Matriz Aumentada Inicial [A | b]',
        line: 'M = [A, b]; % Início do escalonamento'
      });

      // Fase de Eliminação Descendente
      for (let k = 0; k < n - 1; k++) {
        // 1. Pivoteamento Parcial: busca linha com maior |a_{ik}| para i = k..n-1
        let maxVal = Math.abs(M[k][k]);
        let pivotRow = k;
        for (let p = k + 1; p < n; p++) {
          if (Math.abs(M[p][k]) > maxVal) {
            maxVal = Math.abs(M[p][k]);
            pivotRow = p;
          }
        }

        if (pivotRow !== k) {
          // Troca de linhas
          const temp = M[k];
          M[k] = M[pivotRow];
          M[pivotRow] = temp;

          steps.push({
            step: steps.length,
            type: 'pivot_swap',
            stepIndex: steps.length,
            matrix: JSON.parse(JSON.stringify(M)),
            pivotRow: k,
            pivotCol: k,
            activeRow: pivotRow,
            activeCol: k,
            multiplier: 0,
            description: `Pivoteamento: Troca L_${k + 1} ↔ L_${pivotRow + 1} (máx: ${maxVal.toFixed(2)})`,
            line: `M([${k + 1}, ${pivotRow + 1}], :) = M([${pivotRow + 1}, ${k + 1}], :);`
          });
        } else {
          steps.push({
            step: steps.length,
            type: 'pivot_selected',
            stepIndex: steps.length,
            matrix: JSON.parse(JSON.stringify(M)),
            pivotRow: k,
            pivotCol: k,
            activeRow: k,
            activeCol: k,
            multiplier: 0,
            description: `Pivô selecionado na coluna ${k + 1}: ${M[k][k].toFixed(2)}`,
            line: `pivo = M(${k + 1}, ${k + 1});`
          });
        }

        const pivo = M[k][k];
        if (Math.abs(pivo) < 1e-14) {
          throw new Error(`Matriz singular ou quase singular no pivô (${k}, ${k})`);
        }

        // 2. Eliminação dos termos abaixo do pivô
        for (let i = k + 1; i < n; i++) {
          const m_ik = M[i][k] / pivo;
          for (let j = k; j <= n; j++) {
            M[i][j] -= m_ik * M[k][j];
          }
          M[i][k] = 0.0; // Supressão analítica de resíduos de ponto flutuante IEEE 754

          steps.push({
            step: steps.length,
            type: 'eliminate',
            stepIndex: steps.length,
            matrix: JSON.parse(JSON.stringify(M)),
            pivotRow: k,
            pivotCol: k,
            activeRow: i,
            activeCol: k,
            multiplier: m_ik,
            description: `Zera M[${i + 1}, ${k + 1}]: L_${i + 1} ← L_${i + 1} - (${m_ik.toFixed(4)}) · L_${k + 1}`,
            line: `m = M(${i + 1}, ${k + 1}) / pivo; M(${i + 1}, :) -= m * M(${k + 1}, :);`
          });
        }
      }

      // 3. Retrossubstituição
      const x = new Array(n).fill(0);
      for (let i = n - 1; i >= 0; i--) {
        let soma = 0;
        for (let j = i + 1; j < n; j++) {
          soma += M[i][j] * x[j];
        }
        x[i] = (M[i][n] - soma) / M[i][i];

        steps.push({
          step: steps.length,
          type: 'back_sub',
          stepIndex: steps.length,
          matrix: JSON.parse(JSON.stringify(M)),
          pivotRow: i,
          pivotCol: i,
          activeRow: i,
          activeCol: n,
          currentSolution: [...x],
          description: `Retrossubstituição x_${i + 1} = (${M[i][n].toFixed(2)} - ${soma.toFixed(2)}) / ${M[i][i].toFixed(2)} = ${x[i].toFixed(4)}`,
          line: `x(${i + 1}) = (M(${i + 1}, end) - soma) / M(${i + 1}, ${i + 1});`
        });
      }

      // 4. Norma Euclidiana do Resíduo ||r||_2 = ||Ax - b||_2
      const residual = A.map((row, i) => row.reduce((sum, aij, j) => sum + aij * x[j], 0) - b[i]);
      const residualNorm = Math.sqrt(residual.reduce((sum, ri) => sum + ri * ri, 0));

      return {
        key: 'gauss_pivoteamento',
        title: 'ELIMINAÇÃO DE GAUSS C/ PIVOTEAMENTO',
        category: 'Sistemas de Equações Lineares',
        file: 'sistemas_lineares/gauss_pivoteamento.m',
        order: 'O(2/3 n³)',
        stability: 'MÁXIMA (Pivoteamento Parcial |m_ik| ≤ 1)',
        formula: 'L_i \\leftarrow L_i - m_{ik} L_k, \\quad m_{ik} = \\frac{a_{ik}}{a_{kk}}',
        formulaSub: 'Triangularização Superior e Retrossubstituição',
        solution: x,
        residual,
        residualNorm,
        n,
        steps
      };
    },

    /**
     * Método de Gauss-Jordan com Normalização Unitária Diagonal
     */
    solveGaussJordan(A_in, b_in) {
      // Sistema do repositório
      const A = A_in ? JSON.parse(JSON.stringify(A_in)) : [
        [55, 0, 1],
        [75, 1, -1],
        [60, -1, 0]
      ];
      const b = b_in ? [...b_in] : [499.55, 687.75, 528.60];
      const n = A.length;

      let M = A.map((row, i) => [...row, b[i]]);
      const steps = [];

      steps.push({
        step: 0,
        type: 'initial',
        stepIndex: 0,
        matrix: JSON.parse(JSON.stringify(M)),
        pivotRow: -1,
        pivotCol: -1,
        activeRow: -1,
        activeCol: -1,
        multiplier: 0,
        description: 'Matriz Aumentada Inicial [A | b]',
        line: 'M = [A, b];'
      });

      for (let k = 0; k < n; k++) {
        // Pivoteamento parcial
        let maxVal = Math.abs(M[k][k]);
        let pivotRow = k;
        for (let p = k + 1; p < n; p++) {
          if (Math.abs(M[p][k]) > maxVal) {
            maxVal = Math.abs(M[p][k]);
            pivotRow = p;
          }
        }

        if (pivotRow !== k) {
          const temp = M[k];
          M[k] = M[pivotRow];
          M[pivotRow] = temp;

          steps.push({
            step: steps.length,
            type: 'pivot_swap',
            stepIndex: steps.length,
            matrix: JSON.parse(JSON.stringify(M)),
            pivotRow: k,
            pivotCol: k,
            activeRow: pivotRow,
            activeCol: k,
            multiplier: 0,
            description: `Pivoteamento: Troca L_${k + 1} ↔ L_${pivotRow + 1}`,
            line: `M([${k + 1}, ${pivotRow + 1}], :) = M([${pivotRow + 1}, ${k + 1}], :);`
          });
        }

        // Normalização da linha pivô para ter 1.0 na diagonal
        const pivo = M[k][k];
        for (let j = k; j <= n; j++) {
          M[k][j] /= pivo;
        }

        steps.push({
          step: steps.length,
          type: 'normalize_row',
          stepIndex: steps.length,
          matrix: JSON.parse(JSON.stringify(M)),
          pivotRow: k,
          pivotCol: k,
          activeRow: k,
          activeCol: k,
          multiplier: 1 / pivo,
          description: `Normalização Unitária: L_${k + 1} ← L_${k + 1} / ${pivo.toFixed(2)}`,
          line: `M(${k + 1}, :) = M(${k + 1}, :) / pivo;`
        });

        // Eliminação tanto acima quanto abaixo da diagonal
        for (let i = 0; i < n; i++) {
          if (i !== k) {
            const fator = M[i][k];
            for (let j = k; j <= n; j++) {
              M[i][j] -= fator * M[k][j];
            }
            M[i][k] = 0.0; // Supressão analítica de resíduos IEEE 754

            steps.push({
              step: steps.length,
              type: 'eliminate',
              stepIndex: steps.length,
              matrix: JSON.parse(JSON.stringify(M)),
              pivotRow: k,
              pivotCol: k,
              activeRow: i,
              activeCol: k,
              multiplier: fator,
              description: `Zera posição (${i + 1}, ${k + 1}): L_${i + 1} ← L_${i + 1} - (${fator.toFixed(4)}) · L_${k + 1}`,
              line: `M(${i + 1}, :) -= M(${i + 1}, ${k + 1}) * M(${k + 1}, :);`
            });
          }
        }
      }

      const x = M.map((row) => row[n]);

      // Norma Euclidiana do Resíduo ||r||_2 = ||Ax - b||_2
      const residual = A.map((row, i) => row.reduce((sum, aij, j) => sum + aij * x[j], 0) - b[i]);
      const residualNorm = Math.sqrt(residual.reduce((sum, ri) => sum + ri * ri, 0));

      steps.push({
        step: steps.length,
        type: 'final_identity',
        stepIndex: steps.length,
        matrix: JSON.parse(JSON.stringify(M)),
        pivotRow: -1,
        pivotCol: -1,
        activeRow: -1,
        activeCol: -1,
        currentSolution: [...x],
        description: `Forma [I | x] Atingida! Solução desacoplada direta.`,
        line: 'x = M(:, end); % Leitura imediata da solução'
      });

      return {
        key: 'gauss_jordan',
        title: 'MÉTODO DE GAUSS-JORDAN',
        category: 'Sistemas de Equações Lineares',
        file: 'sistemas_lineares/gauss_jordan.m',
        order: 'O(n³)',
        stability: 'DIAGONALIZAÇÃO UNITÁRIA [I | x]',
        formula: '[A \\mid \\mathbf{b}] \\xrightarrow{\\text{elim}} [I_n \\mid \\mathbf{x}]',
        formulaSub: 'Eliminação bilateral (acima e abaixo)',
        solution: x,
        residual,
        residualNorm,
        n,
        steps
      };
    }
  };

  // =========================================================================
  // 5. FAMÍLIA 3: QUADRATURA NUMÉRICA (Newton-Cotes)
  // =========================================================================
  const QuadratureSolvers = {
    /**
     * Regra do Trapézio Composta:
     * I ≈ (h/2) * [f(x0) + 2 * sum(f(xi)) + f(xn)]
     */
    solveTrapezioComposto(params = {}) {
      const f = params.f || ((x) => Math.sin(x) + 1.2);
      const a = params.a !== undefined ? params.a : 0.0;
      const b = params.b !== undefined ? params.b : Math.PI;
      const n = params.n || 8;
      const h = (b - a) / n;

      const steps = [];
      let cumulativeArea = 0;

      for (let i = 0; i < n; i++) {
        const x0 = a + i * h;
        const x1 = a + (i + 1) * h;
        const y0 = f(x0);
        const y1 = f(x1);
        const sliceArea = (h / 2.0) * (y0 + y1);
        cumulativeArea += sliceArea;

        steps.push({
          step: i,
          i,
          x0,
          x1,
          y0,
          y1,
          h,
          sliceArea,
          cumulativeArea,
          status: i === n - 1 ? 'CONVERGED' : 'STEPPING',
          line: `A_${i + 1} = (h/2)*(f(x_${i}) + f(x_${i + 1})); % Trapézio ${i + 1}`
        });
      }

      return {
        key: 'trapezio_composto',
        title: 'REGRA DO TRAPÉZIO COMPOSTA',
        category: 'Integração Numérica',
        file: 'integracao_numerica/trapezio_composto.m',
        order: 'O(h²) SEGUNDA ORDEM',
        stability: 'INCONDICIONAL PARA f CONTÍNUA',
        formula: 'I_{trap} = \\frac{h}{2} \\left[ f(x_0) + 2 \\sum f(x_k) + f(x_n) \\right]',
        formulaSub: 'Aproximação por polinômios de 1º grau',
        funcLabel: 'f(x) = sin(x) + 1.2 // [0, π]',
        a,
        b,
        n,
        h,
        totalArea: cumulativeArea,
        f,
        steps
      };
    },

    /**
     * Regra de Simpson 1/3 Composta:
     * I ≈ (h/3) * [f(x0) + 4 * sum(ímpares) + 2 * sum(pares) + f(xn)]
     * Requer obrigatoriamente n par!
     */
    solveSimpson13Composto(params = {}) {
      const f = params.f || ((x) => Math.sin(x) + 1.2);
      const a = params.a !== undefined ? params.a : 0.0;
      const b = params.b !== undefined ? params.b : Math.PI;
      let n = params.n || 8;
      if (n % 2 !== 0) n += 1; // Garante paridade
      const h = (b - a) / n;

      const steps = [];
      let cumulativeArea = 0;
      const numDoubleSlices = n / 2;

      for (let j = 0; j < numDoubleSlices; j++) {
        const x0 = a + (2 * j) * h;
        const x1 = a + (2 * j + 1) * h;
        const x2 = a + (2 * j + 2) * h;
        const y0 = f(x0);
        const y1 = f(x1);
        const y2 = f(x2);

        // Área da parábola de Simpson sobre o subintervalo duplo
        const sliceArea = (h / 3.0) * (y0 + 4.0 * y1 + y2);
        cumulativeArea += sliceArea;

        steps.push({
          step: j,
          j,
          x0,
          x1,
          x2,
          y0,
          y1,
          y2,
          h,
          sliceArea,
          cumulativeArea,
          status: j === numDoubleSlices - 1 ? 'CONVERGED' : 'STEPPING',
          line: `S_${j + 1} = (h/3)*(y0 + 4*y1 + y2); % Parábola ${j + 1}`
        });
      }

      return {
        key: 'simpson_13_composto',
        title: 'REGRA DE SIMPSON 1/3 COMPOSTA',
        category: 'Integração Numérica',
        file: 'integracao_numerica/simpson_13_composto.m',
        order: 'O(h⁴) QUARTA ORDEM',
        stability: 'MÁXIMA EXATIDÃO (n par)',
        formula: 'I_{simp} = \\frac{h}{3} \\left[ f_0 + 4 \\sum_{ímpar} f_k + 2 \\sum_{par} f_k + f_n \\right]',
        formulaSub: 'Interpolação por parábolas de Lagrange (grau 2)',
        funcLabel: 'f(x) = sin(x) + 1.2 // n par = 8',
        a,
        b,
        n,
        h,
        totalArea: cumulativeArea,
        f,
        steps
      };
    }
  };

  // =========================================================================
  // 6. FAMÍLIA 4: EQUAÇÕES DIFERENCIAIS ORDINÁRIAS (PVI)
  // =========================================================================
  const ODESolvers = {
    /**
     * Euler Explícito: y_{i+1} = y_i + h * f(x_i, y_i)
     */
    solveEulerExplicito(params = {}) {
      // PVI clássico de Chapra validado do repositório
      // y' = -2x³ + 12x² - 20x + 8.5, com y(0) = 1.0
      const f = params.f || ((x, y) => -2.0 * Math.pow(x, 3) + 12.0 * Math.pow(x, 2) - 20.0 * x + 8.5);
      const exact = params.exact || ((x) => -0.5 * Math.pow(x, 4) + 4.0 * Math.pow(x, 3) - 10.0 * Math.pow(x, 2) + 8.5 * x + 1.0);
      const x0 = params.x0 !== undefined ? params.x0 : 0.0;
      const xf = params.xf !== undefined ? params.xf : 4.0;
      const y0 = params.y0 !== undefined ? params.y0 : 1.0;
      const h = params.h !== undefined ? params.h : 0.5;

      const steps = [];
      const numSteps = Math.round((xf - x0) / h);
      let xCurr = x0;
      let yCurr = y0;

      for (let i = 0; i < numSteps; i++) {
        const slope = f(xCurr, yCurr);
        const xNext = xCurr + h;
        const yNext = yCurr + h * slope;
        const yExactNext = exact(xNext);
        const yExactCurr = exact(xCurr);
        const errAbs = Math.abs(yExactNext - yNext);
        const errRel = Math.abs(errAbs / (yExactNext || 1)) * 100;

        steps.push({
          step: i,
          i,
          x: xCurr,
          y: yCurr,
          yExact: yExactCurr,
          slope,
          xNext,
          yNext,
          yExactNext,
          errAbs,
          errRel,
          status: i === numSteps - 1 ? 'CONVERGED' : 'STEPPING',
          line: `y(i+1) = y(i) + h * f(x(i), y(i)); % Inclinacao = ${slope.toFixed(4)}`
        });

        xCurr = xNext;
        yCurr = yNext;
      }

      return {
        key: 'euler_explicito',
        title: 'MÉTODO DE EULER EXPLÍCITO',
        category: 'Equações Diferenciais Ordinárias',
        file: 'equacoes_diferenciais/euler_explicito.m',
        order: 'O(h) PRIMEIRA ORDEM (Local O(h²))',
        stability: 'CONDICIONALMENTE ESTÁVEL',
        formula: 'y_{i+1} = y_i + h \\cdot f(x_i, y_i)',
        formulaSub: 'Passo pela reta tangente inicial',
        funcLabel: 'y\' = -2x³ + 12x² - 20x + 8.5 // y(0) = 1.0',
        xRange: [x0, xf],
        yRange: [-1.0, 9.0],
        f,
        exact,
        h,
        steps
      };
    },

    /**
     * Euler Implícito com Newton-Raphson interno:
     * y_{i+1} = y_i + h * f(x_{i+1}, y_{i+1})
     * A-Estável, imune a oscilações em problemas rígidos (stiff)
     */
    solveEulerImplicito(params = {}) {
      // PVI do repositório: y' = 5(1 - y²), y(0) = 0, analítico: tanh(5x)
      // Ou EDO de decaimento rígido: y' = -5(y - x) + 1, y(0) = 0
      const f = params.f || ((x, y) => -5.0 * (y - 1.0));
      const df_dy = params.df_dy || ((x, y) => -5.0);
      const exact = params.exact || ((x) => 1.0 - Math.exp(-5.0 * x));
      const x0 = params.x0 !== undefined ? params.x0 : 0.0;
      const xf = params.xf !== undefined ? params.xf : 1.2;
      const y0 = params.y0 !== undefined ? params.y0 : 0.0;
      const h = params.h !== undefined ? params.h : 0.15;

      const steps = [];
      const numSteps = Math.round((xf - x0) / h);
      let xCurr = x0;
      let yCurr = y0;

      for (let i = 0; i < numSteps; i++) {
        const xNext = xCurr + h;
        // Chute inicial para Newton: Euler explícito u0 = yCurr + h * f(xCurr, yCurr)
        let u = yCurr + h * f(xCurr, yCurr);
        const newtonHistory = [];

        // Loop interno de Newton-Raphson para resolver G(u) = u - yCurr - h*f(xNext, u) = 0
        for (let it = 0; it < 10; it++) {
          const Gu = u - yCurr - h * f(xNext, u);
          const dGu = 1.0 - h * df_dy(xNext, u);
          const deltaU = Gu / dGu;
          u = u - deltaU;
          newtonHistory.push({ it, u, Gu, deltaU });
          if (Math.abs(Gu) < 1e-9 || Math.abs(deltaU) < 1e-9) break;
        }

        const yNext = u;
        const yExactNext = exact(xNext);
        const yExactCurr = exact(xCurr);
        const errAbs = Math.abs(yExactNext - yNext);
        const errRel = Math.abs(errAbs / (yExactNext || 1)) * 100;

        steps.push({
          step: i,
          i,
          x: xCurr,
          y: yCurr,
          yExact: yExactCurr,
          xNext,
          yNext,
          yExactNext,
          errAbs,
          errRel,
          newtonIterations: newtonHistory.length,
          status: i === numSteps - 1 ? 'CONVERGED' : 'STEPPING',
          line: `u = u - G(u)/G'(u); % Newton interno (${newtonHistory.length} it)`
        });

        xCurr = xNext;
        yCurr = yNext;
      }

      return {
        key: 'euler_implicito',
        title: 'EULER IMPLÍCITO C/ NEWTON-RAPHSON',
        category: 'Equações Diferenciais Ordinárias',
        file: 'equacoes_diferenciais/euler_implicito_newton.m',
        order: 'O(h) A-ESTÁVEL (Incondicional)',
        stability: 'INCONDICIONALMENTE ESTÁVEL (Stiff)',
        formula: 'y_{i+1} = y_i + h \\cdot f(x_{i+1}, y_{i+1})',
        formulaSub: 'Solução não linear por Newton-Raphson',
        funcLabel: 'y\' = -5(y - 1) // y(0) = 0',
        xRange: [x0, xf],
        yRange: [-0.1, 1.2],
        f,
        exact,
        h,
        steps
      };
    },

    /**
     * Runge-Kutta de 2ª Ordem / Método de Heun (preditor-corretor trapezoidal):
     * k1 = f(x_i, y_i)
     * y_pred = y_i + h * k1
     * k2 = f(x_i + h, y_pred)
     * y_{i+1} = y_i + (h/2)*(k1 + k2)
     */
    solveRK2Heun(params = {}) {
      // PVI do repositório: y' = -2y + 4x, y(0) = 1.0, h = 0.2
      // Solução analítica exata: y(x) = 2x - 1 + 2*exp(-2x)
      const f = params.f || ((x, y) => -2.0 * y + 4.0 * x);
      const exact = params.exact || ((x) => 2.0 * x - 1.0 + 2.0 * Math.exp(-2.0 * x));
      const x0 = params.x0 !== undefined ? params.x0 : 0.0;
      const xf = params.xf !== undefined ? params.xf : 2.0;
      const y0 = params.y0 !== undefined ? params.y0 : 1.0;
      const h = params.h !== undefined ? params.h : 0.2;

      const steps = [];
      const numSteps = Math.round((xf - x0) / h);
      let xCurr = x0;
      let yCurr = y0;

      for (let i = 0; i < numSteps; i++) {
        const k1 = f(xCurr, yCurr);
        const yPred = yCurr + h * k1;
        const xNext = xCurr + h;
        const k2 = f(xNext, yPred);
        const slopeAvg = (k1 + k2) / 2.0;
        const yNext = yCurr + h * slopeAvg;
        const yExactNext = exact(xNext);
        const yExactCurr = exact(xCurr);
        const errAbs = Math.abs(yExactNext - yNext);
        const errRel = Math.abs(errAbs / (yExactNext || 1)) * 100;

        steps.push({
          step: i,
          i,
          x: xCurr,
          y: yCurr,
          yExact: yExactCurr,
          k1,
          yPred,
          k2,
          slopeAvg,
          xNext,
          yNext,
          yExactNext,
          errAbs,
          errRel,
          status: i === numSteps - 1 ? 'CONVERGED' : 'STEPPING',
          line: `y(i+1) = y(i) + (h/2)*(k1 + k2); % Heun preditor-corretor`
        });

        xCurr = xNext;
        yCurr = yNext;
      }

      return {
        key: 'rk2_heun',
        title: 'RUNGE-KUTTA 2ª ORDEM (HEUN)',
        category: 'Equações Diferenciais Ordinárias',
        file: 'equacoes_diferenciais/rk2_individual.m',
        order: 'O(h²) SEGUNDA ORDEM (Local O(h³))',
        stability: 'PREDITOR-CORRETOR TRAPEZOIDAL',
        formula: 'y_{i+1} = y_i + \\frac{h}{2} \\left[ f(x_i, y_i) + f(x_{i+1}, y_{pred}) \\right]',
        formulaSub: 'Média de inclinações nas extremidades',
        funcLabel: 'y\' = -2y + 4x // y(0) = 1.0',
        xRange: [x0, xf],
        yRange: [0.0, 4.0],
        f,
        exact,
        h,
        steps
      };
    }
  };

  // =========================================================================
  // 7. REGISTRO & MANIFESTO DE MÉTODOS
  // =========================================================================
  const MethodsManifest = {
    // 01 // ZEROS DE FUNÇÕES
    newton: {
      key: 'newton',
      code: 'ALG.01',
      name: 'NEWTON-RAPHSON',
      family: 'scalar_root',
      file: 'raizes_equacoes/newton_raphson.m',
      solver: ScalarRootSolvers.solveNewton
    },
    bissecao: {
      key: 'bissecao',
      code: 'ALG.02',
      name: 'MÉTODO DA BISSEÇÃO',
      family: 'scalar_root',
      file: 'raizes_equacoes/bissecao.m',
      solver: ScalarRootSolvers.solveBissecao
    },
    secante: {
      key: 'secante',
      code: 'ALG.03',
      name: 'MÉTODO DA SECANTE',
      family: 'scalar_root',
      file: 'raizes_equacoes/secante.m',
      solver: ScalarRootSolvers.solveSecante
    },
    falsa_posicao: {
      key: 'falsa_posicao',
      code: 'ALG.04',
      name: 'FALSA POSIÇÃO',
      family: 'scalar_root',
      file: 'raizes_equacoes/falsa_posicao.m',
      solver: ScalarRootSolvers.solveFalsaPosicao
    },
    ponto_fixo: {
      key: 'ponto_fixo',
      code: 'ALG.05',
      name: 'PONTO FIXO (PICARD-BANACH)',
      family: 'scalar_root',
      file: 'raizes_equacoes/ponto_fixo.m',
      solver: ScalarRootSolvers.solvePontoFixo
    },

    // 02 // SISTEMAS LINEARES
    gauss_pivoteamento: {
      key: 'gauss_pivoteamento',
      code: 'ALG.06',
      name: 'GAUSS C/ PIVOTEAMENTO',
      family: 'linear_system',
      file: 'sistemas_lineares/gauss_pivoteamento.m',
      solver: LinearSystemSolvers.solveGaussPivoteamento
    },
    gauss_jordan: {
      key: 'gauss_jordan',
      code: 'ALG.07',
      name: 'MÉTODO DE GAUSS-JORDAN',
      family: 'linear_system',
      file: 'sistemas_lineares/gauss_jordan.m',
      solver: LinearSystemSolvers.solveGaussJordan
    },

    // 03 // INTEGRAÇÃO NUMÉRICA
    trapezio_composto: {
      key: 'trapezio_composto',
      code: 'ALG.08',
      name: 'TRAPÉZIO COMPOSTO',
      family: 'quadrature',
      file: 'integracao_numerica/trapezio_composto.m',
      solver: QuadratureSolvers.solveTrapezioComposto
    },
    simpson_13_composto: {
      key: 'simpson_13_composto',
      code: 'ALG.09',
      name: 'SIMPSON 1/3 COMPOSTO',
      family: 'quadrature',
      file: 'integracao_numerica/simpson_13_composto.m',
      solver: QuadratureSolvers.solveSimpson13Composto
    },

    // 04 // EQUAÇÕES DIFERENCIAIS
    euler_explicito: {
      key: 'euler_explicito',
      code: 'ALG.10',
      name: 'EULER EXPLÍCITO',
      family: 'ode',
      file: 'equacoes_diferenciais/euler_explicito.m',
      solver: ODESolvers.solveEulerExplicito
    },
    euler_implicito: {
      key: 'euler_implicito',
      code: 'ALG.11',
      name: 'EULER IMPLÍCITO (NEWTON)',
      family: 'ode',
      file: 'equacoes_diferenciais/euler_implicito_newton.m',
      solver: ODESolvers.solveEulerImplicito
    },
    rk2_heun: {
      key: 'rk2_heun',
      code: 'ALG.12',
      name: 'RUNGE-KUTTA 2 (HEUN)',
      family: 'ode',
      file: 'equacoes_diferenciais/rk2_individual.m',
      solver: ODESolvers.solveRK2Heun
    }
  };

  // =========================================================================
  // 8. RENDERIZADOR CRT CANVAS VETORIAL
  // =========================================================================
  const CRTRenderer = {
    /**
     * Limpa a tela e desenha o retículo de medição do osciloscópio
     */
    drawGrid(ctx, width, height) {
      ctx.save();
      ctx.clearRect(0, 0, width, height);

      // Fundo escuro com gradiente radial de fósforo
      const grad = ctx.createRadialGradient(
        width * 0.45, height * 0.5, width * 0.05,
        width * 0.5, height * 0.5, width * 0.75
      );
      grad.addColorStop(0, '#0a1c2e');
      grad.addColorStop(0.7, '#020509');
      grad.addColorStop(1, '#000204');
      ctx.fillStyle = grad;
      ctx.fillRect(0, 0, width, height);

      // Linhas da grade reticular
      ctx.strokeStyle = CRTTheme.reticleLine;
      ctx.lineWidth = 1;
      const step = 35;

      ctx.beginPath();
      for (let x = step; x < width; x += step) {
        ctx.moveTo(x, 0);
        ctx.lineTo(x, height);
      }
      for (let y = step; y < height; y += step) {
        ctx.moveTo(0, y);
        ctx.lineTo(width, y);
      }
      ctx.stroke();

      ctx.restore();
    },

    /**
     * Eixos táticos cartesianos
     */
    drawAxes(ctx, originX, originY, width, height) {
      ctx.save();
      ctx.strokeStyle = CRTTheme.axisLine;
      ctx.lineWidth = 2.0;

      ctx.beginPath();
      // Eixo X
      ctx.moveTo(25, originY);
      ctx.lineTo(width - 25, originY);
      // Eixo Y
      ctx.moveTo(originX, 20);
      ctx.lineTo(originX, height - 20);
      ctx.stroke();

      // Setas dos eixos
      ctx.fillStyle = CRTTheme.axisLine;
      ctx.beginPath();
      ctx.moveTo(width - 20, originY);
      ctx.lineTo(width - 28, originY - 5);
      ctx.lineTo(width - 28, originY + 5);
      ctx.fill();

      ctx.beginPath();
      ctx.moveTo(originX, 15);
      ctx.lineTo(originX - 5, 23);
      ctx.lineTo(originX + 5, 23);
      ctx.fill();

      // Rótulos de eixos
      ctx.font = 'bold 13px "JetBrains Mono", monospace';
      ctx.fillStyle = CRTTheme.textDim;
      ctx.fillText('x', width - 20, originY + 22);
      ctx.fillText('y', originX + 10, 24);

      ctx.restore();
    },

    /**
     * Efeito sutil de scanline CRT
     */
    drawScanlines(ctx, width, height) {
      ctx.save();
      ctx.fillStyle = 'rgba(0, 0, 0, 0.15)';
      for (let y = 0; y < height; y += 4) {
        ctx.fillRect(0, y, width, 1.5);
      }
      ctx.restore();
    },

    // -----------------------------------------------------------------------
    // RENDERIZAÇÃO: RAÍZES DE EQUAÇÕES ESCALARES
    // -----------------------------------------------------------------------
    renderScalarRoot(ctx, methodData, stepIndex, animProgress, canvasWidth, canvasHeight) {
      this.drawGrid(ctx, canvasWidth, canvasHeight);

      const margin = 50;
      const xMin = methodData.xRange[0];
      const xMax = methodData.xRange[1];
      const yMin = methodData.yRange[0];
      const yMax = methodData.yRange[1];

      const toScreenX = (x) => MathUtils.mapRange(x, xMin, xMax, margin, canvasWidth - margin);
      const toScreenY = (y) => MathUtils.mapRange(y, yMin, yMax, canvasHeight - margin, margin);

      const originX = toScreenX(0);
      const originY = toScreenY(0);

      this.drawAxes(ctx, MathUtils.clamp(originX, margin, canvasWidth - margin), originY, canvasWidth, canvasHeight);

      // Curva Analítica f(x)
      ctx.save();
      ctx.strokeStyle = CRTTheme.cyan;
      ctx.lineWidth = 3.5;
      ctx.shadowColor = CRTTheme.cyan;
      ctx.shadowBlur = 12;

      ctx.beginPath();
      const numPts = 120;
      for (let i = 0; i <= numPts; i++) {
        const px = xMin + (i / numPts) * (xMax - xMin);
        const py = methodData.f(px);
        const sx = toScreenX(px);
        const sy = toScreenY(py);
        if (i === 0) ctx.moveTo(sx, sy);
        else ctx.lineTo(sx, sy);
      }
      ctx.stroke();
      ctx.restore();

      // Ponto Fixo: desenha também a reta identidade y = x e curva g(x)
      if (methodData.key === 'ponto_fixo' && methodData.g) {
        ctx.save();
        ctx.strokeStyle = '#5a738e';
        ctx.lineWidth = 1.5;
        ctx.setLineDash([4, 4]);
        ctx.beginPath();
        ctx.moveTo(toScreenX(xMin), toScreenY(xMin));
        ctx.lineTo(toScreenX(xMax), toScreenY(xMax));
        ctx.stroke();

        ctx.strokeStyle = CRTTheme.amber;
        ctx.lineWidth = 2.5;
        ctx.setLineDash([]);
        ctx.beginPath();
        for (let i = 0; i <= numPts; i++) {
          const px = xMin + (i / numPts) * (xMax - xMin);
          const py = methodData.g(px);
          const sx = toScreenX(px);
          const sy = toScreenY(py);
          if (i === 0) ctx.moveTo(sx, sy);
          else ctx.lineTo(sx, sy);
        }
        ctx.stroke();
        ctx.restore();
      }

      // Alvo Raiz Teórica x*
      if (methodData.target !== undefined) {
        const targetX = toScreenX(methodData.target);
        ctx.save();
        const pulse = 1.0 + 0.25 * Math.sin(Date.now() * 0.006);
        ctx.fillStyle = CRTTheme.green;
        ctx.shadowColor = CRTTheme.green;
        ctx.shadowBlur = 14 * pulse;
        ctx.beginPath();
        ctx.arc(targetX, originY, 6.5 * pulse, 0, Math.PI * 2);
        ctx.fill();

        ctx.font = 'bold 14px "JetBrains Mono", monospace';
        ctx.fillText(`ALVO x* [≈ ${methodData.target.toFixed(4)}]`, targetX - 40, originY - 18);
        ctx.restore();
      }

      // Detalhes do passo ativo
      const cur = methodData.steps[stepIndex];
      if (!cur) {
        this.drawScanlines(ctx, canvasWidth, canvasHeight);
        return;
      }

      ctx.save();
      if (methodData.key === 'newton') {
        const sx0 = toScreenX(cur.x);
        const sy0 = originY;
        const sfx0 = toScreenY(cur.fx);
        const sxNext = toScreenX(cur.xNext);

        // Ponto x_k
        ctx.fillStyle = CRTTheme.amber;
        ctx.shadowColor = CRTTheme.amber;
        ctx.shadowBlur = 10;
        ctx.beginPath();
        ctx.arc(sx0, sy0, 7, 0, Math.PI * 2);
        ctx.fill();

        ctx.font = 'bold 15px "JetBrains Mono", monospace';
        ctx.fillText(`x_${cur.k}`, sx0 - 12, sy0 + 24);

        // Feixe vertical até f(x_k)
        const vertProg = Math.min(1.0, animProgress * 2.0);
        ctx.strokeStyle = CRTTheme.amber;
        ctx.lineWidth = 2.0;
        ctx.setLineDash([5, 5]);
        ctx.beginPath();
        ctx.moveTo(sx0, sy0);
        ctx.lineTo(sx0, sy0 + (sfx0 - sy0) * vertProg);
        ctx.stroke();
        ctx.setLineDash([]);

        if (vertProg >= 0.8) {
          ctx.fillStyle = '#ffffff';
          ctx.beginPath();
          ctx.arc(sx0, sfx0, 6, 0, Math.PI * 2);
          ctx.fill();
        }

        // Reta tangente descendo até x_{k+1}
        if (animProgress > 0.4) {
          const tangProg = Math.min(1.0, (animProgress - 0.4) / 0.6);
          ctx.strokeStyle = CRTTheme.red;
          ctx.lineWidth = 3.0;
          ctx.shadowColor = CRTTheme.red;
          ctx.shadowBlur = 14;
          ctx.beginPath();
          ctx.moveTo(sx0, sfx0);
          ctx.lineTo(sx0 + (sxNext - sx0) * tangProg, sfx0 + (sy0 - sfx0) * tangProg);
          ctx.stroke();

          if (tangProg >= 0.9) {
            ctx.fillStyle = CRTTheme.green;
            ctx.shadowColor = CRTTheme.green;
            ctx.shadowBlur = 12;
            ctx.beginPath();
            ctx.arc(sxNext, sy0, 7.5, 0, Math.PI * 2);
            ctx.fill();
            ctx.fillText(`x_${cur.k + 1}`, sxNext - 12, sy0 - 16);
          }
        }
      } else if (methodData.key === 'bissecao' || methodData.key === 'falsa_posicao') {
        const sa = toScreenX(cur.a);
        const sb = toScreenX(cur.b);
        const sxm = toScreenX(methodData.key === 'bissecao' ? cur.xm : cur.xr);

        // Faixa sombreada do intervalo ativo [a, b]
        ctx.fillStyle = 'rgba(0, 240, 255, 0.15)';
        ctx.fillRect(sa, margin, sb - sa, canvasHeight - 2 * margin);

        // Linhas de fronteira a e b
        ctx.strokeStyle = CRTTheme.cyan;
        ctx.lineWidth = 2.5;
        ctx.beginPath();
        ctx.moveTo(sa, margin);
        ctx.lineTo(sa, canvasHeight - margin);
        ctx.moveTo(sb, margin);
        ctx.lineTo(sb, canvasHeight - margin);
        ctx.stroke();

        ctx.fillStyle = CRTTheme.cyan;
        ctx.font = 'bold 14px "JetBrains Mono", monospace';
        const isNarrow = Math.abs(sb - sa) < 70;
        ctx.fillText(`a=${cur.a.toFixed(3)}`, sa - 25, isNarrow ? (margin - 22) : (margin - 10));
        ctx.fillText(`b=${cur.b.toFixed(3)}`, sb - 25, isNarrow ? (margin - 6) : (margin - 10));

        // Falsa Posição: desenha a reta secante unindo (a, f(a)) a (b, f(b))
        if (methodData.key === 'falsa_posicao') {
          const sfa = toScreenY(cur.fa);
          const sfb = toScreenY(cur.fb);
          ctx.strokeStyle = CRTTheme.red;
          ctx.lineWidth = 2.5;
          ctx.shadowColor = CRTTheme.red;
          ctx.shadowBlur = 10;
          ctx.beginPath();
          ctx.moveTo(sa, sfa);
          ctx.lineTo(sb, sfb);
          ctx.stroke();
        }

        // Ponto de corte x_m ou x_r
        if (animProgress > 0.25) {
          ctx.strokeStyle = CRTTheme.amber;
          ctx.lineWidth = 3.0;
          ctx.beginPath();
          ctx.moveTo(sxm, margin + 10);
          ctx.lineTo(sxm, canvasHeight - margin - 10);
          ctx.stroke();

          ctx.fillStyle = CRTTheme.amber;
          const label = methodData.key === 'bissecao' ? `x_m = ${cur.xm.toFixed(4)}` : `x_r = ${cur.xr.toFixed(4)}`;
          ctx.fillText(label, sxm - 30, canvasHeight - margin + 22);
        }
      } else if (methodData.key === 'secante') {
        const sx0 = toScreenX(cur.xPrev);
        const sfx0 = toScreenY(cur.fPrev);
        const sx1 = toScreenX(cur.x);
        const sfx1 = toScreenY(cur.fx);
        const sxNext = toScreenX(cur.xNext);

        // Pontos anteriores
        ctx.fillStyle = CRTTheme.amber;
        ctx.beginPath();
        ctx.arc(sx0, sfx0, 6, 0, Math.PI * 2);
        ctx.arc(sx1, sfx1, 6, 0, Math.PI * 2);
        ctx.fill();

        // Reta secante cruzando y = 0
        ctx.strokeStyle = CRTTheme.red;
        ctx.lineWidth = 2.5;
        ctx.shadowColor = CRTTheme.red;
        ctx.shadowBlur = 12;
        ctx.beginPath();
        ctx.moveTo(sx0, sfx0);
        ctx.lineTo(sxNext, originY);
        ctx.stroke();

        ctx.fillStyle = CRTTheme.green;
        ctx.beginPath();
        ctx.arc(sxNext, originY, 7, 0, Math.PI * 2);
        ctx.fill();
        ctx.fillText(`x_${cur.k + 1}`, sxNext - 12, originY - 14);
      } else if (methodData.key === 'ponto_fixo') {
        const sx = toScreenX(cur.x);
        const sgx = toScreenY(cur.gx);
        const sxNext = toScreenX(cur.xNext);

        // Diagrama em Teia (Cobweb)
        ctx.strokeStyle = CRTTheme.red;
        ctx.lineWidth = 2.5;
        ctx.beginPath();
        // vertical de (x_k, x_k) até (x_k, g(x_k))
        ctx.moveTo(sx, toScreenY(cur.x));
        ctx.lineTo(sx, sgx);
        // horizontal até reta y = x (g(x_k), g(x_k))
        ctx.lineTo(sxNext, sgx);
        ctx.stroke();

        ctx.fillStyle = CRTTheme.green;
        ctx.beginPath();
        ctx.arc(sxNext, sgx, 7, 0, Math.PI * 2);
        ctx.fill();
      }

      ctx.restore();
      this.drawScanlines(ctx, canvasWidth, canvasHeight);
    },

    // -----------------------------------------------------------------------
    // RENDERIZAÇÃO: SISTEMAS LINEARES (TABULEIRO MATRICIAL CRT)
    // -----------------------------------------------------------------------
    renderMatrixElimination(ctx, methodData, stepIndex, animProgress, canvasWidth, canvasHeight) {
      this.drawGrid(ctx, canvasWidth, canvasHeight);

      const cur = methodData.steps[stepIndex];
      if (!cur) {
        this.drawScanlines(ctx, canvasWidth, canvasHeight);
        return;
      }

      const M = cur.matrix;
      const rows = M.length;
      const cols = M[0].length; // n + 1 colunas (última é b)

      ctx.save();

      // Título do Passo / Operação Ativa
      ctx.font = 'bold 18px "JetBrains Mono", monospace';
      ctx.fillStyle = CRTTheme.cyan;
      ctx.shadowColor = CRTTheme.cyan;
      ctx.shadowBlur = 8;
      ctx.fillText(`ETAPA [${stepIndex + 1}/${methodData.steps.length}]: ${cur.description}`, 40, 42);
      ctx.shadowBlur = 0;

      // Código executando
      ctx.font = '14px "JetBrains Mono", monospace';
      ctx.fillStyle = CRTTheme.green;
      ctx.fillText(`LINHA: ${cur.line}`, 40, 70);

      // Caixa da Matriz Aumentada
      const tableWidth = Math.min(canvasWidth - 80, cols * 130);
      const cellWidth = tableWidth / cols;
      const cellHeight = 65;
      const startX = (canvasWidth - tableWidth) / 2;
      const startY = 110;

      // Colchetes externos CRT da matriz aumentada
      ctx.strokeStyle = CRTTheme.cyan;
      ctx.lineWidth = 3.5;
      ctx.beginPath();
      // Colchete esquerdo
      ctx.moveTo(startX - 12, startY - 10);
      ctx.lineTo(startX - 24, startY - 10);
      ctx.lineTo(startX - 24, startY + rows * cellHeight + 10);
      ctx.lineTo(startX - 12, startY + rows * cellHeight + 10);
      // Colchete direito
      const endX = startX + tableWidth;
      ctx.moveTo(endX + 12, startY - 10);
      ctx.lineTo(endX + 24, startY - 10);
      ctx.lineTo(endX + 24, startY + rows * cellHeight + 10);
      ctx.lineTo(endX + 12, startY + rows * cellHeight + 10);
      ctx.stroke();

      // Barra vertical divisória entre [A] e [b]
      const divX = startX + (cols - 1) * cellWidth;
      ctx.strokeStyle = CRTTheme.amber;
      ctx.lineWidth = 2.5;
      ctx.setLineDash([6, 6]);
      ctx.beginPath();
      ctx.moveTo(divX, startY - 5);
      ctx.lineTo(divX, startY + rows * cellHeight + 5);
      ctx.stroke();
      ctx.setLineDash([]);

      // Desenho das Células
      for (let r = 0; r < rows; r++) {
        for (let c = 0; c < cols; c++) {
          const cx = startX + c * cellWidth;
          const cy = startY + r * cellHeight;
          const val = M[r][c];

          const isPivot = r === cur.pivotRow && c === cur.pivotCol;
          const isActiveCell = r === cur.activeRow && c === cur.activeCol;
          const isPivotRow = r === cur.pivotRow;
          const isActiveRow = r === cur.activeRow;

          // Fundo da célula
          if (isPivot) {
            ctx.fillStyle = 'rgba(255, 183, 0, 0.28)';
          } else if (isActiveCell) {
            ctx.fillStyle = 'rgba(0, 255, 102, 0.28)';
          } else if (isPivotRow) {
            ctx.fillStyle = 'rgba(255, 183, 0, 0.08)';
          } else if (isActiveRow) {
            ctx.fillStyle = 'rgba(0, 240, 255, 0.08)';
          } else {
            ctx.fillStyle = 'rgba(9, 18, 29, 0.6)';
          }
          ctx.fillRect(cx + 4, cy + 4, cellWidth - 8, cellHeight - 8);

          // Borda da célula
          if (isPivot) {
            ctx.strokeStyle = CRTTheme.amber;
            ctx.lineWidth = 2.5;
            ctx.strokeRect(cx + 4, cy + 4, cellWidth - 8, cellHeight - 8);
          } else if (isActiveCell) {
            ctx.strokeStyle = CRTTheme.green;
            ctx.lineWidth = 2.5;
            ctx.strokeRect(cx + 4, cy + 4, cellWidth - 8, cellHeight - 8);
          }

          // Rótulo da variável no cabeçalho superior
          if (r === 0) {
            ctx.fillStyle = CRTTheme.textDim;
            ctx.font = '12px "JetBrains Mono", monospace';
            const colLabel = c === cols - 1 ? 'VETOR b' : `x_${c + 1}`;
            ctx.fillText(colLabel, cx + cellWidth * 0.35, startY - 18);
          }

          // Valor numérico
          ctx.font = isPivot || isActiveCell ? 'bold 18px "JetBrains Mono", monospace' : '16px "JetBrains Mono", monospace';
          ctx.fillStyle = isPivot ? CRTTheme.amber : isActiveCell ? CRTTheme.green : CRTTheme.text;
          const textVal = MathUtils.formatFloat(val, 3);
          ctx.fillText(textVal, cx + 18, cy + 38);
        }
      }

      // Painel inferior: exibição da solução e resíduo quando disponível
      if (cur.currentSolution) {
        ctx.fillStyle = CRTTheme.green;
        ctx.font = 'bold 16px "JetBrains Mono", monospace';
        const solText = cur.currentSolution.map((s, idx) => `x_${idx + 1} = ${s.toFixed(4)}`).join('   |   ');
        ctx.fillText(`SOLUÇÃO DO SISTEMA:   ${solText}`, 40, startY + rows * cellHeight + 42);

        const resText = methodData.residualNorm !== undefined ? MathUtils.formatFloat(methodData.residualNorm, 6) : '0.000000e+00';
        ctx.fillStyle = CRTTheme.cyan;
        ctx.font = '13px "JetBrains Mono", monospace';
        ctx.fillText(`RESÍDUO EUCLIDIANO ||Ax - b||₂: ${resText}   |   CONFORMIDADE IEEE 754: NOMINAL`, 40, startY + rows * cellHeight + 66);
      }

      ctx.restore();
      this.drawScanlines(ctx, canvasWidth, canvasHeight);
    },

    // -----------------------------------------------------------------------
    // RENDERIZAÇÃO: QUADRATURA NUMÉRICA (Newton-Cotes)
    // -----------------------------------------------------------------------
    renderQuadrature(ctx, methodData, stepIndex, animProgress, canvasWidth, canvasHeight) {
      this.drawGrid(ctx, canvasWidth, canvasHeight);

      const margin = 55;
      const a = methodData.a;
      const b = methodData.b;
      const xMin = a - 0.2;
      const xMax = b + 0.2;
      const yMin = 0.0;
      const yMax = 3.0;

      const toScreenX = (x) => MathUtils.mapRange(x, xMin, xMax, margin, canvasWidth - margin);
      const toScreenY = (y) => MathUtils.mapRange(y, yMin, yMax, canvasHeight - margin, margin);

      this.drawAxes(ctx, toScreenX(0), toScreenY(0), canvasWidth, canvasHeight);

      // Fatias de Trapézio ou Parábolas de Simpson completadas até o passo atual
      ctx.save();

      if (methodData.key === 'trapezio_composto') {
        for (let i = 0; i <= stepIndex && i < methodData.steps.length; i++) {
          const cur = methodData.steps[i];
          const sx0 = toScreenX(cur.x0);
          const sx1 = toScreenX(cur.x1);
          const sy0 = toScreenY(cur.y0);
          const sy1 = toScreenY(cur.y1);
          const groundY = toScreenY(0);

          const isCurrent = i === stepIndex;

          // Preenchimento poligonal do trapézio
          ctx.fillStyle = isCurrent ? 'rgba(0, 255, 102, 0.35)' : 'rgba(0, 240, 255, 0.2)';
          ctx.beginPath();
          ctx.moveTo(sx0, groundY);
          ctx.lineTo(sx0, sy0);
          ctx.lineTo(sx1, sy1);
          ctx.lineTo(sx1, groundY);
          ctx.closePath();
          ctx.fill();

          // Aresta do teto (segmento de reta linear)
          ctx.strokeStyle = isCurrent ? CRTTheme.green : CRTTheme.cyan;
          ctx.lineWidth = isCurrent ? 3.0 : 2.0;
          ctx.beginPath();
          ctx.moveTo(sx0, sy0);
          ctx.lineTo(sx1, sy1);
          ctx.stroke();

          // Linhas verticais separadoras
          ctx.strokeStyle = 'rgba(26, 50, 77, 0.9)';
          ctx.lineWidth = 1.5;
          ctx.beginPath();
          ctx.moveTo(sx0, groundY);
          ctx.lineTo(sx0, sy0);
          ctx.moveTo(sx1, groundY);
          ctx.lineTo(sx1, sy1);
          ctx.stroke();
        }
      } else if (methodData.key === 'simpson_13_composto') {
        for (let j = 0; j <= stepIndex && j < methodData.steps.length; j++) {
          const cur = methodData.steps[j];
          const sx0 = toScreenX(cur.x0);
          const sx1 = toScreenX(cur.x1);
          const sx2 = toScreenX(cur.x2);
          const groundY = toScreenY(0);

          const isCurrent = j === stepIndex;

          // Interpolação parabólica P(x) = c0 + c1*(x-x0) + c2*(x-x0)*(x-x1)
          const h = cur.h;
          ctx.fillStyle = isCurrent ? 'rgba(255, 183, 0, 0.35)' : 'rgba(0, 240, 255, 0.2)';
          ctx.beginPath();
          ctx.moveTo(sx0, groundY);
          ctx.lineTo(sx0, toScreenY(cur.y0));

          const numSub = 24;
          for (let s = 1; s <= numSub; s++) {
            const px = cur.x0 + (s / numSub) * (2 * h);
            // Avaliação polinomial da parábola pelos pontos (x0, y0), (x1, y1), (x2, y2)
            const L0 = ((px - cur.x1) * (px - cur.x2)) / (2 * h * h);
            const L1 = -((px - cur.x0) * (px - cur.x2)) / (h * h);
            const L2 = ((px - cur.x0) * (px - cur.x1)) / (2 * h * h);
            const py = cur.y0 * L0 + cur.y1 * L1 + cur.y2 * L2;
            ctx.lineTo(toScreenX(px), toScreenY(py));
          }

          ctx.lineTo(sx2, groundY);
          ctx.closePath();
          ctx.fill();

          // Aresta parabólica em destaque
          ctx.strokeStyle = isCurrent ? CRTTheme.amber : CRTTheme.cyan;
          ctx.lineWidth = isCurrent ? 3.5 : 2.0;
          ctx.beginPath();
          for (let s = 0; s <= numSub; s++) {
            const px = cur.x0 + (s / numSub) * (2 * h);
            const L0 = ((px - cur.x1) * (px - cur.x2)) / (2 * h * h);
            const L1 = -((px - cur.x0) * (px - cur.x2)) / (h * h);
            const L2 = ((px - cur.x0) * (px - cur.x1)) / (2 * h * h);
            const py = cur.y0 * L0 + cur.y1 * L1 + cur.y2 * L2;
            if (s === 0) ctx.moveTo(toScreenX(px), toScreenY(py));
            else ctx.lineTo(toScreenX(px), toScreenY(py));
          }
          ctx.stroke();
        }
      }

      // Curva Analítica Real f(x)
      ctx.strokeStyle = '#ffffff';
      ctx.lineWidth = 2.0;
      ctx.setLineDash([4, 4]);
      ctx.beginPath();
      for (let px = a; px <= b; px += 0.05) {
        const sx = toScreenX(px);
        const sy = toScreenY(methodData.f(px));
        if (px === a) ctx.moveTo(sx, sy);
        else ctx.lineTo(sx, sy);
      }
      ctx.stroke();
      ctx.setLineDash([]);

      // Telemetria da Área Acumulada
      const activeStep = methodData.steps[stepIndex];
      if (activeStep) {
        ctx.font = 'bold 16px "JetBrains Mono", monospace';
        ctx.fillStyle = CRTTheme.green;
        ctx.fillText(`ÁREA ACUMULADA: ∫ f(x)dx ≈ ${activeStep.cumulativeArea.toFixed(6)}`, margin + 20, margin + 25);
        ctx.fillStyle = CRTTheme.amber;
        ctx.fillText(`FATIA ATUAL: ΔA = ${activeStep.sliceArea.toFixed(6)}`, margin + 20, margin + 48);
      }

      ctx.restore();
      this.drawScanlines(ctx, canvasWidth, canvasHeight);
    },

    // -----------------------------------------------------------------------
    // RENDERIZAÇÃO: EQUAÇÕES DIFERENCIAIS ORDINÁRIAS (PVI)
    // -----------------------------------------------------------------------
    renderODEStepper(ctx, methodData, stepIndex, animProgress, canvasWidth, canvasHeight) {
      this.drawGrid(ctx, canvasWidth, canvasHeight);

      const margin = 55;
      const xMin = methodData.xRange[0];
      const xMax = methodData.xRange[1];
      const yMin = methodData.yRange[0];
      const yMax = methodData.yRange[1];

      const toScreenX = (x) => MathUtils.mapRange(x, xMin, xMax, margin, canvasWidth - margin);
      const toScreenY = (y) => MathUtils.mapRange(y, yMin, yMax, canvasHeight - margin, margin);

      this.drawAxes(ctx, toScreenX(0), toScreenY(0), canvasWidth, canvasHeight);

      // Curva Analítica Exata
      ctx.save();
      ctx.strokeStyle = '#3b82f6';
      ctx.lineWidth = 2.0;
      ctx.setLineDash([6, 6]);
      ctx.beginPath();
      const numPts = 100;
      for (let i = 0; i <= numPts; i++) {
        const px = xMin + (i / numPts) * (xMax - xMin);
        const py = methodData.exact(px);
        const sx = toScreenX(px);
        const sy = toScreenY(py);
        if (i === 0) ctx.moveTo(sx, sy);
        else ctx.lineTo(sx, sy);
      }
      ctx.stroke();
      ctx.setLineDash([]);

      ctx.fillStyle = '#3b82f6';
      ctx.font = 'bold 13px "JetBrains Mono", monospace';
      ctx.fillText('SOLUÇÃO EXATA y(x)', canvasWidth - margin - 170, margin + 20);

      // Trajetória Numérica Calculada até o passo atual
      ctx.strokeStyle = CRTTheme.green;
      ctx.lineWidth = 3.5;
      ctx.shadowColor = CRTTheme.green;
      ctx.shadowBlur = 10;
      ctx.beginPath();

      const initialStep = methodData.steps[0];
      ctx.moveTo(toScreenX(initialStep.x), toScreenY(initialStep.y));

      for (let i = 0; i <= stepIndex && i < methodData.steps.length; i++) {
        const st = methodData.steps[i];
        if (i < stepIndex) {
          ctx.lineTo(toScreenX(st.xNext), toScreenY(st.yNext));
        } else {
          // Interpolação suave do passo corrente
          const interX = st.x + (st.xNext - st.x) * animProgress;
          const interY = st.y + (st.yNext - st.y) * animProgress;
          ctx.lineTo(toScreenX(interX), toScreenY(interY));
        }
      }
      ctx.stroke();
      ctx.shadowBlur = 0;

      // Desenho dos nós discretos
      for (let i = 0; i <= stepIndex && i < methodData.steps.length; i++) {
        const st = methodData.steps[i];
        const sx = toScreenX(st.x);
        const sy = toScreenY(st.y);

        ctx.fillStyle = CRTTheme.green;
        ctx.beginPath();
        ctx.arc(sx, sy, 5, 0, Math.PI * 2);
        ctx.fill();

        if (i === stepIndex) {
          // Vetor Tangente / Inclinômetro de Heun ou Euler
          const curSX = toScreenX(st.x);
          const curSY = toScreenY(st.y);
          const nextSX = toScreenX(st.xNext);
          const nextSY = toScreenY(st.yNext);

          ctx.strokeStyle = CRTTheme.amber;
          ctx.lineWidth = 2.0;
          ctx.beginPath();
          ctx.moveTo(curSX, curSY);
          ctx.lineTo(nextSX, nextSY);
          ctx.stroke();

          // Nó de avanço
          ctx.fillStyle = CRTTheme.amber;
          ctx.shadowColor = CRTTheme.amber;
          ctx.shadowBlur = 12;
          ctx.beginPath();
          ctx.arc(nextSX, nextSY, 6, 0, Math.PI * 2);
          ctx.fill();
        }
      }

      // HUD de telemetria do nó ativo
      const curStep = methodData.steps[stepIndex];
      if (curStep) {
        ctx.font = 'bold 15px "JetBrains Mono", monospace';
        ctx.fillStyle = CRTTheme.green;
        ctx.fillText(`NÓ x_${curStep.i + 1} = ${curStep.xNext.toFixed(2)}  //  y_num = ${curStep.yNext.toFixed(6)}`, margin + 20, margin + 20);
        ctx.fillStyle = CRTTheme.amber;
        ctx.fillText(`y_exato = ${curStep.yExactNext.toFixed(6)}  //  Erro = ${curStep.errRel.toFixed(4)}%`, margin + 20, margin + 44);
      }

      ctx.restore();
      this.drawScanlines(ctx, canvasWidth, canvasHeight);
    }
  };

  // =========================================================================
  // 9. MOTOR INTERATIVO & MÁQUINA DE ESTADOS (NumericalEngine)
  // =========================================================================
  class NumericalEngine {
    constructor() {
      this.canvas = null;
      this.ctx = null;
      this.activeMethodKey = 'newton';
      this.activeMethodData = null;

      this.currentStep = 0;
      this.animProgress = 0.0;
      this.isPlaying = true;
      this.speed = 1.0;
      this.frame = 0;
      this.animationFrameId = null;

      this.subscribers = new Set();
      this.keyListenersAttached = false;
      this._boundOnKeyDown = this._onKeyDown.bind(this);
      this._boundLoop = this._loop.bind(this);
    }

    /**
     * Inicializa o motor com um elemento canvas HTML (ou mock em testes)
     */
    init(canvasElement, initialMethod = 'newton') {
      if (canvasElement) {
        this.canvas = canvasElement;
        this.ctx = canvasElement.getContext('2d');
      }
      this.loadMethod(initialMethod);
      this.attachKeyboard();
      this.start();
      return this;
    }

    /**
     * Carrega e resolve a família matemática do método selecionado
     */
    loadMethod(methodKey) {
      const entry = MethodsManifest[methodKey];
      if (!entry) {
        throw new Error(`Método numérico desconhecido: "${methodKey}"`);
      }

      this.activeMethodKey = methodKey;
      this.activeMethodData = entry.solver();
      this.currentStep = 0;
      this.animProgress = 0.0;
      this.frame = 0;

      this._notifySubscribers('methodLoaded');
      this.render();
      return this;
    }

    /**
     * Loop principal de animação a 60 FPS
     */
    _loop() {
      if (this.isPlaying && this.activeMethodData) {
        const stepRate = 0.015 * this.speed;
        this.animProgress += stepRate;
        this.frame = (this.frame + 1) % 60000;

        const totalSteps = this.activeMethodData.steps.length;
        if (this.animProgress >= 1.0) {
          this.animProgress = 0.0;
          this.currentStep = (this.currentStep + 1) % totalSteps;
          this._notifySubscribers('stepChanged');
        }
      }

      this.render();

      if (typeof requestAnimationFrame !== 'undefined') {
        this.animationFrameId = requestAnimationFrame(this._boundLoop);
      }
    }

    /**
     * Inicia o loop de animação
     */
    start() {
      if (this.animationFrameId) cancelAnimationFrame(this.animationFrameId);
      if (typeof requestAnimationFrame !== 'undefined') {
        this.animationFrameId = requestAnimationFrame(this._boundLoop);
      }
    }

    /**
     * Renderiza o quadro atual no canvas
     */
    render() {
      if (!this.ctx || !this.canvas || !this.activeMethodData) return;

      const w = this.canvas.width;
      const h = this.canvas.height;
      const manifest = MethodsManifest[this.activeMethodKey];

      switch (manifest.family) {
        case 'scalar_root':
          CRTRenderer.renderScalarRoot(this.ctx, this.activeMethodData, this.currentStep, this.animProgress, w, h);
          break;
        case 'linear_system':
          CRTRenderer.renderMatrixElimination(this.ctx, this.activeMethodData, this.currentStep, this.animProgress, w, h);
          break;
        case 'quadrature':
          CRTRenderer.renderQuadrature(this.ctx, this.activeMethodData, this.currentStep, this.animProgress, w, h);
          break;
        case 'ode':
          CRTRenderer.renderODEStepper(this.ctx, this.activeMethodData, this.currentStep, this.animProgress, w, h);
          break;
        default:
          CRTRenderer.drawGrid(this.ctx, w, h);
          break;
      }
    }

    // -----------------------------------------------------------------------
    // CONTROLES DE REPRODUÇÃO & NAVEGAÇÃO
    // -----------------------------------------------------------------------
    play() {
      this.isPlaying = true;
      this._notifySubscribers('stateChanged');
    }

    pause() {
      this.isPlaying = false;
      this._notifySubscribers('stateChanged');
    }

    togglePlay() {
      this.isPlaying = !this.isPlaying;
      this._notifySubscribers('stateChanged');
      return this.isPlaying;
    }

    stepForward() {
      if (!this.activeMethodData) return;
      this.pause();
      const totalSteps = this.activeMethodData.steps.length;
      this.currentStep = (this.currentStep + 1) % totalSteps;
      this.animProgress = 0.99;
      this.render();
      this._notifySubscribers('stepChanged');
    }

    stepBackward() {
      if (!this.activeMethodData) return;
      this.pause();
      const totalSteps = this.activeMethodData.steps.length;
      this.currentStep = (this.currentStep - 1 + totalSteps) % totalSteps;
      this.animProgress = 0.99;
      this.render();
      this._notifySubscribers('stepChanged');
    }

    reset() {
      this.currentStep = 0;
      this.animProgress = 0.0;
      this.frame = 0;
      this.render();
      this._notifySubscribers('stepChanged');
    }

    setSpeed(speedMultiplier) {
      this.speed = MathUtils.clamp(speedMultiplier, 0.1, 5.0);
      this._notifySubscribers('speedChanged');
    }

    seek(progressPercent) {
      if (!this.activeMethodData) return;
      this.pause();
      const totalSteps = this.activeMethodData.steps.length;
      const stepFloat = (MathUtils.clamp(progressPercent, 0, 100) / 100) * totalSteps;
      this.currentStep = Math.min(totalSteps - 1, Math.floor(stepFloat));
      this.animProgress = stepFloat - this.currentStep;
      this.render();
      this._notifySubscribers('stepChanged');
    }

    // -----------------------------------------------------------------------
    // ATALHOS DE TECLADO (HOTKEYS DE TERMINAL)
    // -----------------------------------------------------------------------
    _onKeyDown(e) {
      // Ignora atalhos se o foco estiver em um input de texto
      if (e.target && (e.target.tagName === 'INPUT' || e.target.tagName === 'TEXTAREA')) {
        return;
      }

      if (e.code === 'Space') {
        e.preventDefault();
        this.togglePlay();
      } else if (e.key === ']' || e.key === 'ArrowRight') {
        e.preventDefault();
        this.stepForward();
      } else if (e.key === '[' || e.key === 'ArrowLeft') {
        e.preventDefault();
        this.stepBackward();
      } else if (e.key === 'r' || e.key === 'R') {
        e.preventDefault();
        this.reset();
      } else if (e.key === '1') {
        this.loadMethod('newton');
      } else if (e.key === '2') {
        this.loadMethod('bissecao');
      } else if (e.key === '3') {
        this.loadMethod('secante');
      } else if (e.key === '4') {
        this.loadMethod('falsa_posicao');
      } else if (e.key === '5') {
        this.loadMethod('ponto_fixo');
      } else if (e.key === '6') {
        this.loadMethod('gauss_pivoteamento');
      } else if (e.key === '7') {
        this.loadMethod('gauss_jordan');
      } else if (e.key === '8') {
        this.loadMethod('trapezio_composto');
      } else if (e.key === '9') {
        this.loadMethod('rk2_heun');
      }
    }

    attachKeyboard(target = null) {
      const listenerTarget = target || (typeof window !== 'undefined' ? window : null);
      if (listenerTarget && !this.keyListenersAttached) {
        listenerTarget.addEventListener('keydown', this._boundOnKeyDown);
        this.keyListenersAttached = true;
      }
    }

    detachKeyboard(target = null) {
      const listenerTarget = target || (typeof window !== 'undefined' ? window : null);
      if (listenerTarget && this.keyListenersAttached) {
        listenerTarget.removeEventListener('keydown', this._boundOnKeyDown);
        this.keyListenersAttached = false;
      }
    }

    // -----------------------------------------------------------------------
    // SISTEMA DE SUBSCRIÇÃO / EVENTOS PARA O FRONTEND
    // -----------------------------------------------------------------------
    subscribe(callback) {
      this.subscribers.add(callback);
      // Notifica imediatamente o estado atual
      this._emitTo(callback, 'initial');
      return () => this.subscribers.delete(callback);
    }

    _emitTo(cb, eventType) {
      if (!this.activeMethodData) return;
      cb({
        event: eventType,
        methodKey: this.activeMethodKey,
        manifest: MethodsManifest[this.activeMethodKey],
        methodData: this.activeMethodData,
        currentStep: this.currentStep,
        totalSteps: this.activeMethodData.steps.length,
        stepData: this.activeMethodData.steps[this.currentStep],
        isPlaying: this.isPlaying,
        animProgress: this.animProgress,
        speed: this.speed,
        frame: this.frame
      });
    }

    _notifySubscribers(eventType) {
      for (const cb of this.subscribers) {
        try {
          this._emitTo(cb, eventType);
        } catch (err) {
          console.error('[NumericalEngine] Erro no callback do assinante:', err);
        }
      }
    }

    destroy() {
      if (this.animationFrameId) cancelAnimationFrame(this.animationFrameId);
      this.detachKeyboard();
      this.subscribers.clear();
      this.canvas = null;
      this.ctx = null;
    }
  }

  // =========================================================================
  // 10. EXPORTAÇÃO UNIVERSAL (BROWSER / NODE.JS)
  // =========================================================================
  const EngineAPI = {
    CRTTheme,
    MathUtils,
    ScalarRootSolvers,
    LinearSystemSolvers,
    QuadratureSolvers,
    ODESolvers,
    MethodsManifest,
    CRTRenderer,
    NumericalEngine,
    createEngine: (canvas, initialMethod) => new NumericalEngine().init(canvas, initialMethod)
  };

  // Exportação Browser
  if (typeof window !== 'undefined') {
    window.NumericalEngineAPI = EngineAPI;
    window.NumericalEngine = NumericalEngine;
    window.CRTTheme = CRTTheme;
    window.CRTRenderer = CRTRenderer;
    window.MethodsManifest = MethodsManifest;
  }

  // Exportação GlobalThis
  if (typeof globalThis !== 'undefined') {
    globalThis.NumericalEngineAPI = EngineAPI;
  }

  // Exportação CommonJS (Node.js)
  if (typeof module !== 'undefined' && module.exports) {
    module.exports = EngineAPI;
  }

})(typeof globalThis !== 'undefined' ? globalThis : typeof window !== 'undefined' ? window : this);
