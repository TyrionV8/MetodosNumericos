% ==============================================================================
% ÁLGEBRA LINEAR E MÉTODOS NUMÉRICOS: OPERAÇÕES ARITMÉTICAS E MATRICIAIS
% ==============================================================================
%
% 1. O QUE CALCULA:
%    Este script unifica e formaliza os exercícios clássicos de introdução à
%    computação científica numérica (Lista 1: Questões 1A, 1B e Questão 2),
%    cobrindo dois pilares fundamentais do cálculo computacional:
%
%    PARTE 1: Expressões Aritméticas Encadeadas, Potenciação e Radiciação
%    - Resolução da Questão 1A:
%      F_A = A + [ (A * B) / C ] * [ ( A * B * (A + D)^2 ) / sqrt(A * B) ]
%    - Resolução da Questão 1B:
%      F_B = D * E^(D/2) + [ (A*D + C*D) / (25/A + 35/B) ] / (A + B + C + D)
%      Avalia precedência de operadores algébricos, potências fracionárias e
%      estabilidade numérica em operações de divisão com denominadores compostos.
%
%    PARTE 2: Álgebra Linear Matricial (Questão 2 Expandida)
%    - Geração de Matriz Quadrada Inteira de ordem n = 6 no intervalo [1, 9].
%    - Transposição Matricial Hermitiana e Clássica: M_transposta = M^T.
%    - Análise de Regularidade e Invertibilidade: Cálculo do Determinante det(M)
%      e Número de Condicionamento cond(M).
%    - Inversão Matricial Numérica: M_inversa = inv(M).
%    - Validação Formal do Teorema Fundamental da Inversa:
%      Verificação do produto interno M * M^(-1) ≈ Identidade I_n,
%      com quantificação do resíduo de erro matricial máximo (norma infinito).
%    - Decomposição Espectral Adicional: Cálculo do traço e espectro de autovalores.
%
% 2. O QUE ENTREGA:
%    - Apresentação passo a passo no console das frações parciais e resultado final da Parte 1.
%    - Exibição visual formatada da matriz original, sua transposta e sua inversa na Parte 2.
%    - Tabela de verificação da identidade com cálculo do erro residual numérico de máquina.
%
% 3. O QUE PRECISA TER (PRÉ-REQUISITOS E PARÂMETROS):
%    - Parâmetros escalares inteiros da lista: A, B, C, D, E.
%    - Dimensão da matriz quadrada (n = 6) e limite superior de inteiros aleatórios (valor máximo = 9).
%
% 4. EXEMPLO DE USO E RESULTADO ESPERADO:
%    - Executar diretamente no terminal:
%      >> operacoes_matrizes
%    - Resultado esperado:
%      F_A ≈ 3415520.4578
%      F_B ≈ 26873860.7620
%      Erro residual ||M * inv(M) - I|| < 1e-14 (precisão quase exata de máquina IEEE 754).
% ==============================================================================

% Limpeza do ambiente de execução
clc;        % Limpa a tela do console de comandos
clear all;  % Remove todas as variáveis ativas do espaço de trabalho (workspace)
close all;  % Fecha todas as janelas de figuras abertas anteriormente

% ==============================================================================
% PARTE 1: AVALIAÇÃO DE EXPRESSÕES ARITMÉTICAS ESCALARES COMPLEXAS
% ==============================================================================

% Definição dos valores escalares inteiros de entrada da lista
parametro_escalar_a = 22.0;
parametro_escalar_b = 183.0;
parametro_escalar_c = 108.0;
parametro_escalar_d = 16.0;
parametro_escalar_e = 6.0;

% ------------------------------------------------------------------------------
% Questão 1A: Combinação de produtos, potências ao quadrado e raiz quadrada
% ------------------------------------------------------------------------------

% Primeira fração parcial: Div1 = (A * B) / C
fracao_parcial_divisao_1 = (parametro_escalar_a * parametro_escalar_b) / parametro_escalar_c;

% Numerador da segunda fração: Cima = A * B * (A + D)^2
fator_soma_ad_ao_quadrado = (parametro_escalar_a + parametro_escalar_d) ^ 2.0;
numerador_fracao_parcial_2 = parametro_escalar_a * parametro_escalar_b * fator_soma_ad_ao_quadrado;

% Denominador da segunda fração: Baixo = sqrt(A * B)
denominador_fracao_parcial_2 = sqrt(parametro_escalar_a * parametro_escalar_b);

% Segunda fração parcial: Div2 = Numerador / Denominador
fracao_parcial_divisao_2 = numerador_fracao_parcial_2 / denominador_fracao_parcial_2;

% Resultado final da Questão 1A: Funcao_1A = A + (Div1 * Div2)
resultado_final_questao_1a = parametro_escalar_a + (fracao_parcial_divisao_1 * fracao_parcial_divisao_2);

% ------------------------------------------------------------------------------
% Questão 1B: Expoentes racionais fracionários e somatório de taxas
% ------------------------------------------------------------------------------

% Parcela exponencial preliminar: S = D * (E ^ (D / 2))
expoente_fracionario_d_sobre_dois = parametro_escalar_d / 2.0; % 16 / 2 = 8
termo_potencia_exponencial_s = parametro_escalar_d * (parametro_escalar_e ^ expoente_fracionario_d_sobre_dois);

% Numerador composto da fração: Um = (A * D) + (C * D)
termo_soma_produtos_numerador = (parametro_escalar_a * parametro_escalar_d) + (parametro_escalar_c * parametro_escalar_d);

% Denominador interno de frações: Outro = (25 / A) + (35 / B)
termo_soma_razoes_fracionarias = (25.0 / parametro_escalar_a) + (35.0 / parametro_escalar_b);

% Razão superior: Cima = Um / Outro
numerador_composto_questao_1b = termo_soma_produtos_numerador / termo_soma_razoes_fracionarias;

% Denominador geral: Baixo = A + B + C + D
denominador_composto_questao_1b = parametro_escalar_a + parametro_escalar_b + parametro_escalar_c + parametro_escalar_d;

% Razão racional secundária: Termo_Racional = Cima / Baixo
termo_racional_secundario_1b = numerador_composto_questao_1b / denominador_composto_questao_1b;

% Resultado final da Questão 1B: Funcao_1B = S + Termo_Racional
resultado_final_questao_1b = termo_potencia_exponencial_s + termo_racional_secundario_1b;

% ==============================================================================
% PARTE 2: ÁLGEBRA LINEAR MATRICIAL (GERAÇÃO, TRANSPOSIÇÃO E INVERSÃO)
% ==============================================================================

% Dimensão da matriz quadrada (número de linhas e colunas)
ordem_matriz_quadrada = 6;

% Valor inteiro máximo admitido nos elementos aleatórios da matriz
limite_superior_inteiro_aleatorio = 9;

% Fixação de semente de gerador pseudoaleatório para garantir reprodutibilidade exata
% Compatível com versões modernas de Octave e MATLAB
if exist('rng', 'file') || exist('rng', 'builtin')
    rng(42);
else
    rand('seed', 42);
end

% Geração da matriz quadrada com elementos inteiros aleatórios no intervalo [1, 9]
matriz_quadrada_original = randi(limite_superior_inteiro_aleatorio, ordem_matriz_quadrada, ordem_matriz_quadrada);

% Operação 1: Transposição Matricial (M^T)
matriz_transposta = matriz_quadrada_original';

% Operação 2: Análise de Invertibilidade (Cálculo do Determinante)
determinante_matriz_original = det(matriz_quadrada_original);

% Operação 3: Inversão Numérica via Eliminação de Gauss-Jordan (M^(-1))
if abs(determinante_matriz_original) < 1e-10
    error('A matriz gerada é singular ou quase-singular (det ≈ 0) e não admite inversa estável.');
else
    matriz_inversa_calculada = inv(matriz_quadrada_original);
end

% Operação 4: Validação do Teorema da Inversa (M * M^(-1) = Identidade)
matriz_identidade_referencia = eye(ordem_matriz_quadrada);
matriz_produto_validacao     = matriz_quadrada_original * matriz_inversa_calculada;

% Cálculo da matriz de resíduo e do erro absoluto máximo (norma do erro)
matriz_diferenca_residuo = abs(matriz_produto_validacao - matriz_identidade_referencia);
erro_residual_maximo_inversao = max(matriz_diferenca_residuo(:));

% Avaliação formal e teste de resíduo com tolerância de precisão numérica
tolerancia_admissivel_maquina = 1.0e-10;
if erro_residual_maximo_inversao < tolerancia_admissivel_maquina
    status_teste_residuo = 'CONFORME (Inversão válida com precisão de máquina)';
else
    status_teste_residuo = 'NÃO CONFORME (Erro residual superior à tolerância)';
end

% Cálculo do número de condicionamento na norma 2
numero_condicionamento_matriz = cond(matriz_quadrada_original);

% Traço da matriz (somatório dos elementos da diagonal principal)
traco_matriz_original = trace(matriz_quadrada_original);

% Espectro de autovalores da matriz
vetor_autovalores_matriz = eig(matriz_quadrada_original);

% ==============================================================================
% APRESENTAÇÃO DOS RESULTADOS NO CONSOLE
% ==============================================================================

fprintf('\n==============================================================================\n');
fprintf('     ÁLGEBRA LINEAR E CÁLCULO NUMÉRICO: OPERAÇÕES ARITMÉTICAS E MATRICIAIS   \n');
fprintf('==============================================================================\n');
fprintf('  PARTE 1: RESULTADOS ARITMÉTICOS ESCALARES (LISTA 1)\n');
fprintf('    - Parâmetros: A=%.0f, B=%.0f, C=%.0f, D=%.0f, E=%.0f\n', ...
        parametro_escalar_a, parametro_escalar_b, parametro_escalar_c, ...
        parametro_escalar_d, parametro_escalar_e);
fprintf('    --------------------------------------------------------------------------\n');
fprintf('    - Questão 1A:\n');
fprintf('        * Fração Parcial 1 [(A*B)/C]:                    %14.4f\n', fracao_parcial_divisao_1);
fprintf('        * Fração Parcial 2 [(A*B*(A+D)^2)/sqrt(A*B)]:   %14.4f\n', fracao_parcial_divisao_2);
fprintf('        * RESULTADO FINAL F_A:                           %14.4f\n', resultado_final_questao_1a);
fprintf('    --------------------------------------------------------------------------\n');
fprintf('    - Questão 1B:\n');
fprintf('        * Parcela Exponencial S [D * E^(D/2)]:           %14.4f\n', termo_potencia_exponencial_s);
fprintf('        * Numerador Composto Cima:                       %14.4f\n', numerador_composto_questao_1b);
fprintf('        * Denominador Composto Baixo:                    %14.4f\n', denominador_composto_questao_1b);
fprintf('        * RESULTADO FINAL F_B:                           %14.4f\n', resultado_final_questao_1b);
fprintf('==============================================================================\n');
fprintf('  PARTE 2: OPERAÇÕES COM MATRIZES QUADRADAS DE ORDEM 6x6\n');
fprintf('==============================================================================\n');
fprintf('  Matriz Original M (6x6):\n');
disp(matriz_quadrada_original);
fprintf('  Matriz Transposta M^T (6x6):\n');
disp(matriz_transposta);
fprintf('  Matriz Inversa M^(-1) (6x6):\n');
disp(matriz_inversa_calculada);
fprintf('------------------------------------------------------------------------------\n');
fprintf('  PROPRIEDADES ALGÉBRICAS E DE CONDICIONAMENTO:\n');
fprintf('    - Determinante det(M):                       %14.4f\n', determinante_matriz_original);
fprintf('    - Traço trace(M):                            %14.4f\n', traco_matriz_original);
fprintf('    - Número de Condicionamento cond(M):         %14.4f\n', numero_condicionamento_matriz);
fprintf('    - Erro Residual Máximo ||M * M^(-1) - I||:   %14.4e\n', erro_residual_maximo_inversao);
fprintf('    - Teste de Resíduo da Inversa:               %s\n', status_teste_residuo);
fprintf('==============================================================================\n\n');
