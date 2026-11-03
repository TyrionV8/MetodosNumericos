% ==============================================================================
% FÍSICO-QUÍMICA / MAGNETOQUÍMICA: DETERMINAÇÃO DO NÚMERO DE ELÉTRONS DESEMPARELHADOS
% ==============================================================================
%
% 1. O QUE CALCULA:
%    Este script resolve o problema inverso fundamental da magnetoquímica:
%    a partir de dados experimentais da Susceptibilidade Magnética Molar (chi_M)
%    e da Temperatura Absoluta (T), determina o número inteiro de elétrons
%    desemparelhados (n) presentes no centro iônico metálico.
%
%    Formulação Matemática:
%    A relação empírica de Curie estabelece o momento magnético efetivo como:
%       mu_eff = constante_proporcionalidade * sqrt( chi_M * T )
%       com constante_proporcionalidade ≈ 2.84, de modo que (2.84)^2 ≈ 8.0656.
%
%    Pelo modelo quântico de spin puro (spin-only):
%       mu_eff = sqrt( n * (n + 2) )
%       ==> mu_eff^2 = n^2 + 2*n
%
%    Igualando as duas relações:
%       n^2 + 2*n = 8.0656 * chi_M * T
%       n^2 + 2*n - C = 0,   onde C = 8.0656 * chi_M * T
%
%    Esta equação algébrica de segundo grau é resolvida numericamente encontrando
%    as raízes do polinômio quadrático P(n) = [1, 2, -C] através da função 'roots'.
%    Descarta-se a raiz espúria negativa e seleciona-se a raiz positiva, que é
%    arredondada ao número inteiro mais próximo, compatível com a natureza quantizada
%    do spin eletrônico.
%
% 2. O QUE ENTREGA:
%    - Raízes numéricas completas do polinômio quadrático.
%    - Identificação e filtragem da raiz fisicamente significativa (n > 0).
%    - Quantidade inteira discreta de elétrons desemparelhados n.
%    - Momento magnético efetivo correspondente em Magnetons de Bohr (mu_B).
%    - Configuração de spin resultante (spin total S e multiplicidade 2S + 1).
%
% 3. O QUE PRECISA TER (PRÉ-REQUISITOS E PARÂMETROS):
%    - Susceptibilidade magnética molar experimental (chi_M em cm^3/mol ou emu/mol).
%    - Temperatura absoluta em Kelvin na qual a medição foi efetuada.
%
% 4. EXEMPLO DE USO E RESULTADO ESPERADO:
%    - Entrada: chi_M = 3.33e-2 cm^3/mol, T = 298.15 K.
%    - Coeficiente independente C ≈ 80.078.
%    - Raízes calculadas: n_1 ≈ +8.004 e n_2 ≈ -10.004.
%    - Resultado esperado: n = 8 elétrons desemparelhados.
% ==============================================================================

% Limpeza do ambiente de execução
clc;        % Limpa a tela do console de comandos
clear all;  % Remove todas as variáveis ativas do espaço de trabalho (workspace)
close all;  % Fecha todas as janelas de figuras abertas anteriormente

% ==============================================================================
% BLOCO DE PARÂMETROS DE ENTRADA EXPERIMENTAIS
% ==============================================================================

% Temperatura termodinâmica absoluta da medição magnetoquímica em Kelvin
temperatura_ensaio_kelvin = 298.15; % Temperatura ambiente padrão de 25 °C

% Susceptibilidade magnética molar experimental medida em unidades CGS (cm^3 / mol ou emu / mol)
susceptibilidade_magnetica_molar_experimental_cgs = 3.33e-2;

% Fator multiplicador fenomenológico empírico clássico (2.84^2)
fator_calibracao_curie = 8.0656; % Derivado de (2.84)^2

% ==============================================================================
% VALIDAÇÃO DOS PARÂMETROS EXPERIMENTAIS DE ENTRADA
% ==============================================================================

if temperatura_ensaio_kelvin <= 0.0
    error('A temperatura absoluta de ensaio deve ser estritamente positiva em Kelvin.');
end
if susceptibilidade_magnetica_molar_experimental_cgs <= 0.0
    error('A susceptibilidade magnética molar experimental deve ser positiva para substâncias paramagnéticas.');
end

% ==============================================================================
% PROCESSAMENTO MATEMÁTICO: MONTAGEM DO POLINÔMIO E RESOLUÇÃO POR ROOTS
% ==============================================================================

% Cálculo do termo independente C da equação quadrática n^2 + 2*n - C = 0
termo_independente_constante_c = susceptibilidade_magnetica_molar_experimental_cgs * ...
                                 temperatura_ensaio_kelvin * ...
                                 fator_calibracao_curie;

% Vetor de coeficientes do polinômio de segundo grau na forma: [a, b, c]
% correspondente a: a * n^2 + b * n + c = 0
coeficiente_quadratico_a = 1.0;
coeficiente_linear_b     = 2.0;
coeficiente_constante_c  = -termo_independente_constante_c;

vetor_coeficientes_polinomio = [coeficiente_quadratico_a, coeficiente_linear_b, coeficiente_constante_c];

% Determinação das raízes do polinômio através do método de autovalores de matriz companheira (roots)
vetor_raizes_polinomio = roots(vetor_coeficientes_polinomio);

% Seleção da raiz positiva fisicamente plausível (número de partículas deve ser não negativo)
raizes_reais_positivas = vetor_raizes_polinomio(vetor_raizes_polinomio > 0);

if isempty(raizes_reais_positivas)
    error('Nenhuma raiz positiva encontrada. Verifique os dados de entrada de susceptibilidade.');
end

% Extração analítica explícita das raízes calculadas
raiz_positiva_fisica     = max(raizes_reais_positivas);
raiz_negativa_descartada = min(vetor_raizes_polinomio);

% Extração do valor contínuo da raiz física
estimativa_continua_eletrons = raiz_positiva_fisica;

% Arredondamento ao número inteiro mais próximo devido à quantização fundamental do elétron
numero_inteiro_eletrons_desemparelhados = round(estimativa_continua_eletrons);

% Resíduo numérico da quantização (diferença entre a medição experimental e o valor inteiro exato)
residuo_quantizacao_experimental = abs(estimativa_continua_eletrons - numero_inteiro_eletrons_desemparelhados);

% Cálculo do momento magnético efetivo teórico associado a este número inteiro de elétrons
momento_magnetico_spin_puro_teorico = sqrt(numero_inteiro_eletrons_desemparelhados * (numero_inteiro_eletrons_desemparelhados + 2.0));

% Número quântico de spin total S = n / 2
numero_quantico_spin_total = numero_inteiro_eletrons_desemparelhados / 2.0;

% Multiplicidade espectroscópica (2S + 1)
multiplicidade_spin = 2.0 * numero_quantico_spin_total + 1.0;

% ==============================================================================
% APRESENTAÇÃO DOS RESULTADOS NO CONSOLE
% ==============================================================================

fprintf('\n==============================================================================\n');
fprintf('     MAGNETOQUÍMICA: DETERMINAÇÃO DO NÚMERO DE ELÉTRONS DESEMPARELHADOS (n)  \n');
fprintf('==============================================================================\n');
fprintf('  DADOS EXPERIMENTAIS FORNECIDOS:\n');
fprintf('    - Temperatura Absoluta (T):                  %8.2f K\n', temperatura_ensaio_kelvin);
fprintf('    - Susceptibilidade Molar (chi_M):            %11.4e cm^3/mol\n', susceptibilidade_magnetica_molar_experimental_cgs);
fprintf('    - Termo Independente do Polinômio (C):       %10.4f\n', termo_independente_constante_c);
fprintf('------------------------------------------------------------------------------\n');
fprintf('  RESOLUÇÃO ALGÉBRICA VIA ROOTS:\n');
fprintf('    - Equação Algébrica Formada:                 1.0*n^2 + 2.0*n - (%.4f) = 0\n', termo_independente_constante_c);
fprintf('    - Raiz 1 (positiva física):                  %+10.5f\n', raiz_positiva_fisica);
fprintf('    - Raiz 2 (negativa descartada):              %+10.5f\n', raiz_negativa_descartada);
fprintf('------------------------------------------------------------------------------\n');
fprintf('  INTERPRETAÇÃO FÍSICA E ESTADO QUÂNTICO:\n');
fprintf('    - Estimativa contínua não arredondada:       %10.5f elétrons\n', estimativa_continua_eletrons);
fprintf('    - Número Inteiro de Elétrons (n):            %d elétrons desemparelhados\n', numero_inteiro_eletrons_desemparelhados);
fprintf('    - Resíduo de Quantização Experimental:       %10.5f\n', residuo_quantizacao_experimental);
fprintf('    - Momento Magnético de Spin Puro (mu_so):    %10.5f mu_B\n', momento_magnetico_spin_puro_teorico);
fprintf('    - Spin Total do Sistema (S = n/2):           %6.1f\n', numero_quantico_spin_total);
fprintf('    - Multiplicidade de Degenerescência (2S+1):  %6.1f\n', multiplicidade_spin);
fprintf('==============================================================================\n\n');
