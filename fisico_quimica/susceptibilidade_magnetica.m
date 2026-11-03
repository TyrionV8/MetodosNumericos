% ==============================================================================
% FÍSICO-QUÍMICA / MAGNETOQUÍMICA: MOMENTO MAGNÉTICO E SUSCEPTIBILIDADE MOLAR
% ==============================================================================
%
% 1. O QUE CALCULA:
%    Este script calcula o Momento Magnético Efetivo de Spin Puro (spin-only) e
%    a Susceptibilidade Magnética Molar (chi_M) para íons de metais de transição
%    e complexos de coordenação paramagnéticos.
%
%    Formulação Matemática:
%    a) Momento Magnético Efetivo de Spin Puro (em Magnetons de Bohr, mu_B):
%       mu_so = sqrt( n * (n + 2) )
%       Onde 'n' representa a quantidade discreta de elétrons desemparelhados
%       (obtida pela regra de Hund e configuração de campo cristalino de orbitais d).
%       Em termos do número quântico de spin total S = n / 2:
%       mu_so = 2 * sqrt( S * (S + 1) )
%
%    b) Susceptibilidade Magnética Molar pela Lei de Curie (chi_M em cm^3/mol):
%       chi_M = ( (mu_so / constante_proporcionalidade)^2 ) / T
%       Onde a constante_proporcionalidade clássica aproximada é 2.84,
%       derivada de sqrt( (3 * k_B) / (N_A * beta^2) ) ≈ 2.828 - 2.84,
%       e T é a temperatura termodinâmica absoluta em Kelvin (K).
%
% 2. O QUE ENTREGA:
%    - Valor do número quântico de spin total S e multiplicidade de spin (2S + 1).
%    - Momento magnético teórico de spin puro em Magnetons de Bohr (mu_B).
%    - Constante de Curie molar calculada (C_molar).
%    - Susceptibilidade magnética molar resultante (chi_M) em unidades CGS (cm^3/mol ou emu/mol).
%
% 3. O QUE PRECISA TER (PRÉ-REQUISITOS E PARÂMETROS):
%    - Quantidade inteira de elétrons desemparelhados n (geralmente entre 1 e 5 para subnível d).
%    - Temperatura absoluta de medição em Kelvin (padrão ambiente: 298.15 K).
%
% 4. EXEMPLO DE USO E RESULTADO ESPERADO:
%    - Para n = 4 (exemplo: íon Fe2+ d6 em campo fraco de spin alto ou Cr2+ d4):
%      Temperatura T = 298.15 K.
%    - Resultado esperado:
%      mu_so = sqrt(4 * 6) = sqrt(24) ≈ 4.898979 mu_B
%      chi_M ≈ 0.009980 cm^3/mol (9.98e-3 emu/mol)
% ==============================================================================

% Limpeza do ambiente de execução
clc;        % Limpa a tela do console de comandos
clear all;  % Remove todas as variáveis ativas do espaço de trabalho (workspace)
close all;  % Fecha todas as janelas de figuras abertas anteriormente

% ==============================================================================
% BLOCO DE PARÂMETROS DE ENTRADA
% ==============================================================================

% Quantidade discreta de elétrons desemparelhados na camada de valência
numero_eletrons_desemparelhados = 4; % Valor inteiro adimensional (n = 1, 2, 3, 4, 5, ...)

% Temperatura absoluta da amostra em escala termodinâmica de Kelvin
temperatura_absoluta_kelvin = 298.15; % Padrão de referência para temperatura ambiente (25 °C)

% Constante fenomenológica empírica da magnetoquímica clássica de Curie
% Relação teórica: mu_eff = 2.84 * sqrt(chi_M * T)
constante_proporcionalidade_curie = 2.84; % Fator dimensional empírico de calibração

% ==============================================================================
% VALIDAÇÃO DOS PARÂMETROS FÍSICOS DE ENTRADA
% ==============================================================================

if numero_eletrons_desemparelhados < 0 || floor(numero_eletrons_desemparelhados) ~= numero_eletrons_desemparelhados
    error('O número de elétrons desemparelhados deve ser um número inteiro maior ou igual a zero.');
end
if temperatura_absoluta_kelvin <= 0.0
    error('A temperatura termodinâmica absoluta deve ser estritamente positiva em Kelvin.');
end

% ==============================================================================
% PROCESSAMENTO MATEMÁTICO E RELAÇÕES DE SPIN
% ==============================================================================

% Número quântico de spin total do sistema eletrônico (cada elétron desemparelhado contribui com s = 1/2)
numero_quantico_spin_total = numero_eletrons_desemparelhados / 2.0;

% Multiplicidade espectroscópica de spin (degenerescência 2S + 1)
multiplicidade_espectroscopica_spin = 2.0 * numero_quantico_spin_total + 1.0;

% Cálculo analítico do momento magnético efetivo de spin puro em Magnetons de Bohr (mu_B)
momento_magnetico_spin_puro_bohr_magnetons = sqrt(numero_eletrons_desemparelhados * (numero_eletrons_desemparelhados + 2.0));

% Constante molar de Curie do composto C = chi_M * T = (mu_eff / 2.84)^2
constante_curie_molar = (momento_magnetico_spin_puro_bohr_magnetons / constante_proporcionalidade_curie) ^ 2;

% Susceptibilidade magnética molar (chi_M) em unidades CGS eletromagnéticas (cm^3 / mol ou emu / mol)
susceptibilidade_magnetica_molar_cgs = constante_curie_molar / temperatura_absoluta_kelvin;

% Conversão para o Sistema Internacional (SI): chi_M_SI = 4 * pi * 1e-6 * chi_M_CGS (m^3 / mol)
fator_conversao_cgs_para_si = 4.0 * pi * 1.0e-6;
susceptibilidade_magnetica_molar_si = susceptibilidade_magnetica_molar_cgs * fator_conversao_cgs_para_si;

% ==============================================================================
% APRESENTAÇÃO DOS RESULTADOS NO CONSOLE
% ==============================================================================

fprintf('\n==============================================================================\n');
fprintf('       MAGNETOQUÍMICA: MOMENTO MAGNÉTICO E SUSCEPTIBILIDADE DE CURIE          \n');
fprintf('==============================================================================\n');
fprintf('  CONDIÇÕES EXPERIMENTAIS E DE CONFIGURAÇÃO ELETRÔNICA:\n');
fprintf('    - Número de elétrons desemparelhados (n):    %d\n', numero_eletrons_desemparelhados);
fprintf('    - Número quântico de spin total (S):        %6.1f\n', numero_quantico_spin_total);
fprintf('    - Multiplicidade de spin (2S + 1):          %6.1f\n', multiplicidade_espectroscopica_spin);
fprintf('    - Temperatura termodinâmica de ensaio (T):   %8.2f K (%5.2f °C)\n', ...
        temperatura_absoluta_kelvin, temperatura_absoluta_kelvin - 273.15);
fprintf('------------------------------------------------------------------------------\n');
fprintf('  GRANDEZAS MAGNÉTICAS CALCULADAS:\n');
fprintf('    - Momento Magnético de Spin Puro (mu_so):   %10.6f mu_B (Magnetons de Bohr)\n', ...
        momento_magnetico_spin_puro_bohr_magnetons);
fprintf('    - Constante Molar de Curie (C):             %10.6f cm^3*K/mol\n', constante_curie_molar);
fprintf('    - Susceptibilidade Molar (chi_M - CGS):      %10.6f cm^3/mol (emu/mol)\n', ...
        susceptibilidade_magnetica_molar_cgs);
fprintf('    - Susceptibilidade Molar (chi_M - SI):       %12.5e m^3/mol\n', ...
        susceptibilidade_magnetica_molar_si);
fprintf('==============================================================================\n\n');
