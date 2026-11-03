% ==============================================================================
% FÍSICO-QUÍMICA INORGÂNICA: EQUAÇÃO DE BORN-MAYER PARA ENERGIA RETICULAR
% ==============================================================================
%
% 1. O QUE CALCULA:
%    Este script calcula a Energia Reticular (U) de uma rede cristalina iônica
%    utilizando o modelo termodinâmico-eletrostático de Born-Mayer.
%
%    Formulação Matemática:
%       U = - ( (N_A * M * |z_+ * z_-| * e^2) / (4 * pi * epsilon_0 * r_0) ) * ( 1 - rho / r_0 )
%
%    Onde:
%       - N_A: Constante de Avogadro (partículas por mol)
%       - M: Constante geométrica de Madelung (adimensional, função da rede cristalina)
%       - z_+, z_-: Cargas formais (valências) dos íons cátion e ânion
%       - e: Carga elementar fundamental do elétron (Coulombs)
%       - epsilon_0: Permissividade elétrica absoluta do vácuo (Farads por metro)
%       - r_0: Distância interiônica de equilíbrio no retículo cristalino (metros)
%       - rho: Parâmetro fenomenológico de repulsão de Born-Mayer (metros)
%
%    A equação divide a energia em duas componentes físicas complementares:
%    1) Termo Eletrostático Atrativo de Coulomb (estabilizador de longo alcance).
%    2) Fator de Correção de Repulsão Quântica de Pauli (desestabilizador de curto alcance).
%
% 2. O QUE ENTREGA:
%    - Exibição estruturada no console com detalhamento das componentes físicas:
%      * Parcela puramente eletrostática atrativa (J/mol e kJ/mol)
%      * Fator adimensional de atenuação repulsiva
%      * Energia reticular líquida final U em Joules/mol e quilojoules/mol (kJ/mol)
%      * Comparação com a convenção termodinâmica padrão (U < 0 para formação espontânea)
%
% 3. O QUE PRECISA TER (PRÉ-REQUISITOS E PARÂMETROS):
%    - Constantes físicas fundamentais (CODATA mais recente).
%    - Raios iônicos de Shannon ou Pauling do cátion e do ânion em picômetros.
%    - Constante de Madelung correspondente à estrutura cristalina (ex.: 1.74756 para NaCl).
%    - Valências formais dos íons constituintes.
%
% 4. EXEMPLO DE USO E RESULTADO ESPERADO:
%    - Entrada padrão configurada para o Cloreto de Sódio (NaCl, Halita):
%      Cátion Na+ (102 pm, z = +1), Ânion Cl- (181 pm, z = -1), M = 1.74756, rho = 34.5 pm.
%    - Resultado esperado:
%      Distância interiônica r_0 = 283.00 pm (2.83e-10 m)
%      Energia Reticular U ≈ -753.35 kJ/mol
% ==============================================================================

% Limpeza do ambiente de execução
clc;        % Limpa a tela do console de comandos
clear all;  % Remove todas as variáveis ativas do espaço de trabalho (workspace)
close all;  % Fecha todas as janelas de figuras abertas anteriormente

% ==============================================================================
% BLOCO DE CONSTANTES FÍSICAS FUNDAMENTAIS (SISTEMA INTERNACIONAL - SI)
% ==============================================================================

% Constante de Avogadro (número de entidades elementares em 1 mol de substância)
constante_avogadro_particulas_por_mol = 6.02214076e23; % Unidade: mol^-1

% Carga elementar fundamental transportada por um único próton ou elétron
carga_elementar_eletron_coulombs = 1.602176634e-19; % Unidade: Coulombs (C)

% Permissividade elétrica absoluta do vácuo (constante dielétrica do espaço livre)
permissividade_eletrica_vacuo_farads_por_metro = 8.8541878128e-12; % Unidade: F/m (ou C^2 / (N * m^2))

% ==============================================================================
% BLOCO DE PARÂMETROS CRISTALOGRÁFICOS E IÔNICOS DE ENTRADA
% ==============================================================================

% Valência e número de oxidação formal do cátion metálico (ex.: Na+ possui carga +1)
valencia_formal_cation_adimensional = +1;

% Valência e número de oxidação formal do ânion haleto/óxido (ex.: Cl- possui carga -1)
valencia_formal_anion_adimensional  = -1;

% Constante geométrica de Madelung (dependente estritamente da geometria espacial do empacotamento)
% Para rede cristalina do tipo cúbica de faces centradas (CFC / NaCl): M = 1.74756
constante_madelung_geometria_rede = 1.74756;

% Raio iônico efetivo do cátion em picômetros (1 pm = 1e-12 metros)
raio_ionico_cation_picometros = 102.0; % Cátion Sódio (Na+)

% Raio iônico efetivo do ânion em picômetros (1 pm = 1e-12 metros)
raio_ionico_anion_picometros  = 181.0; % Ânion Cloreto (Cl-)

% Parâmetro empírico de compressibilidade e repulsão das camadas eletrônicas internas de Born-Mayer
parametro_repulsao_born_picometros = 34.5; % Valor clássico padrão em picômetros (0.345 Angstroms)

% ==============================================================================
% VALIDAÇÃO DOS PARÂMETROS FÍSICOS E CRISTALOGRÁFICOS
% ==============================================================================

if raio_ionico_cation_picometros <= 0.0 || raio_ionico_anion_picometros <= 0.0
    error('Os raios iônicos do cátion e do ânion devem ser estritamente positivos em picômetros.');
end
if parametro_repulsao_born_picometros <= 0.0
    error('O parâmetro fenomenológico de repulsão de Born deve ser estritamente positivo em picômetros.');
end
if parametro_repulsao_born_picometros >= (raio_ionico_cation_picometros + raio_ionico_anion_picometros)
    error('O parâmetro de repulsão rho deve ser estritamente menor que a distância interiônica de equilíbrio r_0.');
end

% ==============================================================================
% PROCESSAMENTO MATEMÁTICO E CONVERSÕES DE UNIDADES
% ==============================================================================

% Distância interiônica de equilíbrio entre centros de núcleos adjacentes em picômetros
distancia_interionica_equilibrio_picometros = raio_ionico_cation_picometros + raio_ionico_anion_picometros;

% Conversão da distância interiônica de picômetros para metros (escala padrão do SI)
distancia_interionica_equilibrio_metros = distancia_interionica_equilibrio_picometros * 1.0e-12;

% Conversão do parâmetro de repulsão de picômetros para metros
parametro_repulsao_born_metros = parametro_repulsao_born_picometros * 1.0e-12;

% Produto em módulo das cargas dos íons interagentes
produto_cargas_ionicas_absoluto = abs(valencia_formal_cation_adimensional * valencia_formal_anion_adimensional);

% Cálculo do numerador do termo eletrostático coulombiano molar
numerador_termo_eletrostatico = constante_avogadro_particulas_por_mol * ...
                                 constante_madelung_geometria_rede * ...
                                 produto_cargas_ionicas_absoluto * ...
                                 (carga_elementar_eletron_coulombs ^ 2);

% Cálculo do denominador do termo eletrostático molar (4 * pi * epsilon_0 * r_0)
denominador_termo_eletrostatico = 4.0 * pi * ...
                                   permissividade_eletrica_vacuo_farads_por_metro * ...
                                   distancia_interionica_equilibrio_metros;

% Energia de atração coulombiana estrita por mol de pares iônicos (Joules por mol)
energia_atrativa_coulombiana_joules_por_mol = numerador_termo_eletrostatico / denominador_termo_eletrostatico;

% Fator adimensional corretivo de atenuação por repulsão eletrônica quântica (1 - rho / r_0)
fator_correcao_repulsao_pauli_adimensional = 1.0 - (parametro_repulsao_born_metros / distancia_interionica_equilibrio_metros);

% Energia Reticular líquida U em Joules por mol (negativa pela convenção de liberação energética)
energia_reticular_total_joules_por_mol = -energia_atrativa_coulombiana_joules_por_mol * fator_correcao_repulsao_pauli_adimensional;

% Conversão da Energia Reticular de Joules/mol para quilojoules por mol (kJ/mol)
energia_reticular_total_quilojoules_por_mol = energia_reticular_total_joules_por_mol / 1000.0;

% ==============================================================================
% APRESENTAÇÃO DOS RESULTADOS DETALHADOS NO CONSOLE
% ==============================================================================

fprintf('\n==============================================================================\n');
fprintf('       FÍSICO-QUÍMICA: CÁLCULO DA ENERGIA RETICULAR (EQUAÇÃO DE BORN-MAYER)   \n');
fprintf('==============================================================================\n');
fprintf('  DADOS CRISTALOGRÁFICOS E IÔNICOS INFORMADOS:\n');
fprintf('    - Carga formal do cátion (z+):               %+2d\n', valencia_formal_cation_adimensional);
fprintf('    - Carga formal do ânion (z-):                %+2d\n', valencia_formal_anion_adimensional);
fprintf('    - Constante de Madelung da rede (M):         %10.5f\n', constante_madelung_geometria_rede);
fprintf('    - Raio do cátion (r+):                       %10.2f pm\n', raio_ionico_cation_picometros);
fprintf('    - Raio do ânion (r-):                        %10.2f pm\n', raio_ionico_anion_picometros);
fprintf('    - Distância interiônica de equilíbrio (r_0): %10.2f pm (%11.4e m)\n', ...
        distancia_interionica_equilibrio_picometros, distancia_interionica_equilibrio_metros);
fprintf('    - Parâmetro de repulsão de Born (rho):       %10.2f pm (%11.4e m)\n', ...
        parametro_repulsao_born_picometros, parametro_repulsao_born_metros);
fprintf('------------------------------------------------------------------------------\n');
fprintf('  COMPONENTES ENERGÉTICAS TERMODINÂMICAS:\n');
fprintf('    - Parcela Eletrostática Atrativa Pura:       %11.4f kJ/mol\n', energia_atrativa_coulombiana_joules_por_mol / 1000.0);
fprintf('    - Fator Adimensional de Correção Repulsiva:  %11.6f (redução de %5.2f%%)\n', ...
        fator_correcao_repulsao_pauli_adimensional, (1.0 - fator_correcao_repulsao_pauli_adimensional) * 100.0);
fprintf('------------------------------------------------------------------------------\n');
fprintf('  RESULTADO FINAL DA ENERGIA RETICULAR (U):\n');
fprintf('    - Energia Reticular Líquida:                %14.2f J/mol\n', energia_reticular_total_joules_por_mol);
fprintf('    - Energia Reticular Líquida:                %14.2f kJ/mol\n', energia_reticular_total_quilojoules_por_mol);
fprintf('==============================================================================\n\n');
