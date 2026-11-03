% ==============================================================================
% FARMACOCINÉTICA: ESTUDO PARAMÉTRICO DO ACÚMULO EM FUNÇÃO DA RAZÃO INTERVALO / MEIA-VIDA
% ==============================================================================
%
% 1. O QUE CALCULA:
%    Este script realiza uma varredura paramétrica sistemática para quantificar como
%    a relação entre a frequência de administração (intervalo entre doses, tau)
%    e a velocidade de depuração do organismo (meia-vida de eliminação, t_1/2)
%    governa o acúmulo de fármaco no estado de equilíbrio terapêutico (steady-state).
%
%    Formulação Matemática:
%    - Para cada ponto de ensaio 'k', define-se um intervalo posológico:
%      tau(k) = k * delta_horas
%    - Simula-se a evolução temporal da concentração plasmática sob injeções repetidas
%      em bolus por um período extenso (120 dias) até que o sistema atinja o regime
%      permanente estático.
%    - Calcula-se a média móvel final da concentração no platô (C_plato) e determina-se
%      o Fator de Acúmulo Normalizado pela Dose:
%      Fator_Acumulo = C_plato / Dose
%    - Compara-se a curva de pontos discretos simulados com a relação hiperbólica teórica:
%      Fator_Acumulo_Teorico(r) = 2 / r,  onde r = tau / t_1/2.
%
% 2. O QUE ENTREGA:
%    - Janela gráfica estilizada com a curva do Fator de Acúmulo no Estado Estacionário
%      em função da razão adimensional (Intervalo Posológico / Meia-Vida):
%      * Marcadores circulares com linha conectora para os pontos numéricos simulados.
%      * Curva contínua para o modelo teórico de conservação de massa em regime estacionário.
%    - Tabela comparativa impressa no console com valores simulados, teóricos e desvios.
%
% 3. O QUE PRECISA TER (PRÉ-REQUISITOS E PARÂMETROS):
%    - Dose unitária de referência (mg).
%    - Meia-vida biológica de eliminação plasmática fixada (horas).
%    - Conjunto discreto de intervalos a serem testados na varredura (ex.: de 12h a 120h).
%    - Período de monitoramento longo o suficiente para que t_total >> 5 * t_1/2.
%
% 4. EXEMPLO DE USO E RESULTADO ESPERADO:
%    - Dose = 100 mg, Meia-vida = 48 h (2 dias), Intervalos de 12 h a 120 h (passo de 12 h).
%    - Resultado esperado:
%      Para r = 0.25 (administração muito frequente), acúmulo massivo ≈ 8 vezes a dose.
%      Para r = 1.00 (administração igual à meia-vida), acúmulo de exatamente 2 vezes a dose.
%      Para r = 2.00 (administração a cada duas meias-vidas), acúmulo de 1 vez a dose.
% ==============================================================================

% Limpeza do ambiente de execução
clc;        % Limpa a tela do console de comandos
clear all;  % Remove todas as variáveis ativas do espaço de trabalho (workspace)
close all;  % Fecha todas as janelas de figuras abertas anteriormente

% ==============================================================================
% BLOCO DE PARÂMETROS GERAIS DA SIMULAÇÃO FARMACOCINÉTICA
% ==============================================================================

% Dose única padronizada administrada a cada ciclo posológico
dose_medicamento_administrada_mg = 100.0; % Unidade: miligramas (mg)

% Tempo de meia-vida biológica de eliminação plasmática constante do composto
tempo_meia_vida_eliminacao_horas = 2.0 * 24.0; % 48 horas (meia-vida de 2 dias completos)

% Duração total do monitoramento contínuo para assegurar convergência ao estado de equilíbrio
tempo_total_simulacao_horas = 120.0 * 24.0; % 2880 horas (120 dias, correspondendo a 60 meias-vidas)

% Janela temporal da média móvel para extração do nível médio estável
janela_media_movel_horas = 12.0 * 24.0; % 288 horas (12 dias para suavização do platô)

% Número total de cenários posológicos avaliados na varredura paramétrica
numero_cenarios_posologicos_varredura = 10;

% Incremento de tempo de intervalo a cada cenário avaliado (em horas)
passo_intervalo_administracao_horas = 12.0;

% ==============================================================================
% VALIDAÇÃO DOS PARÂMETROS DA VARREDURA
% ==============================================================================

if dose_medicamento_administrada_mg <= 0.0
    error('A dose unitária administrada deve ser estritamente positiva em mg.');
end
if tempo_meia_vida_eliminacao_horas <= 0.0
    error('O tempo de meia-vida de eliminação plasmática deve ser estritamente positivo em horas.');
end
if numero_cenarios_posologicos_varredura < 1
    error('O número de cenários posológicos avaliados deve ser no mínimo 1.');
end
if passo_intervalo_administracao_horas <= 0.0
    error('O passo de incremento de intervalo deve ser estritamente positivo em horas.');
end
if janela_media_movel_horas <= 0.0 || janela_media_movel_horas >= tempo_total_simulacao_horas
    error('A janela de média móvel deve ser positiva e menor que o tempo total simulado.');
end

% ==============================================================================
% PRÉ-ALOCAÇÃO PRECISA DE VETORES DE RESULTADOS
% ==============================================================================

% Vetor para armazenar a razão adimensional (tau / t_1/2) de cada cenário avaliado
vetor_razao_intervalo_meia_vida = zeros(1, numero_cenarios_posologicos_varredura);

% Vetor para armazenar o fator de acúmulo medido numericamente (C_plato / Dose)
vetor_fator_acumulo_simulado = zeros(1, numero_cenarios_posologicos_varredura);

% Vetor para armazenar a previsão teórica hiperbólica de acúmulo
vetor_fator_acumulo_teorico = zeros(1, numero_cenarios_posologicos_varredura);

% Metade da janela para indexação centralizada
semi_largura_janela_horas = round(janela_media_movel_horas / 2.0);

% Índice final estável onde a concentração média no platô é extraída
indice_amostra_referencia_estavel = tempo_total_simulacao_horas - janela_media_movel_horas;

% ==============================================================================
% LAÇO PRINCIPAL DA VARREDURA PARAMÉTRICA DOS INTERVALOS DE ADMINISTRAÇÃO
% ==============================================================================

for indice_cenario = 1 : numero_cenarios_posologicos_varredura
    
    % Intervalo posológico fixado para o cenário atual
    intervalo_administracao_horas = indice_cenario * passo_intervalo_administracao_horas;
    
    % Reinicialização explícita do vetor de concentração plasmática para o novo ensaio
    vetor_concentracao_plasmatica = zeros(1, tempo_total_simulacao_horas);
    
    % Contador horário interno para disparar as administrações
    contador_horas_intervalo = 0;
    
    % Integração horária do perfil de concentração para o cenário corrente
    for indice_hora_atual = 1 : (tempo_total_simulacao_horas - 1)
        
        % Eliminação fracionária horária de primeira ordem
        perda_metabolica_horaria = (vetor_concentracao_plasmatica(indice_hora_atual) / tempo_meia_vida_eliminacao_horas) / 2.0;
        vetor_concentracao_plasmatica(indice_hora_atual + 1) = vetor_concentracao_plasmatica(indice_hora_atual) - perda_metabolica_horaria;
        
        % Incremento e verificação da injeção de nova dose
        contador_horas_intervalo = contador_horas_intervalo + 1;
        if contador_horas_intervalo == intervalo_administracao_horas
            contador_horas_intervalo = 0;
            vetor_concentracao_plasmatica(indice_hora_atual + 1) = vetor_concentracao_plasmatica(indice_hora_atual + 1) + dose_medicamento_administrada_mg;
        end
    end
    
    % Cálculo da média móvel centralizada na região final do estado de equilíbrio
    vetor_media_movel_concentracao = zeros(1, tempo_total_simulacao_horas);
    for indice_ponto_tempo = (1 + semi_largura_janela_horas) : indice_amostra_referencia_estavel
        inicio_janela = indice_ponto_tempo - semi_largura_janela_horas;
        fim_janela    = indice_ponto_tempo + semi_largura_janela_horas;
        vetor_media_movel_concentracao(indice_ponto_tempo) = mean(vetor_concentracao_plasmatica(inicio_janela : fim_janela));
    end
    
    % Razão adimensional entre intervalo e meia-vida (tau / t_1/2)
    razao_tau_meia_vida = intervalo_administracao_horas / tempo_meia_vida_eliminacao_horas;
    vetor_razao_intervalo_meia_vida(indice_cenario) = razao_tau_meia_vida;
    
    % Concentração média representativa obtida no platô do regime permanente
    concentracao_plato_estacionario = vetor_media_movel_concentracao(indice_amostra_referencia_estavel);
    
    % Fator de acúmulo simulado normalizado pela dose
    vetor_fator_acumulo_simulado(indice_cenario) = concentracao_plato_estacionario / dose_medicamento_administrada_mg;
    
    % Modelo teórico de conservação de massa em regime permanente (Fator = 2 / r)
    vetor_fator_acumulo_teorico(indice_cenario) = 2.0 / razao_tau_meia_vida;
end

% ==============================================================================
% APRESENTAÇÃO DOS RESULTADOS NO CONSOLE
% ==============================================================================

fprintf('\n==============================================================================\n');
fprintf('     FARMACOCINÉTICA: VARREDURA DA RELAÇÃO INTERVALO / MEIA-VIDA vs ACÚMULO   \n');
fprintf('==============================================================================\n');
fprintf('  Dose Padrão: %6.1f mg | Meia-Vida: %6.1f h (%3.1f dias) | Janela Móvel: %6.1f h\n', ...
        dose_medicamento_administrada_mg, tempo_meia_vida_eliminacao_horas, ...
        tempo_meia_vida_eliminacao_horas / 24.0, janela_media_movel_horas);
fprintf('------------------------------------------------------------------------------\n');
fprintf('  Cenário | Intervalo (h) | Razão tau / t_1/2 | Acúmulo Simulado | Acúmulo Teórico\n');
fprintf('------------------------------------------------------------------------------\n');
for indice_cenario = 1 : numero_cenarios_posologicos_varredura
    fprintf('    %2d    |   %5.1f h     |      %6.4f       |     %7.4f x      |    %7.4f x\n', ...
            indice_cenario, ...
            indice_cenario * passo_intervalo_administracao_horas, ...
            vetor_razao_intervalo_meia_vida(indice_cenario), ...
            vetor_fator_acumulo_simulado(indice_cenario), ...
            vetor_fator_acumulo_teorico(indice_cenario));
end
fprintf('==============================================================================\n\n');

% ==============================================================================
% VISUALIZAÇÃO GRÁFICA DA CURVA DE ACÚMULO FARMACOCINÉTICO
% ==============================================================================

% Geração de uma malha fina teórica contínua para comparação visual
razao_teorica_fina = linspace(min(vetor_razao_intervalo_meia_vida), max(vetor_razao_intervalo_meia_vida), 200);
fator_acumulo_teorico_fino = 2.0 ./ razao_teorica_fina;

figura_acumulo = figure(1);
set(figura_acumulo, 'Name', 'Varredura de Acumulo Farmacocinetico');
hold on;

% Curva analítica teórica contínua
plot(razao_teorica_fina, fator_acumulo_teorico_fino, 'r--', 'LineWidth', 2.0);

% Pontos medidos da simulação numérica horária
plot(vetor_razao_intervalo_meia_vida, vetor_fator_acumulo_simulado, ...
     'bo-', 'LineWidth', 2.0, 'MarkerSize', 8, 'MarkerFaceColor', [0.2, 0.6, 1.0]);

grid on;
title('Fator de Acúmulo Plasmático no Equilíbrio vs Relação Intervalo / Meia-Vida', ...
      'FontSize', 13, 'FontWeight', 'bold');
xlabel('Razão Adimensional entre Intervalo de Dose e Meia-Vida (\tau / t_{1/2})', 'FontSize', 11);
ylabel('Fator de Acúmulo no Platô (\bar{C}_{ss} / Dose)', 'FontSize', 11);

identificador_legenda_acumulo = legend('Modelo Teórico Hiperbólico: 2 / (\tau / t_{1/2})', ...
                                      'Simulação Numérica Discreta de Doses', ...
                                      'Location', 'northeast');
set(identificador_legenda_acumulo, 'FontSize', 10);
hold off;
