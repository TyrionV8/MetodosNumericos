% ==============================================================================
% FARMACOCINÉTICA: MODELO BICOMPARTIMENTAL COM ABSORÇÃO GRADUAL E ELIMINAÇÃO
% ==============================================================================
%
% 1. O QUE CALCULA:
%    Este script simula a cinética plasmática de um fármaco administrado por via
%    extravascular (oral, intramuscular ou subcutânea), onde a substância não
%    entra imediatamente na circulação sistêmica, mas passa por um compartimento
%    de depósito/absorção preliminar.
%
%    Formulação Matemática:
%    - Compartimento de Absorção (Dose Não Absorvida, D_NA):
%      A cada intervalo de administração fixo (tau), o reservatório de depósito
%      recebe uma nova dose:
%      D_NA(t) = D_NA(t) + Dose
%      A droga é transferida ao sangue a uma taxa fracionária horária de absorção k_a:
%      d(D_NA)/dt = -k_a * D_NA
%      Delta_Absorvido = D_NA * k_a * Delta_t
%    - Compartimento Plasmático Sistêmico (C):
%      Recebe o aporte absorvido e sofre metabolização/eliminação de primeira ordem:
%      C(t + 1) = C(t) - ( C(t) / (2 * t_1/2) ) + Delta_Absorvido
%    - Suavização por Média Móvel:
%      Cálculo do nível médio com janela deslizante de W horas (janela_media_movel_horas),
%      resolvendo o bug histórico do script original onde a variável 'horasmedia'
%      não estava definida, causando parada inesperada na execução.
%    - Razão e Fator de Acúmulo no Equilíbrio (Steady-State):
%      R = C_media_estacionaria / Dose
%
% 2. O QUE ENTREGA:
%    - Janela gráfica didática com:
%      * Curva contínua de concentração plasmática instantânea C(t).
%      * Curva da média móvel revelando o platô de equilíbrio estacionário.
%      * Eixo das abscissas formatado em dias de tratamento continuado.
%    - Relatório quantitativo no console exibindo parâmetros cinéticos e taxa de acúmulo.
%
% 3. O QUE PRECISA TER (PRÉ-REQUISITOS E PARÂMETROS):
%    - Dose unitária administrada (mg).
%    - Constante de absorção fracionária horária (k_a).
%    - Intervalo de administração (tau em horas).
%    - Meia-vida biológica de eliminação (t_1/2 em horas).
%    - Período total de monitoramento e largura da janela de média móvel.
%
% 4. EXEMPLO DE USO E RESULTADO ESPERADO:
%    - Dose = 100 mg, Taxa de Absorção = 1/24 h^-1 (absorção lenta ao longo de 1 dia),
%      Intervalo tau = 96 h (4 dias), Meia-vida = 168 h (7 dias), Duração = 90 dias (2160 h).
%    - Correção do bug legado: 'horasmedia' substituída por 'janela_media_movel_horas = 168 h'.
%    - Resultado esperado:
%      Taxa Intervalo / Meia-Vida = 0.5714.
%      Concentração plasmática média no platô ≈ 349.96 mg/L.
%      Fator de acúmulo ≈ 3.50 vezes a dose administrada.
% ==============================================================================

% Limpeza do ambiente de execução
clc;        % Limpa a tela do console de comandos
clear all;  % Remove todas as variáveis ativas do espaço de trabalho (workspace)
close all;  % Fecha todas as janelas de figuras abertas anteriormente

% ==============================================================================
% BLOCO DE PARÂMETROS FARMACOCINÉTICOS E DE REGIME TERAPÊUTICO
% ==============================================================================

% Duração total do monitoramento clínico do paciente em horas
tempo_total_monitoramento_horas = 90.0 * 24.0; % 2160 horas (90 dias de acompanhamento)

% Quantidade em massa do princípio ativo administrada a cada ciclo posológico
dose_medicamento_administrada_mg = 100.0; % Unidade: miligramas (mg)

% Taxa fracionária de absorção do compartimento de depósito para a circulação sistêmica por hora
taxa_absorcao_fracionaria_por_hora = 1.0 / 24.0; % Cinética de liberação lenta (aprox. 4.17% por hora)

% Intervalo de tempo programado entre doses consecutivas
intervalo_administracao_horas = 4.0 * 24.0; % 96 horas (administração a cada 4 dias)

% Meia-vida biológica de eliminação plasmática (depuração metabólica e renal)
tempo_meia_vida_eliminacao_horas = 7.0 * 24.0; % 168 horas (meia-vida biológica de 7 dias)

% Janela temporal de suavização por média móvel (resolvendo a variável 'horasmedia' não declarada no código legado)
janela_media_movel_horas = 7.0 * 24.0; % 168 horas (janela móvel correspondente a 1 semana inteira)

% ==============================================================================
% VALIDAÇÃO DOS PARÂMETROS FARMACOCINÉTICOS DE ENTRADA
% ==============================================================================

if dose_medicamento_administrada_mg <= 0.0
    error('A dose unitária administrada deve ser estritamente positiva em mg.');
end
if taxa_absorcao_fracionaria_por_hora <= 0.0 || taxa_absorcao_fracionaria_por_hora > 1.0
    error('A taxa fracionária de absorção horária k_a deve estar no intervalo (0, 1].');
end
if intervalo_administracao_horas <= 0.0
    error('O intervalo entre administrações deve ser estritamente positivo em horas.');
end
if tempo_meia_vida_eliminacao_horas <= 0.0
    error('O tempo de meia-vida de eliminação plasmática deve ser estritamente positivo em horas.');
end
if janela_media_movel_horas <= 0.0 || janela_media_movel_horas >= tempo_total_monitoramento_horas
    error('A janela de média móvel deve ser positiva e menor que o tempo total de monitoramento.');
end

% ==============================================================================
% INICIALIZAÇÃO DE VARIÁVEIS DE CONTROLE E VETORES DINÂMICOS
% ==============================================================================

% Vetor temporal discreto de amostragem em horas
vetor_tempo_horas = 1.0 : 1.0 : tempo_total_monitoramento_horas;
numero_total_instantes_horas = length(vetor_tempo_horas);

% Vetor temporal convertido para dias corridos para exibição gráfica
vetor_tempo_dias = (vetor_tempo_horas - 1.0) / 24.0;

% Vetor de concentração plasmática sistêmica inicializado em repouso (zero no instante inicial)
vetor_concentracao_plasmatica = zeros(1, numero_total_instantes_horas);

% Quantidade de fármaco retida no compartimento de absorção inicializada com a primeira dose em t = 0
dose_deposito_nao_absorvida_mg = dose_medicamento_administrada_mg;

% Contador de tempo para controle dos ciclos de reinjeção/reingestão
contador_horas_intervalo = 0;

% ==============================================================================
% PROCESSAMENTO DINÂMICO DA FARMACOCINÉTICA HORÁRIA COM ABSORÇÃO E ELIMINAÇÃO
% ==============================================================================

for indice_hora_atual = 1 : (numero_total_instantes_horas - 1)
    
    % Incremento do temporizador de ciclo posológico
    contador_horas_intervalo = contador_horas_intervalo + 1;
    
    % Verificação de chegada de nova administração periódica
    if contador_horas_intervalo == intervalo_administracao_horas
        dose_deposito_nao_absorvida_mg = dose_deposito_nao_absorvida_mg + dose_medicamento_administrada_mg;
        contador_horas_intervalo = 0; % Reinicia o ciclo até o próximo aporte
    end
    
    % Manutenção do estoque de concentração prévia
    vetor_concentracao_plasmatica(indice_hora_atual + 1) = vetor_concentracao_plasmatica(indice_hora_atual);
    
    % Subtração da perda por depuração metabólica e excreção renal de primeira ordem
    perda_metabolica_horaria = (vetor_concentracao_plasmatica(indice_hora_atual) / 2.0) / tempo_meia_vida_eliminacao_horas;
    vetor_concentracao_plasmatica(indice_hora_atual + 1) = vetor_concentracao_plasmatica(indice_hora_atual + 1) - perda_metabolica_horaria;
    
    % Transferência de massa por absorção do reservatório para a corrente sanguínea
    massa_farmaco_absorvida_horaria = dose_deposito_nao_absorvida_mg * taxa_absorcao_fracionaria_por_hora;
    vetor_concentracao_plasmatica(indice_hora_atual + 1) = vetor_concentracao_plasmatica(indice_hora_atual + 1) + massa_farmaco_absorvida_horaria;
    
    % Atualização do saldo remanescente no compartimento de depósito
    dose_deposito_nao_absorvida_mg = dose_deposito_nao_absorvida_mg - massa_farmaco_absorvida_horaria;
end

% ==============================================================================
% PROCESSAMENTO DA MÉDIA MÓVEL (ELIMINAÇÃO DEFINITIVA DO BUG DE 'horasmedia')
% ==============================================================================

% Vetor pré-alocado para os valores da média móvel suavizada
vetor_media_movel_concentracao = zeros(1, numero_total_instantes_horas);

% Metade da largura da janela para indexação centralizada
semi_largura_janela_horas = round(janela_media_movel_horas / 2.0);

% Convolução discreta da média móvel ponto a ponto ao longo do tempo
for indice_ponto_tempo = 1 : (numero_total_instantes_horas - semi_largura_janela_horas)
    if indice_ponto_tempo <= semi_largura_janela_horas
        % Período de aquecimento inicial (janela progressiva cumulativa)
        vetor_media_movel_concentracao(indice_ponto_tempo) = sum(vetor_concentracao_plasmatica(1 : indice_ponto_tempo)) / indice_ponto_tempo;
    else
        % Janela móvel simétrica centralizada
        inicio_janela = indice_ponto_tempo - semi_largura_janela_horas;
        fim_janela    = indice_ponto_tempo + semi_largura_janela_horas;
        vetor_media_movel_concentracao(indice_ponto_tempo) = sum(vetor_concentracao_plasmatica(inicio_janela : fim_janela)) / (2.0 * semi_largura_janela_horas + 1.0);
    end
end

% ==============================================================================
% AVALIAÇÃO QUANTITATIVA DO EQUILÍBRIO ESTACIONÁRIO (STEADY-STATE)
% ==============================================================================

% Razão entre o intervalo de administração e o tempo de meia-vida biológica
taxa_intervalo_por_meia_vida = intervalo_administracao_horas / tempo_meia_vida_eliminacao_horas;

% Cálculo da concentração representativa do platô terapêutico final utilizando a janela corrigida
indice_amostra_final_estavel = numero_total_instantes_horas - semi_largura_janela_horas;
concentracao_media_plato_final = vetor_media_movel_concentracao(indice_amostra_final_estavel);

% Fator de acúmulo normalizado pela dose administrada
fator_acumulo_concentracao_dose_final = concentracao_media_plato_final / dose_medicamento_administrada_mg;

% ==============================================================================
% APRESENTAÇÃO DOS RESULTADOS NO CONSOLE
% ==============================================================================

fprintf('\n==============================================================================\n');
fprintf('     FARMACOCINÉTICA: MODELO COM FASE DE ABSORÇÃO GRADUAL E METABOLISMO       \n');
fprintf('==============================================================================\n');
fprintf('  PARÂMETROS FARMACOLÓGICOS INFORMADOS:\n');
fprintf('    - Dose Administrada a cada Ciclo:             %10.2f mg\n', dose_medicamento_administrada_mg);
fprintf('    - Taxa de Absorção Fracionária Horária (k_a): %10.4f h^-1\n', taxa_absorcao_fracionaria_por_hora);
fprintf('    - Intervalo entre Administrações (tau):       %10.2f h (%4.1f dias)\n', ...
        intervalo_administracao_horas, intervalo_administracao_horas / 24.0);
fprintf('    - Meia-Vida de Eliminação Plasmática (t_1/2):  %10.2f h (%4.1f dias)\n', ...
        tempo_meia_vida_eliminacao_horas, tempo_meia_vida_eliminacao_horas / 24.0);
fprintf('    - Janela da Média Móvel (corrigida):          %10.2f h (%4.1f dias)\n', ...
        janela_media_movel_horas, janela_media_movel_horas / 24.0);
fprintf('------------------------------------------------------------------------------\n');
fprintf('  INDICADORES TERAPÊUTICOS DE ESTADO DE EQUILÍBRIO:\n');
fprintf('    - Taxa Intervalo Posológico / Meia-Vida:      %10.4f\n', taxa_intervalo_por_meia_vida);
fprintf('    - Concentração Plasmática Média no Platô:     %10.2f mg/L\n', concentracao_media_plato_final);
fprintf('    - Fator de Acúmulo no Estado Estacionário:    %10.4f x Dose\n', fator_acumulo_concentracao_dose_final);
fprintf('==============================================================================\n\n');

% ==============================================================================
% VISUALIZAÇÃO GRÁFICA DO PERFIL TEMPORAL
% ==============================================================================

figura_perfil_absorcao = figure(1);
set(figura_perfil_absorcao, 'Name', 'Perfil Farmacocinetico - Modelo com Absorcao');
hold on;

% Curva da concentração instantânea
plot(vetor_tempo_dias, vetor_concentracao_plasmatica, 'Color', [0.1, 0.6, 0.3], 'LineWidth', 1.2);

% Curva da média móvel suavizada do platô
indices_validos_plato = 1 : (numero_total_instantes_horas - semi_largura_janela_horas);
plot(vetor_tempo_dias(indices_validos_plato), ...
     vetor_media_movel_concentracao(indices_validos_plato), ...
     'Color', [0.85, 0.1, 0.1], 'LineWidth', 2.5);

grid on;
title('Cinética Plasmática com Absorção Gradual e Doses Múltiplas Extravasculares', ...
      'FontSize', 13, 'FontWeight', 'bold');
xlabel('Tempo decorrido de tratamento (dias)', 'FontSize', 11);
ylabel('Concentração Plasmática Sistêmica (mg / L)', 'FontSize', 11);

identificador_legenda_absorcao = legend('Concentração Instantânea (Absorção + Eliminação)', ...
                                        'Média Móvel Suavizada (Platô Estacionário)', ...
                                        'Location', 'southeast');
set(identificador_legenda_absorcao, 'FontSize', 10);
hold off;
