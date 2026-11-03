% ==============================================================================
% FARMACOCINÉTICA: MODELO DE DOSES REPETIDAS EM BOLUS INSTANTÂNEO (SEM ABSORÇÃO)
% ==============================================================================
%
% 1. O QUE CALCULA:
%    Este script simula a cinética plasmática de um fármaco administrado por via
%    intravenosa rápida (bolus intravascular instantâneo) em regime posológico de
%    múltiplas doses regulares a intervalos fixos (tau), sem fase prévia de absorção.
%
%    Formulação Matemática:
%    - Decaimento e Eliminação Metabólica de Primeira Ordem:
%      dC/dt = -k_e * C,  onde k_e é a constante de eliminação relacionada à meia-vida:
%      k_e = ln(2) / t_1/2.
%      A discretização horária simplificada do modelo clássico emprega:
%      C(t + 1) = C(t) - ( C(t) / (2 * t_1/2) )
%    - Injeção Instantânea Periódica:
%      A cada instante em que o contador atinge o intervalo de administração (tau),
%      adiciona-se a dose prescrita inteira de forma impulsiva:
%      C(t + 1) = C(t + 1) + Dose
%    - Suavização por Média Móvel Centralizada:
%      Calcula o nível plasmático médio contínuo com janela deslizante de W horas:
%      C_media(t) = (1 / (W + 1)) * somatorio de (t - W/2) até (t + W/2) de C(tau_i)
%    - Fator de Acúmulo no Estado de Equilíbrio (Steady-State):
%      R = C_media_estacionaria / Dose
%
% 2. O QUE ENTREGA:
%    - Janela gráfica com o perfil temporal da concentração plasmática:
%      * Curva serrilhada de pulsos de concentração instantânea (picos e vales).
%      * Curva suave da média móvel temporal evidenciando a ascensão até o platô.
%      * Eixo temporal em dias para facilidade de interpretação médica e farmacológica.
%    - Relatório quantitativo no console com:
%      * Razão entre Intervalo Posológico e Meia-Vida (tau / t_1/2).
%      * Nível plasmático médio atingido no estado estacionário.
%      * Fator de acúmulo farmacocinético no equilíbrio.
%
% 3. O QUE PRECISA TER (PRÉ-REQUISITOS E PARÂMETROS):
%    - Dose unitária administrada (em mg ou unidades arbitrárias).
%    - Intervalo entre administrações (tau em horas ou dias).
%    - Meia-vida biológica de eliminação plasmática do fármaco (t_1/2 em horas).
%    - Duração total do tratamento monitorado.
%
% 4. EXEMPLO DE USO E RESULTADO ESPERADO:
%    - Dose = 250 mg, Intervalo tau = 48 h (2 dias), Meia-vida = 168 h (7 dias),
%      Duração = 3600 h (150 dias), Janela média móvel = 288 h (12 dias).
%    - Resultado esperado:
%      Razão tau / t_1/2 = 0.2857 (alta frequência de administração relativa).
%      Concentração plasmática média no equilíbrio ≈ 1749.5 mg/L.
%      Fator de acúmulo no equilíbrio R ≈ 7.00 vezes a dose administrada.
% ==============================================================================

% Limpeza do ambiente de execução
clc;        % Limpa a tela do console de comandos
clear all;  % Remove todas as variáveis ativas do espaço de trabalho (workspace)
close all;  % Fecha todas as janelas de figuras abertas anteriormente

% ==============================================================================
% BLOCO DE PARÂMETROS FARMACOCINÉTICOS E POSOLÓGICOS DE ENTRADA
% ==============================================================================

% Quantidade em massa do princípio ativo administrada a cada injeção
dose_medicamento_administrada_mg = 250.0; % Unidade: miligramas (mg)

% Intervalo de tempo fixo entre doses consecutivas administradas ao paciente
intervalo_administracao_horas = 2.0 * 24.0; % 48 horas (administração a cada 2 dias)

% Meia-vida biológica de eliminação do fármaco pelo organismo (depuração hepática/renal)
tempo_meia_vida_eliminacao_horas = 7.0 * 24.0; % 168 horas (meia-vida de 7 dias)

% Duração total do período de monitoramento terapêutico clínico
tempo_total_monitoramento_horas = 150.0 * 24.0; % 3600 horas (150 dias de tratamento continuado)

% Largura da janela temporal para cálculo da média móvel suavizadora
janela_media_movel_horas = 12.0 * 24.0; % 288 horas (janela centralizada de 12 dias)

% ==============================================================================
% VALIDAÇÃO DOS PARÂMETROS FARMACOCINÉTICOS DE ENTRADA
% ==============================================================================

if dose_medicamento_administrada_mg <= 0.0
    error('A dose unitária administrada deve ser estritamente positiva em mg.');
end
if intervalo_administracao_horas <= 0.0
    error('O intervalo entre administrações deve ser estritamente positivo em horas.');
end
if tempo_meia_vida_eliminacao_horas <= 0.0
    error('O tempo de meia-vida de eliminação plasmática deve ser estritamente positivo em horas.');
end
if janela_media_movel_horas <= 0.0 || janela_media_movel_horas >= tempo_total_monitoramento_horas
    error('A janela de média móvel deve ser positiva e menor que o tempo total de simulação.');
end

% ==============================================================================
% INICIALIZAÇÃO DE ESTRUTURAS DE DADOS E VETORES TEMPORAIS
% ==============================================================================

% Passo temporal discreto da simulação numérica
passo_temporal_horas = 1.0; % Integração temporal horária

% Vetor discreto contendo todos os instantes temporais em horas
vetor_tempo_horas = 1.0 : passo_temporal_horas : tempo_total_monitoramento_horas;

% Número total de instantes de amostragem discretos
numero_total_instantes_horas = length(vetor_tempo_horas);

% Conversão do eixo temporal para dias para enriquecer a apresentação visual
vetor_tempo_dias = vetor_tempo_horas / 24.0;

% Pré-alocação do vetor de concentração plasmática instantânea inicializado com zeros
vetor_concentracao_plasmatica = zeros(1, numero_total_instantes_horas);

% Contador de tempo transcorrido desde a última dose para disparar a nova administração
contador_horas_intervalo = 0;

% ==============================================================================
% PROCESSAMENTO DINÂMICO DA CONCENTRAÇÃO PLASMÁTICA HORA A HORA
% ==============================================================================

% Laço temporal discreto calculando a evolução da concentração no sangue
for indice_hora_atual = 1 : (numero_total_instantes_horas - 1)
    
    % Fração eliminada no passo temporal de 1 hora pela taxa de depuração metabólica
    decaimento_metabolico_por_hora = (vetor_concentracao_plasmatica(indice_hora_atual) / tempo_meia_vida_eliminacao_horas) / 2.0;
    
    % Atualização preliminar da concentração após a eliminação biológica horária
    vetor_concentracao_plasmatica(indice_hora_atual + 1) = vetor_concentracao_plasmatica(indice_hora_atual) - decaimento_metabolico_por_hora;
    
    % Incremento do contador de horas para sincronização com o ciclo da prescrição médica
    contador_horas_intervalo = contador_horas_intervalo + 1;
    
    % Verificação de disparo de nova dose no término do intervalo programado
    if contador_horas_intervalo == intervalo_administracao_horas
        contador_horas_intervalo = 0; % Reinicia a contagem horária até o próximo pulso
        vetor_concentracao_plasmatica(indice_hora_atual + 1) = vetor_concentracao_plasmatica(indice_hora_atual + 1) + dose_medicamento_administrada_mg;
    end
end

% ==============================================================================
% PROCESSAMENTO DA MÉDIA MÓVEL SUAVIZADA (CURVA DO PLATÔ TERAPÊUTICO)
% ==============================================================================

% Pré-alocação do vetor de média móvel suavizada com zeros
vetor_media_movel_concentracao = zeros(1, numero_total_instantes_horas);

% Raio de abrangência lateral da janela centralizada
semi_largura_janela_horas = round(janela_media_movel_horas / 2.0);

% Início e término válidos para a convolução da média móvel sem extrapolar os limites do vetor
limite_inicial_convolucao = 1 + semi_largura_janela_horas;
limite_final_convolucao   = numero_total_instantes_horas - semi_largura_janela_horas;

% Cálculo da média móvel centralizada ponto a ponto
for indice_ponto_central = limite_inicial_convolucao : limite_final_convolucao
    intervalo_janela_inicio = indice_ponto_central - semi_largura_janela_horas;
    intervalo_janela_fim    = indice_ponto_central + semi_largura_janela_horas;
    
    vetor_media_movel_concentracao(indice_ponto_central) = ...
        mean(vetor_concentracao_plasmatica(intervalo_janela_inicio : intervalo_janela_fim));
end

% ==============================================================================
% AVALIAÇÃO DE PARÂMETROS QUANTITATIVOS DE EQUILÍBRIO (STEADY-STATE)
% ==============================================================================

% Razão dimensional entre o intervalo posológico e o tempo de meia-vida biológica
taxa_intervalo_por_meia_vida = intervalo_administracao_horas / tempo_meia_vida_eliminacao_horas;

% Extração do valor da concentração média móvel representativa no platô estável
indice_referencia_equilibrio = numero_total_instantes_horas - semi_largura_janela_horas;
concentracao_media_equilibrio_mg = vetor_media_movel_concentracao(indice_referencia_equilibrio);

% Fator adimensional de acúmulo no estado estacionário normalizado pela dose
fator_acumulo_concentracao_dose = concentracao_media_equilibrio_mg / dose_medicamento_administrada_mg;

% ==============================================================================
% APRESENTAÇÃO DOS RESULTADOS NO CONSOLE
% ==============================================================================

fprintf('\n==============================================================================\n');
fprintf('     FARMACOCINÉTICA: MODELO DE DOSES REPETIDAS EM BOLUS (SEM ABSORÇÃO)      \n');
fprintf('==============================================================================\n');
fprintf('  PARÂMETROS FARMACOLÓGICOS DO TRATAMENTO:\n');
fprintf('    - Dose Administrada a cada ciclo:            %10.2f mg\n', dose_medicamento_administrada_mg);
fprintf('    - Intervalo entre Administrações (tau):      %10.2f h (%4.1f dias)\n', ...
        intervalo_administracao_horas, intervalo_administracao_horas / 24.0);
fprintf('    - Meia-Vida de Eliminação Plasmática (t_1/2): %10.2f h (%4.1f dias)\n', ...
        tempo_meia_vida_eliminacao_horas, tempo_meia_vida_eliminacao_horas / 24.0);
fprintf('    - Janela da Média Móvel Suavizadora:         %10.2f h (%4.1f dias)\n', ...
        janela_media_movel_horas, janela_media_movel_horas / 24.0);
fprintf('    - Período Total Monitorado:                  %10.2f h (%4.1f dias)\n', ...
        tempo_total_monitoramento_horas, tempo_total_monitoramento_horas / 24.0);
fprintf('------------------------------------------------------------------------------\n');
fprintf('  INDICADORES TERAPÊUTICOS DE ESTADO DE EQUILÍBRIO:\n');
fprintf('    - Taxa Intervalo Posológico / Meia-Vida:     %10.4f\n', taxa_intervalo_por_meia_vida);
fprintf('    - Concentração Plasmática Média no Platô:    %10.2f mg/L\n', concentracao_media_equilibrio_mg);
fprintf('    - Fator de Acúmulo no Estado Estacionário:   %10.4f x Dose\n', fator_acumulo_concentracao_dose);
fprintf('==============================================================================\n\n');

% ==============================================================================
% VISUALIZAÇÃO GRÁFICA DO PERFIL TEMPORAL
% ==============================================================================

figura_perfil_farmacocinetico = figure(1);
set(figura_perfil_farmacocinetico, 'Name', 'Perfil Farmacocinetico - Bolus sem Absorcao');
hold on;

% Plotagem da curva serrilhada de pulsos instantâneos
plot(vetor_tempo_dias, vetor_concentracao_plasmatica, 'Color', [0.1, 0.4, 0.8], 'LineWidth', 1.2);

% Plotagem da curva suave da média móvel
indices_validos_media = limite_inicial_convolucao : limite_final_convolucao;
plot(vetor_tempo_dias(indices_validos_media), ...
     vetor_media_movel_concentracao(indices_validos_media), ...
     'Color', [0.85, 0.2, 0.1], 'LineWidth', 2.5);

grid on;
title('Perfil de Concentração Plasmática em Regime de Doses Repetidas (Bolus IV)', ...
      'FontSize', 13, 'FontWeight', 'bold');
xlabel('Tempo decorrido de tratamento (dias)', 'FontSize', 11);
ylabel('Concentração Plasmática (mg / L)', 'FontSize', 11);

identificador_legenda_farmaco = legend('Concentração Instantânea (Pulsos e Eliminação)', ...
                                       'Média Móvel Suavizada (Platô Estacionário)', ...
                                       'Location', 'southeast');
set(identificador_legenda_farmaco, 'FontSize', 10);
hold off;
