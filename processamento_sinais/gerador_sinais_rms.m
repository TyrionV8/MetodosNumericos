% ==============================================================================
% PROCESSAMENTO DE SINAIS: GERADOR DE FORMAS DE ONDA E CÁLCULO DE VALOR EFICAZ (RMS)
% ==============================================================================
%
% 1. O QUE CALCULA:
%    Este script realiza a síntese temporal e a análise estatística de três formas
%    de onda fundamentais da engenharia elétrica e acústica:
%    - Onda Senoidal pura com nível contínuo (offset DC): x(t) = V_dc + A * sen(2*pi*f*t)
%    - Onda Quadrada simétrica periódica: x(t) = V_dc + A * sinal(sen(2*pi*f*t))
%    - Onda Triangular periódica simétrica: x(t) = V_dc + A * tri(2*pi*f*t)
%
%    Para cada sinal, calcula-se numericamente:
%    a) O Valor Médio Temporal (Componente Contínua DC):
%       V_medio = (1 / T) * integral de 0 a T de x(t) dt ≈ (1 / N) * somatorio(x_i)
%    b) O Valor Eficaz ou RMS (Root Mean Square):
%       V_rms = sqrt( (1 / T) * integral de 0 a T de x(t)^2 dt ) ≈ sqrt( (1 / N) * somatorio(x_i^2) )
%
%    Compara analiticamente os valores medidos discretizados com as fórmulas exatas:
%    - Senoide:    V_medio_teo = V_dc | V_rms_teo = sqrt( V_dc^2 + (A^2)/2 )
%    - Quadrada:   V_medio_teo = V_dc | V_rms_teo = sqrt( V_dc^2 + A^2 )
%    - Triangular: V_medio_teo = V_dc | V_rms_teo = sqrt( V_dc^2 + (A^2)/3 )
%
%    Garante conformidade com o Teorema de Nyquist-Shannon (taxa de amostragem
%    f_amostragem >> 2 * f_maxima) para mitigar distorções por aliasing (serrilhamento).
%
% 2. O QUE ENTREGA:
%    - 4 janelas gráficas didáticas:
%      * Figura 1: Forma de onda senoidal com linha de média, linha RMS e linha zero.
%      * Figura 2: Forma de onda quadrada com referências de média e RMS.
%      * Figura 3: Forma de onda triangular com referências de média e RMS.
%      * Figura 4: Comparação simultânea com as três ondas sobrepostas (handle L4 corrigido).
%    - Tabela comparativa impressa no console com valores teóricos, medidos e erros relativos.
%
% 3. O QUE PRECISA TER (PRÉ-REQUISITOS E PARÂMETROS):
%    - Nenhuma toolbox obrigatória (utiliza formulações analíticas universais com
%      compatibilidade total para GNU Octave e MATLAB).
%    - Parâmetros configuráveis de frequência, amplitude, offset e janela de amostragem.
%
% 4. EXEMPLO DE USO E RESULTADO ESPERADO:
%    - Executar diretamente no terminal do Octave/MATLAB:
%      >> gerador_sinais_rms
%    - Resultado esperado:
%      Exibição dos 4 gráficos estilizados e tabela quantitativa no terminal
%      comprovando erro inferior a 0.5% entre medições discretas e teoria analítica.
% ==============================================================================

% Limpeza do ambiente de execução
clc;        % Limpa a tela do console de comandos
clear all;  % Remove todas as variáveis ativas do espaço de trabalho (workspace)
close all;  % Fecha todas as janelas de figuras abertas anteriormente

% ==============================================================================
% BLOCO DE PARÂMETROS DE ENTRADA E CONFIGURAÇÃO DE AMOSTRAGEM
% ==============================================================================

% Parâmetros do sinal senoidal
frequencia_sinal_senoidal_hertz   = 1000.0; % Frequência fundamental da senoide em Hertz (1 kHz)
amplitude_sinal_senoidal_volts     = 2.5;    % Amplitude de pico da senoide em Volts
offset_tensao_dc_senoidal_volts    = 0.5;    % Deslocamento de nível contínuo (tensão DC) em Volts

% Parâmetros do sinal quadrado
frequencia_sinal_quadrado_hertz    = 1000.0; % Frequência fundamental da onda quadrada em Hertz
amplitude_sinal_quadrado_volts     = 2.0;    % Amplitude do degrau da onda quadrada em Volts
offset_tensao_dc_quadrado_volts    = 1.0;    % Nível contínuo da onda quadrada em Volts

% Parâmetros do sinal triangular
frequencia_sinal_triangular_hertz  = 1000.0; % Frequência fundamental da onda triangular em Hertz
amplitude_sinal_triangular_volts   = 3.0;    % Amplitude de pico da onda triangular em Volts
offset_tensao_dc_triangular_volts  = 0.0;    % Nível contínuo da onda triangular em Volts

% Janela temporal de observação (duração em segundos)
tempo_inicial_segundos             = 0.0;    % Instante temporal inicial da simulação
tempo_final_segundos               = 0.003;  % Instante temporal final (3 milissegundos = 3 períodos de 1 kHz)

% Validação dos parâmetros físicos de entrada
if frequencia_sinal_senoidal_hertz <= 0 || frequencia_sinal_quadrado_hertz <= 0 || frequencia_sinal_triangular_hertz <= 0
    error('Todas as frequências dos sinais devem ser estritamente positivas em Hertz.');
end
if tempo_final_segundos <= tempo_inicial_segundos
    error('O tempo final de simulação deve ser estritamente maior que o tempo inicial.');
end

% Determinação da frequência máxima para o critério de Nyquist
frequencia_maxima_sinais_hertz     = max([frequencia_sinal_senoidal_hertz, ...
                                          frequencia_sinal_quadrado_hertz, ...
                                          frequencia_sinal_triangular_hertz]);

% Taxa de Nyquist estrita teórica (frequência mínima para evitar aliasing)
frequencia_nyquist_minima_hertz    = 2.0 * frequencia_maxima_sinais_hertz;

% Fator de sobreamostragem (oversampling) para visualização analógica suave e alta precisão de integração
fator_sobreamostragem_suave        = 200; % 200 amostras discretas para cada período de oscilação

% Taxa de amostragem efetiva adotada (muito superior ao limite de Nyquist)
taxa_amostragem_efetiva_hertz      = fator_sobreamostragem_suave * frequencia_maxima_sinais_hertz;

% Número total de pontos de amostragem no intervalo temporal observado
duracao_total_janela_segundos      = tempo_final_segundos - tempo_inicial_segundos;
numero_pontos_amostras             = round(duracao_total_janela_segundos * taxa_amostragem_efetiva_hertz);

% Garantia de que o vetor temporal contenha no mínimo 1000 amostras para gráficos contínuos
if numero_pontos_amostras < 1000
    numero_pontos_amostras = 1000;
end

% Vetor temporal linearmente espaçado representando o fluxo contínuo do tempo
vetor_tempo_segundos = linspace(tempo_inicial_segundos, tempo_final_segundos, numero_pontos_amostras);

% Vetor horizontal com nível zero para servir como referência nos eixos gráficos
vetor_linha_zero_referencia = zeros(1, numero_pontos_amostras);

% ==============================================================================
% PROCESSAMENTO NUMÉRICO 1: SÍNTESE E ANÁLISE DA ONDA SENOIDAL
% ==============================================================================

% Geração matemática da onda senoidal com componente contínua
sinal_onda_senoidal_volts = offset_tensao_dc_senoidal_volts + ...
    amplitude_sinal_senoidal_volts * sin(2.0 * pi * frequencia_sinal_senoidal_hertz * vetor_tempo_segundos);

% Cálculo do valor médio numérico medido no conjunto de amostras
valor_medio_medido_senoidal_volts = mean(sinal_onda_senoidal_volts);

% Cálculo analítico do valor médio teórico da senoide
valor_medio_teorico_senoidal_volts = offset_tensao_dc_senoidal_volts;

% Cálculo do valor eficaz RMS numérico medido via raiz do valor quadrático médio
valor_rms_medido_senoidal_volts = sqrt(mean(sinal_onda_senoidal_volts .^ 2));

% Cálculo analítico teórico do valor eficaz RMS da senoide com nível DC
valor_rms_teorico_senoidal_volts = sqrt(offset_tensao_dc_senoidal_volts^2 + (amplitude_sinal_senoidal_volts^2) / 2.0);

% Vetores de linhas de referência constantes para exibição no gráfico da senoide
vetor_linha_medio_senoidal = valor_medio_medido_senoidal_volts * ones(1, numero_pontos_amostras);
vetor_linha_rms_senoidal   = valor_rms_medido_senoidal_volts   * ones(1, numero_pontos_amostras);

% ==============================================================================
% PROCESSAMENTO NUMÉRICO 2: SÍNTESE E ANÁLISE DA ONDA QUADRADA
% ==============================================================================

% Geração analítica e universal da onda quadrada simétrica (50% de ciclo ativo)
% sign(sin(...)) gera +1 no semiciclo positivo e -1 no semiciclo negativo
sinal_onda_quadrada_volts = offset_tensao_dc_quadrado_volts + ...
    amplitude_sinal_quadrado_volts * sign(sin(2.0 * pi * frequencia_sinal_quadrado_hertz * vetor_tempo_segundos));

% Cálculo do valor médio numérico medido no conjunto de amostras da onda quadrada
valor_medio_medido_quadrado_volts = mean(sinal_onda_quadrada_volts);

% Cálculo analítico do valor médio teórico da onda quadrada
valor_medio_teorico_quadrado_volts = offset_tensao_dc_quadrado_volts;

% Cálculo do valor eficaz RMS numérico medido da onda quadrada
valor_rms_medido_quadrado_volts = sqrt(mean(sinal_onda_quadrada_volts .^ 2));

% Cálculo analítico teórico do valor eficaz RMS da onda quadrada com offset DC
valor_rms_teorico_quadrado_volts = sqrt(offset_tensao_dc_quadrado_volts^2 + amplitude_sinal_quadrado_volts^2);

% Vetores de linhas de referência para exibição gráfica
vetor_linha_medio_quadrado = valor_medio_medido_quadrado_volts * ones(1, numero_pontos_amostras);
vetor_linha_rms_quadrado   = valor_rms_medido_quadrado_volts   * ones(1, numero_pontos_amostras);

% ==============================================================================
% PROCESSAMENTO NUMÉRICO 3: SÍNTESE E ANÁLISE DA ONDA TRIANGULAR
% ==============================================================================

% Síntese universal da onda triangular simétrica variando entre -1 e +1
% Utiliza a função periódica linear modular: 2 * |2 * (f*t - floor(f*t + 0.5))| - 1
fase_normalizada_triangular = frequencia_sinal_triangular_hertz * vetor_tempo_segundos;
onda_triangular_normalizada = 2.0 * abs(2.0 * (fase_normalizada_triangular - floor(fase_normalizada_triangular + 0.5))) - 1.0;

% Aplicação da amplitude de pico e da componente contínua
sinal_onda_triangular_volts = offset_tensao_dc_triangular_volts + ...
    amplitude_sinal_triangular_volts * onda_triangular_normalizada;

% Cálculo do valor médio numérico medido da onda triangular
valor_medio_medido_triangular_volts = mean(sinal_onda_triangular_volts);

% Cálculo analítico do valor médio teórico da onda triangular
valor_medio_teorico_triangular_volts = offset_tensao_dc_triangular_volts;

% Cálculo do valor eficaz RMS numérico medido da onda triangular
valor_rms_medido_triangular_volts = sqrt(mean(sinal_onda_triangular_volts .^ 2));

% Cálculo analítico teórico do valor eficaz RMS da onda triangular com offset DC
valor_rms_teorico_triangular_volts = sqrt(offset_tensao_dc_triangular_volts^2 + (amplitude_sinal_triangular_volts^2) / 3.0);

% Vetores de linhas de referência para exibição gráfica
vetor_linha_medio_triangular = valor_medio_medido_triangular_volts * ones(1, numero_pontos_amostras);
vetor_linha_rms_triangular   = valor_rms_medido_triangular_volts   * ones(1, numero_pontos_amostras);

% ==============================================================================
% EXIBIÇÃO DE RESULTADOS QUANTITATIVOS NO CONSOLE
% ==============================================================================

fprintf('\n==============================================================================\n');
fprintf('  RELATÓRIO DE SÍNTESE DE SINAIS E AVALIAÇÃO DE VALOR MÉDIO E EFICAZ (RMS)    \n');
fprintf('==============================================================================\n');
fprintf('  Parâmetros Gerais de Amostragem:\n');
fprintf('    - Frequência Máxima do Sistema:        %10.2f Hz\n', frequencia_maxima_sinais_hertz);
fprintf('    - Limite Mínimo de Nyquist (2*f_max):  %10.2f Hz\n', frequencia_nyquist_minima_hertz);
fprintf('    - Taxa Efetiva de Amostragem Adotada:  %10.2f Hz\n', taxa_amostragem_efetiva_hertz);
fprintf('    - Número Total de Amostras Geradas:    %10d pontos\n', numero_pontos_amostras);
fprintf('------------------------------------------------------------------------------\n');
fprintf('  Tabela Comparativa: Teórico vs Medido Discreto\n');
fprintf('------------------------------------------------------------------------------\n');
fprintf('  Sinal        | Parâmetro   | Teórico (V) | Medido (V)  | Erro Relativo (%%)\n');
fprintf('------------------------------------------------------------------------------\n');

% Cálculo de erros relativos e absolutos com nomes de variáveis plenamente descritivos
erro_relativo_medio_senoidal_percentual = abs(valor_medio_medido_senoidal_volts - valor_medio_teorico_senoidal_volts) / max(abs(valor_medio_teorico_senoidal_volts), 1e-9) * 100.0;
erro_relativo_rms_senoidal_percentual   = abs(valor_rms_medido_senoidal_volts - valor_rms_teorico_senoidal_volts) / valor_rms_teorico_senoidal_volts * 100.0;

erro_relativo_medio_quadrado_percentual = abs(valor_medio_medido_quadrado_volts - valor_medio_teorico_quadrado_volts) / max(abs(valor_medio_teorico_quadrado_volts), 1e-9) * 100.0;
erro_relativo_rms_quadrado_percentual   = abs(valor_rms_medido_quadrado_volts - valor_rms_teorico_quadrado_volts) / valor_rms_teorico_quadrado_volts * 100.0;

erro_absoluto_medio_triangular_volts    = abs(valor_medio_medido_triangular_volts - valor_medio_teorico_triangular_volts);
erro_relativo_rms_triangular_percentual = abs(valor_rms_medido_triangular_volts - valor_rms_teorico_triangular_volts) / valor_rms_teorico_triangular_volts * 100.0;

fprintf('  Senoidal     | Valor Médio | %11.4f | %11.4f | %14.4f\n', valor_medio_teorico_senoidal_volts, valor_medio_medido_senoidal_volts, erro_relativo_medio_senoidal_percentual);
fprintf('  Senoidal     | Valor RMS   | %11.4f | %11.4f | %14.4f\n', valor_rms_teorico_senoidal_volts, valor_rms_medido_senoidal_volts, erro_relativo_rms_senoidal_percentual);
fprintf('  Quadrada     | Valor Médio | %11.4f | %11.4f | %14.4f\n', valor_medio_teorico_quadrado_volts, valor_medio_medido_quadrado_volts, erro_relativo_medio_quadrado_percentual);
fprintf('  Quadrada     | Valor RMS   | %11.4f | %11.4f | %14.4f\n', valor_rms_teorico_quadrado_volts, valor_rms_medido_quadrado_volts, erro_relativo_rms_quadrado_percentual);
fprintf('  Triangular   | Valor Médio | %11.4f | %11.4f | %14.4f (abs)\n', valor_medio_teorico_triangular_volts, valor_medio_medido_triangular_volts, erro_absoluto_medio_triangular_volts);
fprintf('  Triangular   | Valor RMS   | %11.4f | %11.4f | %14.4f\n', valor_rms_teorico_triangular_volts, valor_rms_medido_triangular_volts, erro_relativo_rms_triangular_percentual);
fprintf('==============================================================================\n\n');

% ==============================================================================
% VISUALIZAÇÃO GRÁFICA DOS RESULTADOS
% ==============================================================================

% Escala de tempo convertida para milissegundos para facilitar leitura humana
vetor_tempo_milissegundos = vetor_tempo_segundos * 1000.0;

% ------------------------------------------------------------------------------
% Figura 1: Forma de Onda Senoidal com Linhas de Referência
% ------------------------------------------------------------------------------
figura_senoide = figure(1);
set(figura_senoide, 'Name', 'Onda Senoidal - Analise de RMS e Media');
hold on;
plot(vetor_tempo_milissegundos, sinal_onda_senoidal_volts, 'r-', 'LineWidth', 2.0);
plot(vetor_tempo_milissegundos, vetor_linha_medio_senoidal, 'k--', 'LineWidth', 1.8);
plot(vetor_tempo_milissegundos, vetor_linha_rms_senoidal,   'b-.', 'LineWidth', 1.8);
plot(vetor_tempo_milissegundos, vetor_linha_zero_referencia, ':', 'Color', [0.4, 0.4, 0.4], 'LineWidth', 1.0);
grid on;
title('Sinal Senoidal com Offset Contínuo e Níveis Notáveis', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Tempo decorrido (milissegundos)', 'FontSize', 12);
ylabel('Tensão instantânea (Volts)', 'FontSize', 12);

texto_legenda_medio_senoide = sprintf('Valor Médio = %.4f V', valor_medio_medido_senoidal_volts);
texto_legenda_rms_senoide   = sprintf('Valor RMS = %.4f V', valor_rms_medido_senoidal_volts);
identificador_legenda_figura_1 = legend('Onda Senoidal', ...
                                        texto_legenda_medio_senoide, ...
                                        texto_legenda_rms_senoide, ...
                                        'Referência 0 V', ...
                                        'Location', 'northeast');
set(identificador_legenda_figura_1, 'FontSize', 11);
hold off;

% ------------------------------------------------------------------------------
% Figura 2: Forma de Onda Quadrada com Linhas de Referência
% ------------------------------------------------------------------------------
figura_quadrada = figure(2);
set(figura_quadrada, 'Name', 'Onda Quadrada - Analise de RMS e Media');
hold on;
plot(vetor_tempo_milissegundos, sinal_onda_quadrada_volts, 'm-', 'LineWidth', 2.0);
plot(vetor_tempo_milissegundos, vetor_linha_medio_quadrado, 'k--', 'LineWidth', 1.8);
plot(vetor_tempo_milissegundos, vetor_linha_rms_quadrado,   'b-.', 'LineWidth', 1.8);
plot(vetor_tempo_milissegundos, vetor_linha_zero_referencia, ':', 'Color', [0.4, 0.4, 0.4], 'LineWidth', 1.0);
grid on;
title('Sinal de Onda Quadrada Periódica com Níveis de Tensão', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Tempo decorrido (milissegundos)', 'FontSize', 12);
ylabel('Tensão instantânea (Volts)', 'FontSize', 12);

texto_legenda_medio_quadrada = sprintf('Valor Médio = %.4f V', valor_medio_medido_quadrado_volts);
texto_legenda_rms_quadrada   = sprintf('Valor RMS = %.4f V', valor_rms_medido_quadrado_volts);
identificador_legenda_figura_2 = legend('Onda Quadrada', ...
                                        texto_legenda_medio_quadrada, ...
                                        texto_legenda_rms_quadrada, ...
                                        'Referência 0 V', ...
                                        'Location', 'northeast');
set(identificador_legenda_figura_2, 'FontSize', 11);
hold off;

% ------------------------------------------------------------------------------
% Figura 3: Forma de Onda Triangular com Linhas de Referência
% ------------------------------------------------------------------------------
figura_triangular = figure(3);
set(figura_triangular, 'Name', 'Onda Triangular - Analise de RMS e Media');
hold on;
plot(vetor_tempo_milissegundos, sinal_onda_triangular_volts, 'b-', 'LineWidth', 2.0);
plot(vetor_tempo_milissegundos, vetor_linha_medio_triangular, 'k--', 'LineWidth', 1.8);
plot(vetor_tempo_milissegundos, vetor_linha_rms_triangular,   'r-.', 'LineWidth', 1.8);
plot(vetor_tempo_milissegundos, vetor_linha_zero_referencia, ':', 'Color', [0.4, 0.4, 0.4], 'LineWidth', 1.0);
grid on;
title('Sinal Triangular Simétrico com Níveis de Tensão', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Tempo decorrido (milissegundos)', 'FontSize', 12);
ylabel('Tensão instantânea (Volts)', 'FontSize', 12);

texto_legenda_medio_triangular = sprintf('Valor Médio = %.4f V', valor_medio_medido_triangular_volts);
texto_legenda_rms_triangular   = sprintf('Valor RMS = %.4f V', valor_rms_medido_triangular_volts);
identificador_legenda_figura_3 = legend('Onda Triangular', ...
                                        texto_legenda_medio_triangular, ...
                                        texto_legenda_rms_triangular, ...
                                        'Referência 0 V', ...
                                        'Location', 'northeast');
set(identificador_legenda_figura_3, 'FontSize', 11);
hold off;

% ------------------------------------------------------------------------------
% Figura 4: Comparação Simultânea das Três Ondas Sobrepostas
% (Correção do handle de legenda L4 que apontava incorretamente para L3 no legado)
% ------------------------------------------------------------------------------
figura_sobreposta = figure(4);
set(figura_sobreposta, 'Name', 'Comparacao Sobreposta de Sinais Periodicos');
hold on;
plot(vetor_tempo_milissegundos, sinal_onda_senoidal_volts,   'r-', 'LineWidth', 2.0);
plot(vetor_tempo_milissegundos, sinal_onda_quadrada_volts,   'g-', 'LineWidth', 2.0);
plot(vetor_tempo_milissegundos, sinal_onda_triangular_volts, 'b-', 'LineWidth', 2.0);
plot(vetor_tempo_milissegundos, vetor_linha_zero_referencia, ':', 'Color', [0.3, 0.3, 0.3], 'LineWidth', 1.2);
grid on;
title('Sobreposição Comparativa das Formas de Onda Sintetizadas', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Tempo decorrido (milissegundos)', 'FontSize', 12);
ylabel('Tensão instantânea (Volts)', 'FontSize', 12);

% Criação e configuração do identificador de legenda correto para a Figura 4
identificador_legenda_figura_4 = legend('Onda Senoidal (1 kHz)', ...
                                        'Onda Quadrada (1 kHz)', ...
                                        'Onda Triangular (1 kHz)', ...
                                        'Referência 0 V', ...
                                        'Location', 'northeast');
% Aplicação da formatação de fonte no handle correto L4
set(identificador_legenda_figura_4, 'FontSize', 12);
hold off;
