% =========================================================================
% UNIVERSIDADE DE BRASÍLIA / CÁLCULO NUMÉRICO / COMPUTAÇÃO CIENTÍFICA
% SCRIPT: rk2_individual.m
% OBJETIVO: Resolução de EDO Escalar via Runge-Kutta 2ª Ordem (Heun)
% AUTOR ORIGINAL: Davi Monteiro Alves
% REFATORAÇÃO: Engenharia Matemática e Especialista em Computação Científica
% =========================================================================
%
% 1. O QUE CALCULA:
%    Aproxima a solucao de um Problema de Valor Inicial (PVI) para uma
%    Equacao Diferencial Ordinaria (EDO) escalar de primeira ordem:
%              dy/dx = f(x, y),   com condicao inicial y(x0) = y0
%
%    Fundamentacao teorica:
%    O metodo de Runge-Kutta de 2a ordem (Heun) e um algoritmo preditor-corretor
%    baseado na expansao em Serie de Taylor ate termos de segunda ordem:
%          y(x + h) = y(x) + h * y'(x) + (h^2 / 2) * y''(x) + O(h^3)
%    Pela regra da cadeia, y''(x) = df/dx + (df/dy) * f.
%    O metodo calcula duas inclinacoes estrategicas no intervalo [x_i, x_i + h]:
%       1) Inclinacao no inicio do subintervalo (Euler):
%             k1 = f(x_i, y_i)
%       2) Estimativa preditora no final do subintervalo:
%             y_pred = y_i + h * k1
%       3) Inclinacao no final do subintervalo avaliada com o valor predito:
%             k2 = f(x_i + h, y_pred)
%       4) Atualizacao corretora ponderada (Regra Trapezoidal das inclinacoes):
%             y_{i+1} = y_i + (h / 2) * (k1 + k2)
%    O erro local de truncamento por passo e O(h^3), e o erro global e O(h^2).
%
% 2. O QUE ENTREGA:
%    - Vetor com os pontos da malha discretizada x.
%    - Vetor com as aproximacoes calculadas da solucao y(x).
%    - Vetor com os valores analiticos exatos (quando fornecidos para comparacao).
%    - Vetores de erro absoluto e erro relativo percentual a cada passo.
%    - Vetores com os coeficientes de inclinacao k1 e k2 gerados.
%    - Relatorio tabular impresso no terminal com alinhamento rigoroso.
%    - Grafico comparativo da solucao numerica versus solucao exata e evolucao do erro.
%
% 3. O QUE PRECISA TER (PARAMETROS DE ENTRADA):
%    - funcao_derivada_f: funcao anonima f(x, y) que expressa dy/dx.
%    - coordenada_inicial_x0: ponto inicial do dominio de integracao.
%    - coordenada_final_xf: ponto final do dominio de integracao.
%    - condicao_inicial_y0: valor inicial da solucao em x0, y(x0) = y0.
%    - tamanho_passo_h: largura uniforme de cada passo de integracao.
%
% 4. EXEMPLO DE USO E RESULTADO ESPERADO:
%    PVI: dy/dx = -2*y + 4*x,   com y(0) = 1, h = 0.1 no intervalo [0, 2]
%    Solucao analitica exata: y(x) = 2*x - 1 + 2*exp(-2*x)
%    Em x = 0.1:
%       k1 = f(0, 1) = -2(1) + 4(0) = -2.0000
%       y_pred = 1 + 0.1*(-2.0000) = 0.8000
%       k2 = f(0.1, 0.8) = -2(0.8) + 4(0.1) = -1.2000
%       y(0.1) = 1 + 0.05*(-2.0000 + -1.2000) = 0.840000 (Exato = 0.837462)
%       Erro absoluto = 0.002538 (Erro relativo = 0.303 %)
% =========================================================================

% -------------------------------------------------------------------------
% ETAPA 0: LIMPEZA DO AMBIENTE DE EXECUCAO
% -------------------------------------------------------------------------
% Limpa o terminal de saida
clc;
% Esvazia o espaco de variaveis da memoria
clear all;
% Fecha janelas graficas abertas
close all;

% -------------------------------------------------------------------------
% ETAPA 1: DEFINICAO DOS PARAMETROS DE ENTRADA DO PVI
% -------------------------------------------------------------------------
% Funcao da equacao diferencial ordinaria dy/dx = f(x, y)
funcao_derivada_f = @(coordenada_x, coordenada_y) (-2.0 * coordenada_y + 4.0 * coordenada_x);

% Solucao analitica exata para afericao de precisao e calculo do erro real
funcao_solucao_analitica_exata = @(coordenada_x) (2.0 * coordenada_x - 1.0 + 2.0 * exp(-2.0 * coordenada_x));

% Ponto inicial do intervalo de integracao temporal/espacial
coordenada_inicial_x0 = 0.0;

% Limite superior do intervalo de integracao
coordenada_final_xf = 2.0;

% Condicao inicial conhecida da solucao y(x0) = y0
condicao_inicial_y0 = 1.0;

% Tamanho constante do incremento na variavel independente (passo h)
tamanho_passo_h = 0.1;

% -------------------------------------------------------------------------
% ETAPA 2: VALIDACAO PREVIA DE CONSISTENCIA MATEMATICA
% -------------------------------------------------------------------------
% Verifica se a coordenada final e estritamente superior a coordenada inicial
if coordenada_final_xf <= coordenada_inicial_x0
    error('Erro de validacao: A coordenada final xf deve ser estritamente maior que a inicial x0.');
end

% Verifica se o passo de integracao h e estritamente positivo
if tamanho_passo_h <= 0.0
    error('Erro de validacao: O tamanho do passo h deve ser estritamente positivo (h > 0).');
end

% Calcula e valida o numero total de subintervalos
numero_total_passos_discretizacao = round((coordenada_final_xf - coordenada_inicial_x0) / tamanho_passo_h);
if numero_total_passos_discretizacao < 1
    error('Erro de validacao: O dominio e passo configurados geram menos de um passo de integracao.');
end

% -------------------------------------------------------------------------
% ETAPA 3: ALOCACAO DE ESTRUTURAS DE MEMORIA E INICIALIZACAO
% -------------------------------------------------------------------------
% Vetor coluna para armazenar os nos da malha independente x
vetor_malha_independente_x = zeros(numero_total_passos_discretizacao + 1, 1);

% Vetor coluna para armazenar os valores aproximados de y(x) via RK2
vetor_solucao_aproximada_y = zeros(numero_total_passos_discretizacao + 1, 1);

% Vetor coluna para armazenar a solucao exata de referencia
vetor_solucao_analitica_exata_y = zeros(numero_total_passos_discretizacao + 1, 1);

% Vetor coluna para armazenar o erro absoluto em cada no da discretizacao
vetor_erro_absoluto_pontual = zeros(numero_total_passos_discretizacao + 1, 1);

% Vetor coluna para armazenar o erro relativo percentual normalizado
vetor_erro_relativo_percentual = zeros(numero_total_passos_discretizacao + 1, 1);

% Vetores de inclinacao k1 e k2 para registro e documentacao numerica
vetor_inclinacao_inicial_k1 = zeros(numero_total_passos_discretizacao, 1);
vetor_inclinacao_preditora_k2 = zeros(numero_total_passos_discretizacao, 1);

% Define as condicoes iniciais no primeiro elemento dos vetores
vetor_malha_independente_x(1) = coordenada_inicial_x0;
vetor_solucao_aproximada_y(1) = condicao_inicial_y0;
vetor_solucao_analitica_exata_y(1) = funcao_solucao_analitica_exata(coordenada_inicial_x0);
vetor_erro_absoluto_pontual(1) = 0.0;
vetor_erro_relativo_percentual(1) = 0.0;

% -------------------------------------------------------------------------
% ETAPA 4: LACO DE INTEGRACAO NUMERICA VIA RUNGE-KUTTA DE SEGUNDA ORDEM (HEUN)
% -------------------------------------------------------------------------
% Executa a marcha temporal passo a passo
for indice_passo_integracao = 1:numero_total_passos_discretizacao

    % Ponto independente no inicio do subintervalo atual
    coordenada_x_corrente = vetor_malha_independente_x(indice_passo_integracao);

    % Estimativa da variavel dependente no inicio do subintervalo atual
    valor_y_corrente = vetor_solucao_aproximada_y(indice_passo_integracao);

    % 1. Calcula a primeira inclinacao k1 no ponto de partida (x_i, y_i)
    inclinacao_inicial_k1 = funcao_derivada_f(coordenada_x_corrente, valor_y_corrente);
    vetor_inclinacao_inicial_k1(indice_passo_integracao) = inclinacao_inicial_k1;

    % 2. Executa a estimativa preditora provisoria pelo metodo de Euler no ponto avancado
    coordenada_x_avancada = coordenada_inicial_x0 + (indice_passo_integracao * tamanho_passo_h);
    estimativa_preditora_y = valor_y_corrente + tamanho_passo_h * inclinacao_inicial_k1;

    % 3. Calcula a segunda inclinacao k2 no ponto predito (x_i + h, y_pred)
    inclinacao_preditora_k2 = funcao_derivada_f(coordenada_x_avancada, estimativa_preditora_y);
    vetor_inclinacao_preditora_k2(indice_passo_integracao) = inclinacao_preditora_k2;

    % 4. Calcula a media das duas inclinacoes (Regra Trapezoidal para integracao da derivada)
    inclinacao_media_corretora = (inclinacao_inicial_k1 + inclinacao_preditora_k2) / 2.0;

    % 5. Atualiza o estado da solucao no proximo no da malha
    vetor_malha_independente_x(indice_passo_integracao + 1) = coordenada_x_avancada;
    vetor_solucao_aproximada_y(indice_passo_integracao + 1) = valor_y_corrente + tamanho_passo_h * inclinacao_media_corretora;

    % 6. Avaliacao comparativa com a solucao exata
    valor_exato_avancado = funcao_solucao_analitica_exata(coordenada_x_avancada);
    vetor_solucao_analitica_exata_y(indice_passo_integracao + 1) = valor_exato_avancado;

    % Calcula o erro absoluto pontual
    erro_absoluto_passo = abs(valor_exato_avancado - vetor_solucao_aproximada_y(indice_passo_integracao + 1));
    vetor_erro_absoluto_pontual(indice_passo_integracao + 1) = erro_absoluto_passo;

    % Calcula o erro relativo percentual
    if abs(valor_exato_avancado) > 1e-12
        vetor_erro_relativo_percentual(indice_passo_integracao + 1) = (erro_absoluto_passo / abs(valor_exato_avancado)) * 100.0;
    else
        vetor_erro_relativo_percentual(indice_passo_integracao + 1) = erro_absoluto_passo * 100.0;
    end

% Conclui a iteracao do passo
end

% -------------------------------------------------------------------------
% ETAPA 5: EXIBICAO FORMATADA DOS DADOS NUMERICOS EM TABELA
% -------------------------------------------------------------------------
fprintf('========================================================================================================\n');
fprintf('     METODO DE RUNGE-KUTTA DE 2a ORDEM (HEUN / EULER APERFEICOADO) - EDO ESCALAR INDIVIDUAL            \n');
fprintf('========================================================================================================\n');
fprintf(' Parametros: x0 = %.4f | xf = %.4f | y0 = %.4f | Passo h = %.4f | Total de Passos = %d\n', ...
        coordenada_inicial_x0, coordenada_final_xf, condicao_inicial_y0, ...
        tamanho_passo_h, numero_total_passos_discretizacao);
fprintf('--------------------------------------------------------------------------------------------------------\n');
fprintf('  Passo  |      x      |    y (RK2)    |   y (Exato)   |     k1      |     k2      |  Erro Absoluto |  Erro (%%) \n');
fprintf('--------------------------------------------------------------------------------------------------------\n');

% Exibe o ponto inicial
fprintf('  %4d   |  %9.4f  |  %11.6f  |  %11.6f  |  ---------  |  ---------  |  %12.6e  |  %8.4f %%\n', ...
        0, vetor_malha_independente_x(1), vetor_solucao_aproximada_y(1), ...
        vetor_solucao_analitica_exata_y(1), vetor_erro_absoluto_pontual(1), vetor_erro_relativo_percentual(1));

% Exibe os passos calculados
for indice_exibicao_passo = 1:numero_total_passos_discretizacao
    fprintf('  %4d   |  %9.4f  |  %11.6f  |  %11.6f  |  %9.4f  |  %9.4f  |  %12.6e  |  %8.4f %%\n', ...
            indice_exibicao_passo, ...
            vetor_malha_independente_x(indice_exibicao_passo + 1), ...
            vetor_solucao_aproximada_y(indice_exibicao_passo + 1), ...
            vetor_solucao_analitica_exata_y(indice_exibicao_passo + 1), ...
            vetor_inclinacao_inicial_k1(indice_exibicao_passo), ...
            vetor_inclinacao_preditora_k2(indice_exibicao_passo), ...
            vetor_erro_absoluto_pontual(indice_exibicao_passo + 1), ...
            vetor_erro_relativo_percentual(indice_exibicao_passo + 1));
end
fprintf('========================================================================================================\n');

% -------------------------------------------------------------------------
% ETAPA 6: GERACAO DE FIGURAS GRAFICAS COMPARATIVAS
% -------------------------------------------------------------------------
figura_grafico_rk2 = figure('Name', 'Aproximacao via Runge-Kutta 2a Ordem (Heun)', 'NumberTitle', 'off');

% Subplot superior: Comparacao de trajetoria da solucao
subplot(2, 1, 1);
plot(vetor_malha_independente_x, vetor_solucao_analitica_exata_y, 'k-', 'LineWidth', 2.0);
hold on;
plot(vetor_malha_independente_x, vetor_solucao_aproximada_y, 'b--s', 'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerFaceColor', 'c');
grid on;
title('Solucao de EDO Escalar: Runge-Kutta 2a Ordem (Heun) vs. Solucao Exata', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Variavel Independente (x)', 'FontSize', 11);
ylabel('Variavel Dependente y(x)', 'FontSize', 11);
legend('Solucao Analitica Exata', 'Aproximacao RK2 (Heun)', 'Location', 'NorthEast');
hold off;

% Subplot inferior: Perfil de Erro Absoluto
subplot(2, 1, 2);
plot(vetor_malha_independente_x, vetor_erro_absoluto_pontual, 'r-^', 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', 'r');
grid on;
title('Comportamento do Erro Absoluto Global ao Longo da Integracao', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Variavel Independente (x)', 'FontSize', 11);
ylabel('Erro Absoluto |y_{exato} - y_{RK2}|', 'FontSize', 11);
legend('Erro Global O(h^2)', 'Location', 'NorthEast');

% =========================================================================
% FIM DO SCRIPT rk2_individual.m
% =========================================================================
