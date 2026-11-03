% =========================================================================
% UNIVERSIDADE DE BRASÍLIA / CÁLCULO NUMÉRICO / COMPUTAÇÃO CIENTÍFICA
% SCRIPT: rk2_heun_sistema.m
% OBJETIVO: Resolução de Sistema de EDOs via Runge-Kutta 2ª Ordem (Heun)
% AUTOR ORIGINAL: Davi Monteiro Alves
% REFATORAÇÃO: Engenharia Matemática e Especialista em Computação Científica
% =========================================================================
%
% 1. O QUE CALCULA:
%    Resolve numericamente um sistema acoplado de 2 Equacoes Diferenciais
%    Ordinarias (EDOs) de primeira ordem com Problema de Valor Inicial (PVI):
%              dy/dx = f(x, y, z),   com condicao inicial y(x0) = y0
%              dz/dx = g(x, y, z),   com condicao inicial z(x0) = z0
%
%    A formulacao utiliza o algoritmo de Heun (Runge-Kutta de 2a ordem):
%    - ETAPA PREDITORA (Euler Simples nos pontos iniciais):
%          k_{y1} = f(x_i, y_i, z_i)
%          k_{z1} = g(x_i, y_i, z_i)
%          y_{pred} = y_i + h * k_{y1}
%          z_{pred} = z_i + h * k_{z1}
%    - ETAPA CORRETORA (Regra Trapezoidal usando as inclinacoes no ponto predito):
%          k_{y2} = f(x_i + h, y_{pred}, z_{pred})
%          k_{z2} = g(x_i + h, y_{pred}, z_{pred})
%          y_{i+1} = y_i + (h / 2) * (k_{y1} + k_{y2})
%          z_{i+1} = z_i + (h / 2) * (k_{z1} + k_{z2})
%    O erro de truncamento local e de ordem O(h^3) e o erro global acumulado e O(h^2).
%
% 2. O QUE ENTREGA:
%    - Vetor da malha independente x discretizada.
%    - Vetores com as trajetorias das solucoes acopladas y(x) e z(x).
%    - Tabela completa de passos com coordenadas e coeficientes k1 e k2.
%    - Grafico com as curvas de evolucao simultanea de y(x) e z(x).
%    - Grafico do espaco de fase (plano de fase y vs. z).
%
% 3. O QUE PRECISA TER (PARAMETROS DE ENTRADA):
%    - funcao_derivada_y_f: funcao anonima f(x, y, z) para a taxa dy/dx.
%    - funcao_derivada_z_g: funcao anonima g(x, y, z) para a taxa dz/dx.
%    - coordenada_inicial_x0: valor inicial da variavel independente.
%    - coordenada_final_xn: valor final da variavel independente.
%    - condicao_inicial_y0: valor inicial da variavel y em x0.
%    - condicao_inicial_z0: valor inicial da variavel z em x0.
%    - tamanho_passo_h: incremento uniforme entre passos sucessivos da malha.
%
% 4. EXEMPLO DE USO E RESULTADO ESPERADO:
%    Sistema: dy/dx = z*x - y,   dz/dx = y*x + z
%    Condicoes: x0 = 0.0, y(0) = 0.5, z(0) = 1.0, xn = 1.2, h = 0.1
%    Em x = 0.1:
%       k_{y1} = -0.5000, k_{z1} = 1.0000
%       y_{pred} = 0.4500, z_{pred} = 1.1000
%       k_{y2} = (1.1000*0.1) - 0.4500 = -0.3400
%       k_{z2} = (0.4500*0.1) + 1.1000 = 1.1450
%       y(0.1) = 0.5 + 0.05*(-0.5000 - 0.3400) = 0.4580
%       z(0.1) = 1.0 + 0.05*(1.0000 + 1.1450) = 1.10725
% =========================================================================

% -------------------------------------------------------------------------
% ETAPA 0: LIMPEZA DO AMBIENTE DE EXECUCAO
% -------------------------------------------------------------------------
% Limpa as mensagens residuais na janela de comando
clc;
% Elimina variaveis preexistentes no espaco de trabalho
clear all;
% Fecha todas as figuras abertas
close all;

% -------------------------------------------------------------------------
% ETAPA 1: DEFINICAO DOS PARAMETROS DO SISTEMA DE EDOS ACOPLADAS
% -------------------------------------------------------------------------
% Equacao diferencial para a primeira variavel dependente: dy/dx = z*x - y
funcao_derivada_y_f = @(coordenada_x, variavel_y, variavel_z) (variavel_z .* coordenada_x - variavel_y);

% Equacao diferencial para a segunda variavel dependente: dz/dx = y*x + z
funcao_derivada_z_g = @(coordenada_x, variavel_y, variavel_z) (variavel_y .* coordenada_x + variavel_z);

% Ponto inicial da variavel independente (tempo ou posicao)
coordenada_inicial_x0 = 0.0;

% Limite superior da integracao (ponto final do dominio)
coordenada_final_xn = 1.2;

% Condicao inicial da primeira funcao dependente y(x0) = y0
condicao_inicial_y0 = 0.5;

% Condicao inicial da segunda funcao dependente z(x0) = z0
condicao_inicial_z0 = 1.0;

% Incremento uniforme da malha temporal/espacial (tamanho do passo h)
tamanho_passo_h = 0.1;

% -------------------------------------------------------------------------
% ETAPA 2: VALIDACAO PREVIA DE CONSISTENCIA MATEMATICA
% -------------------------------------------------------------------------
% Verifica se o limite final e estritamente superior ao ponto inicial
if coordenada_final_xn <= coordenada_inicial_x0
    error('Erro de validacao: A coordenada final xn deve ser estritamente maior que a inicial x0.');
end

% Verifica se o incremento de passo h e estritamente positivo
if tamanho_passo_h <= 0.0
    error('Erro de validacao: O tamanho do passo h deve ser estritamente positivo (h > 0).');
end

% Valida se a malha gera ao menos um passo de discretizacao
numero_total_passos_discretizacao = round((coordenada_final_xn - coordenada_inicial_x0) / tamanho_passo_h);
if numero_total_passos_discretizacao < 1
    error('Erro de validacao: O dominio e passo configurados geram menos de um passo de integracao.');
end

% -------------------------------------------------------------------------
% ETAPA 3: ALOCACAO DE ESTRUTURAS DE DADOS E INICIALIZACAO
% -------------------------------------------------------------------------
% Pre-alocacao de memoria para a coordenada independente x
vetor_malha_independente_x = zeros(numero_total_passos_discretizacao + 1, 1);

% Pre-alocacao de memoria para os valores da solucao y(x)
vetor_solucao_variavel_y = zeros(numero_total_passos_discretizacao + 1, 1);

% Pre-alocacao de memoria para os valores da solucao z(x)
vetor_solucao_variavel_z = zeros(numero_total_passos_discretizacao + 1, 1);

% Pre-alocacao dos vetores de coeficientes de inclinacao para fins de relatorio
vetor_inclinacao_inicial_y_k1 = zeros(numero_total_passos_discretizacao, 1);
vetor_inclinacao_inicial_z_k1 = zeros(numero_total_passos_discretizacao, 1);
vetor_inclinacao_preditora_y_k2 = zeros(numero_total_passos_discretizacao, 1);
vetor_inclinacao_preditora_z_k2 = zeros(numero_total_passos_discretizacao, 1);

% Inicializa o primeiro no da malha com a coordenada inicial
vetor_malha_independente_x(1) = coordenada_inicial_x0;

% Atribui o valor da condicao inicial y(x0)
vetor_solucao_variavel_y(1) = condicao_inicial_y0;

% Atribui o valor da condicao inicial z(x0)
vetor_solucao_variavel_z(1) = condicao_inicial_z0;

% -------------------------------------------------------------------------
% ETAPA 4: LACO DE INTEGRACAO NUMERICA - RUNGE-KUTTA DE 2a ORDEM (HEUN)
% -------------------------------------------------------------------------
% Itera para avancar a solucao acoplada do passo i para o passo i+1
for indice_passo_integracao = 1:numero_total_passos_discretizacao

    % Variavel independente atual x_i
    coordenada_x_atual = vetor_malha_independente_x(indice_passo_integracao);

    % Variavel dependente y_i no inicio do subintervalo
    valor_y_atual = vetor_solucao_variavel_y(indice_passo_integracao);

    % Variavel dependente z_i no inicio do subintervalo
    valor_z_atual = vetor_solucao_variavel_z(indice_passo_integracao);

    % 1. Calcula as inclinacoes iniciais no ponto (x_i, y_i, z_i)
    inclinacao_inicial_y_k1 = funcao_derivada_y_f(coordenada_x_atual, valor_y_atual, valor_z_atual);
    inclinacao_inicial_z_k1 = funcao_derivada_z_g(coordenada_x_atual, valor_y_atual, valor_z_atual);

    % Registra as inclinacoes k1 calculadas
    vetor_inclinacao_inicial_y_k1(indice_passo_integracao) = inclinacao_inicial_y_k1;
    vetor_inclinacao_inicial_z_k1(indice_passo_integracao) = inclinacao_inicial_z_k1;

    % 2. Passo preditor de Euler: estima a solucao no proximo ponto x_{i+1}
    coordenada_x_avancada = coordenada_inicial_x0 + (indice_passo_integracao * tamanho_passo_h);
    estimativa_preditora_y = valor_y_atual + tamanho_passo_h * inclinacao_inicial_y_k1;
    estimativa_preditora_z = valor_z_atual + tamanho_passo_h * inclinacao_inicial_z_k1;

    % 3. Calcula as inclinacoes no ponto predito (x_{i+1}, y_pred, z_pred)
    inclinacao_preditora_y_k2 = funcao_derivada_y_f(coordenada_x_avancada, estimativa_preditora_y, estimativa_preditora_z);
    inclinacao_preditora_z_k2 = funcao_derivada_z_g(coordenada_x_avancada, estimativa_preditora_y, estimativa_preditora_z);

    % Registra as inclinacoes k2 calculadas
    vetor_inclinacao_preditora_y_k2(indice_passo_integracao) = inclinacao_preditora_y_k2;
    vetor_inclinacao_preditora_z_k2(indice_passo_integracao) = inclinacao_preditora_z_k2;

    % 4. Passo corretor: media aritmetica das inclinacoes (Regra Trapezoidal)
    inclinacao_media_ponderada_y = (inclinacao_inicial_y_k1 + inclinacao_preditora_y_k2) / 2.0;
    inclinacao_media_ponderada_z = (inclinacao_inicial_z_k1 + inclinacao_preditora_z_k2) / 2.0;

    % Atualiza a coordenada independente para o novo ponto
    vetor_malha_independente_x(indice_passo_integracao + 1) = coordenada_x_avancada;

    % Atualiza o estado da variavel y para o instante i+1
    vetor_solucao_variavel_y(indice_passo_integracao + 1) = valor_y_atual + tamanho_passo_h * inclinacao_media_ponderada_y;

    % Atualiza o estado da variavel z para o instante i+1
    vetor_solucao_variavel_z(indice_passo_integracao + 1) = valor_z_atual + tamanho_passo_h * inclinacao_media_ponderada_z;

% Fim do laco de integracao do sistema
end

% -------------------------------------------------------------------------
% ETAPA 5: APRESENTACAO DOS DADOS COMPUTADOS EM TABELA
% -------------------------------------------------------------------------
fprintf('====================================================================================================\n');
fprintf('     SISTEMA ACOPLADO DE EDOS - METODO DE RUNGE-KUTTA DE SEGUNDA ORDEM (HEUN PREDITOR-CORRETOR)    \n');
fprintf('====================================================================================================\n');
fprintf(' Parametros: x0 = %.4f | xn = %.4f | y0 = %.4f | z0 = %.4f | Passo h = %.4f | Passos = %d\n', ...
        coordenada_inicial_x0, coordenada_final_xn, condicao_inicial_y0, condicao_inicial_z0, ...
        tamanho_passo_h, numero_total_passos_discretizacao);
fprintf('----------------------------------------------------------------------------------------------------\n');
fprintf('  Passo  |      x      |    y(x)    |    z(x)    |   k_{y1}   |   k_{y2}   |   k_{z1}   |   k_{z2}   \n');
fprintf('----------------------------------------------------------------------------------------------------\n');

% Exibe o ponto inicial (sem inclinacoes anteriores)
fprintf('  %4d   |  %9.4f  |  %8.5f  |  %8.5f  |  --------  |  --------  |  --------  |  --------  \n', ...
        0, vetor_malha_independente_x(1), vetor_solucao_variavel_y(1), vetor_solucao_variavel_z(1));

% Exibe as demais linhas da tabela com os coeficientes calculados
for indice_exibicao_passo = 1:numero_total_passos_discretizacao
    fprintf('  %4d   |  %9.4f  |  %8.5f  |  %8.5f  |  %8.4f  |  %8.4f  |  %8.4f  |  %8.4f  \n', ...
            indice_exibicao_passo, ...
            vetor_malha_independente_x(indice_exibicao_passo + 1), ...
            vetor_solucao_variavel_y(indice_exibicao_passo + 1), ...
            vetor_solucao_variavel_z(indice_exibicao_passo + 1), ...
            vetor_inclinacao_inicial_y_k1(indice_exibicao_passo), ...
            vetor_inclinacao_preditora_y_k2(indice_exibicao_passo), ...
            vetor_inclinacao_inicial_z_k1(indice_exibicao_passo), ...
            vetor_inclinacao_preditora_z_k2(indice_exibicao_passo));
end
fprintf('====================================================================================================\n');

% -------------------------------------------------------------------------
% ETAPA 6: RENDERIZACAO GRAFICA DAS TRAJETORIAS TEMPORAIS E DO ESPACO DE FASE
% -------------------------------------------------------------------------
figura_sistema_edos = figure('Name', 'Solucao de Sistema de EDOs via RK2 (Heun)', 'NumberTitle', 'off');

% Subplot 1: Evolucao temporal/espacial das duas variaveis dependentes y(x) e z(x)
subplot(2, 1, 1);
plot(vetor_malha_independente_x, vetor_solucao_variavel_y, 'b-o', 'LineWidth', 1.8, 'MarkerSize', 5, 'MarkerFaceColor', 'b');
hold on;
plot(vetor_malha_independente_x, vetor_solucao_variavel_z, 'r-s', 'LineWidth', 1.8, 'MarkerSize', 5, 'MarkerFaceColor', 'r');
grid on;
title('Evolucao das Solucoes Acopladas y(x) e z(x) via Metodo de Heun (RK2)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Variavel Independente (x)', 'FontSize', 11);
ylabel('Variaveis de Estado (y, z)', 'FontSize', 11);
legend('Variavel y(x)', 'Variavel z(x)', 'Location', 'NorthWest');
hold off;

% Subplot 2: Retrato de fase (Trajetoria de y em funcao de z)
subplot(2, 1, 2);
plot(vetor_solucao_variavel_y, vetor_solucao_variavel_z, 'm-d', 'LineWidth', 1.8, 'MarkerSize', 5, 'MarkerFaceColor', 'm');
grid on;
title('Plano de Fase: Trajetoria de z(x) em relacao a y(x)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Variavel de Estado y(x)', 'FontSize', 11);
ylabel('Variavel de Estado z(x)', 'FontSize', 11);
legend('Trajetoria no Espaco de Estados', 'Location', 'NorthWest');

% =========================================================================
% FIM DO SCRIPT rk2_heun_sistema.m
% =========================================================================
