% =========================================================================
% UNIVERSIDADE DE BRASÍLIA / CÁLCULO NUMÉRICO / COMPUTAÇÃO CIENTÍFICA
% SCRIPT: euler_explicito.m
% OBJETIVO: Resolução Numérica de EDOs via Método de Euler Explícito
% AUTOR ORIGINAL: Davi Monteiro Alves
% REFATORAÇÃO: Engenharia Matemática e Especialista em Computação Científica
% =========================================================================
%
% 1. O QUE CALCULA:
%    Aproxima numericamente a solução de um Problema de Valor Inicial (PVI)
%    para uma Equação Diferencial Ordinária (EDO) de primeira ordem da forma:
%              dy/dx = f(x, y),   com condição inicial y(x0) = y0
%    A formulação matemática baseia-se na expansão em Série de Taylor de
%    primeira ordem truncada em torno do ponto atual x_i:
%              y(x_i + h) = y(x_i) + h * y'(x_i) + O(h^2)
%    Substituindo a derivada ordinária pela inclinação da reta tangente:
%              y_{i+1} = y_i + h * f(x_i, y_i)
%    O erro local de truncamento por passo é de ordem O(h^2), enquanto o
%    erro global acumulado ao longo de todo o intervalo é de ordem O(h).
%
% 2. O QUE ENTREGA:
%    - Vetor com os pontos da malha espacial/temporal discretizada (x).
%    - Vetor com as aproximações numéricas da função resposta y(x).
%    - Vetor com a solução analítica exata (quando disponível para comparação).
%    - Vetores com erro absoluto e erro relativo percentual a cada passo.
%    - Tabela comparativa detalhada impressa no console com formatação precisa.
%    - Gráfico comparativo entre a solução numérica de Euler e a solução exata.
%
% 3. O QUE PRECISA TER (PARÂMETROS DE ENTRADA):
%    - funcao_derivada_f: função anônima f(x, y) que expressa a derivada dy/dx.
%    - coordenada_inicial_x0: valor inicial da variável independente x.
%    - coordenada_final_xf: valor final da variável independente x (xf > x0).
%    - condicao_inicial_y0: valor inicial da variável dependente y em x0.
%    - tamanho_passo_h: espaçamento fixo e positivo entre pontos consecutivos da malha.
%
% 4. EXEMPLO DE USO E RESULTADO ESPERADO:
%    PVI: dy/dx = -2*x^3 + 12*x^2 - 20*x + 8.5,  com y(0) = 1 em [0, 4]
%    Solução analítica exata: y(x) = -0.5*x^4 + 4*x^3 - 10*x^2 + 8.5*x + 1
%    Exemplo com h = 0.5 (passo único didático):
%      Em x = 0.5: Euler calcula y(0.5) = 1 + 0.5*(8.5) = 5.2500 (Exato = 3.21875)
%    Configuração padrão adotada no script (h = 0.1, 40 passos):
%      Em x = 0.1: y_Euler = 1.850000, y_Exato = 1.753950, Erro Abs = 9.605e-02
%      Em x = 0.5: y_Euler = 3.590000, y_Exato = 3.218750, Erro Abs = 3.7125e-01
% =========================================================================

% -------------------------------------------------------------------------
% ETAPA 0: LIMPEZA DO AMBIENTE DE EXECUCAO
% -------------------------------------------------------------------------
% Limpa as mensagens da janela de comando do terminal
clc;
% Elimina todas as variaveis da memoria para evitar conflitos de escopo
clear all;
% Fecha todas as janelas de figuras graficas abertas anteriormente
close all;

% -------------------------------------------------------------------------
% ETAPA 1: DEFINICAO DOS PARAMETROS DE ENTRADA DO PVI
% -------------------------------------------------------------------------
% Funcao anonima f(x, y) que governa a taxa de variacao instantanea dy/dx
funcao_derivada_f = @(coordenada_x, coordenada_y) (-2.0 * (coordenada_x .^ 3) + 12.0 * (coordenada_x .^ 2) - 20.0 * coordenada_x + 8.5);

% Solucao analitica exata obtida por integracao direta para avaliacao de precisao
funcao_solucao_analitica_exata = @(coordenada_x) (-0.5 * (coordenada_x .^ 4) + 4.0 * (coordenada_x .^ 3) - 10.0 * (coordenada_x .^ 2) + 8.5 * coordenada_x + 1.0);

% Coordenada independente inicial (ponto inicial do intervalo de integracao)
coordenada_inicial_x0 = 0.0;

% Coordenada independente final (limite superior do intervalo de integracao)
coordenada_final_xf = 4.0;

% Condicao inicial conhecida da funcao dependente no ponto inicial y(x0) = y0
condicao_inicial_y0 = 1.0;

% Tamanho do passo de discretizacao da malha (passo de integracao h)
tamanho_passo_h = 0.1;

% -------------------------------------------------------------------------
% ETAPA 2: VALIDACAO PREVIA DE CONSISTENCIA MATEMATICA
% -------------------------------------------------------------------------
% Verifica se o limite final da variavel independente e estritamente maior que o inicial
if coordenada_final_xf <= coordenada_inicial_x0
    error('Erro de validacao: A coordenada final xf deve ser estritamente maior que a inicial x0.');
end

% Verifica se o tamanho do passo h e estritamente positivo
if tamanho_passo_h <= 0.0
    error('Erro de validacao: O tamanho do passo h deve ser estritamente positivo (h > 0).');
end

% Calcula e valida o numero total de subintervalos na discretizacao
numero_total_passos_discretizacao = round((coordenada_final_xf - coordenada_inicial_x0) / tamanho_passo_h);
if numero_total_passos_discretizacao < 1
    error('Erro de validacao: O dominio e passo configurados geram menos de um passo de integracao.');
end

% -------------------------------------------------------------------------
% ETAPA 3: PRE-ALOCACAO DE MEMORIA E INICIALIZACAO DOS VETORES
% -------------------------------------------------------------------------
% Pre-aloca o vetor coluna que contera todos os nos discretos da malha x
vetor_malha_espacial_x = zeros(numero_total_passos_discretizacao + 1, 1);

% Pre-aloca o vetor coluna que contera a solucao numerica aproximada y(x)
vetor_solucao_aproximada_y = zeros(numero_total_passos_discretizacao + 1, 1);

% Pre-aloca o vetor coluna que contera os valores da solucao analitica exata
vetor_solucao_analitica_exata_y = zeros(numero_total_passos_discretizacao + 1, 1);

% Pre-aloca o vetor coluna que armazenara os erros absolutos pontuais
vetor_erro_absoluto_pontual = zeros(numero_total_passos_discretizacao + 1, 1);

% Pre-aloca o vetor coluna que armazenara os erros relativos percentuais
vetor_erro_relativo_percentual = zeros(numero_total_passos_discretizacao + 1, 1);

% Atribui o ponto inicial da malha na primeira posicao do vetor de coordenadas
vetor_malha_espacial_x(1) = coordenada_inicial_x0;

% Atribui a condicao inicial fornecida na primeira posicao da solucao numerica
vetor_solucao_aproximada_y(1) = condicao_inicial_y0;

% Avalia a solucao exata analitica no ponto inicial da malha
vetor_solucao_analitica_exata_y(1) = funcao_solucao_analitica_exata(coordenada_inicial_x0);

% No ponto inicial, o erro absoluto e nulo pois o valor inicial e exato
vetor_erro_absoluto_pontual(1) = abs(vetor_solucao_analitica_exata_y(1) - vetor_solucao_aproximada_y(1));

% No ponto inicial, o erro relativo percentual e identicamente zero
vetor_erro_relativo_percentual(1) = 0.0;

% -------------------------------------------------------------------------
% ETAPA 4: INTEGRACAO NUMERICA PELO METODO DE EULER EXPLICITO
% -------------------------------------------------------------------------
% Itera sequencialmente sobre cada subintervalo da malha temporal/espacial
for indice_passo_integracao = 1:numero_total_passos_discretizacao

    % Recupera a coordenada independente x do ponto corrente da iteracao
    coordenada_x_corrente = vetor_malha_espacial_x(indice_passo_integracao);

    % Recupera o valor aproximado da variavel dependente y no ponto corrente
    coordenada_y_corrente = vetor_solucao_aproximada_y(indice_passo_integracao);

    % Calcula a inclinacao da reta tangente local avaliando a EDO f(x_i, y_i)
    inclinacao_tangente_euler = funcao_derivada_f(coordenada_x_corrente, coordenada_y_corrente);

    % Avanca a coordenada x para o proximo ponto da malha sem acumulo de erro de arredondamento
    vetor_malha_espacial_x(indice_passo_integracao + 1) = coordenada_inicial_x0 + (indice_passo_integracao * tamanho_passo_h);

    % Aplica a formula de extrapolacao linear de Euler para estimar y_{i+1}
    % y_{i+1} = y_i + h * f(x_i, y_i)
    vetor_solucao_aproximada_y(indice_passo_integracao + 1) = coordenada_y_corrente + tamanho_passo_h * inclinacao_tangente_euler;

    % Recupera a nova coordenada x recem-calculada para avaliacoes comparativas
    coordenada_x_seguinte = vetor_malha_espacial_x(indice_passo_integracao + 1);

    % Calcula o valor teorico exato da solucao no novo ponto x_{i+1}
    valor_solucao_exata_seguinte = funcao_solucao_analitica_exata(coordenada_x_seguinte);

    % Armazena a solucao analitica correspondente no vetor comparativo
    vetor_solucao_analitica_exata_y(indice_passo_integracao + 1) = valor_solucao_exata_seguinte;

    % Calcula o erro absoluto pontual entre a solucao exata e a aproximacao numerica
    erro_absoluto_atual = abs(valor_solucao_exata_seguinte - vetor_solucao_aproximada_y(indice_passo_integracao + 1));
    vetor_erro_absoluto_pontual(indice_passo_integracao + 1) = erro_absoluto_atual;

    % Calcula o erro relativo percentual normalizado pelo valor teorico exato
    if abs(valor_solucao_exata_seguinte) > 1e-12
        vetor_erro_relativo_percentual(indice_passo_integracao + 1) = (erro_absoluto_atual / abs(valor_solucao_exata_seguinte)) * 100.0;
    else
        % Evita divisao por zero caso a solucao analitica cruze o eixo horizontal
        vetor_erro_relativo_percentual(indice_passo_integracao + 1) = erro_absoluto_atual * 100.0;
    end

% Conclui o laco iterativo do metodo de passo simples
end

% -------------------------------------------------------------------------
% ETAPA 5: EXIBICAO DOS RESULTADOS TABELADOS NO TERMINAL
% -------------------------------------------------------------------------
% Imprime cabecalho informativo com descricao do metodo e parâmetros
fprintf('========================================================================================\n');
fprintf('         METODO DE EULER EXPLICITO PARA EQUACOES DIFERENCIAIS ORDINARIAS (PVI)         \n');
fprintf('========================================================================================\n');
fprintf(' Parametros: x0 = %.4f | xf = %.4f | y0 = %.4f | Passo h = %.4f | Passos = %d\n', ...
        coordenada_inicial_x0, coordenada_final_xf, condicao_inicial_y0, ...
        tamanho_passo_h, numero_total_passos_discretizacao);
fprintf('----------------------------------------------------------------------------------------\n');
fprintf('  Passo  |      x      |   y (Euler)   |   y (Exato)   |  Erro Absoluto | Erro Relativo (%%)\n');
fprintf('----------------------------------------------------------------------------------------\n');

% Itera imprimindo cada linha da tabela com formatacao numerica estrita
for indice_exibicao_tabela = 1:(numero_total_passos_discretizacao + 1)
    fprintf('  %4d   |  %9.4f  |  %11.6f  |  %11.6f  |  %12.6e  |  %14.4f %%\n', ...
            indice_exibicao_tabela - 1, ...
            vetor_malha_espacial_x(indice_exibicao_tabela), ...
            vetor_solucao_aproximada_y(indice_exibicao_tabela), ...
            vetor_solucao_analitica_exata_y(indice_exibicao_tabela), ...
            vetor_erro_absoluto_pontual(indice_exibicao_tabela), ...
            vetor_erro_relativo_percentual(indice_exibicao_tabela));
end
fprintf('========================================================================================\n');

% -------------------------------------------------------------------------
% ETAPA 6: VISUALIZACAO GRAFICA COMPARATIVA DOS RESULTADOS
% -------------------------------------------------------------------------
% Cria uma nova figura para exibicao visual das curvas
figura_comparativa_solucoes = figure('Name', 'Solucao de EDO via Metodo de Euler Explicito', 'NumberTitle', 'off');

% Subplot superior: comparacao direta entre curvas numerica e analitica
subplot(2, 1, 1);
% Plota a solucao exata analitica como linha continua preta
plot(vetor_malha_espacial_x, vetor_solucao_analitica_exata_y, 'k-', 'LineWidth', 2.0);
hold on;
% Plota a solucao aproximada do metodo de Euler com marcadores circulares vermelhos
plot(vetor_malha_espacial_x, vetor_solucao_aproximada_y, 'ro--', 'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerFaceColor', 'r');
grid on;
title('Solucao Numerica de Euler Explicito vs. Solucao Analitica Exata', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Variavel Independente (x)', 'FontSize', 11);
ylabel('Variavel Dependente y(x)', 'FontSize', 11);
legend('Solucao Analitica Exata', 'Aproximacao de Euler Explicito', 'Location', 'NorthWest');
hold off;

% Subplot inferior: evolucao do erro absoluto ao longo do dominio
subplot(2, 1, 2);
% Plota a evolucao do erro absoluto em escala convencional
plot(vetor_malha_espacial_x, vetor_erro_absoluto_pontual, 'b-s', 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', 'b');
grid on;
title('Evolucao do Erro Absoluto de Truncamento Global |y_{exato} - y_{euler}|', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Variavel Independente (x)', 'FontSize', 11);
ylabel('Erro Absoluto Pontual', 'FontSize', 11);
legend('Erro Absoluto |E(x)|', 'Location', 'NorthWest');

% =========================================================================
% FIM DO SCRIPT euler_explicito.m
% =========================================================================
