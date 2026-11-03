% =========================================================================
% UNIVERSIDADE DE BRASÍLIA / CÁLCULO NUMÉRICO / COMPUTAÇÃO CIENTÍFICA
% SCRIPT: diferenca_progressiva_3pontos.m
% OBJETIVO: Diferenciação numérica pela Fórmula Progressiva de Três Pontos O(h^2)
% AUTOR ORIGINAL: Davi Monteiro Alves
% REFATORAÇÃO: Engenharia Matemática e Especialista em Computação Científica
% =========================================================================
%
% -------------------------------------------------------------------------
% 1. O QUE CALCULA:
% -------------------------------------------------------------------------
%   Aproximação numérica da primeira derivada f'(x) de uma função contínua
%   ou conjunto de dados discretos utilizando a Fórmula de Diferença Finita
%   Progressiva de Três Pontos (Forward Finite Difference de 2ª ordem).
%
%   Formulação Matemática e Dedução via Séries de Taylor:
%     Expandindo f(x + h) e f(x + 2h) em torno do ponto x_0:
%       (1) f(x_0 + h)  = f(x_0) + h*f'(x_0) + (h^2/2)*f''(x_0) + (h^3/6)*f'''(x_0) + O(h^4)
%       (2) f(x_0 + 2h) = f(x_0) + 2h*f'(x_0) + 2h^2*f''(x_0) + (4h^3/3)*f'''(x_0) + O(h^4)
%
%     Multiplicando a equação (1) por 4 e subtraindo a equação (2), cancela-se
%     o termo de segunda derivada f''(x_0):
%       4*f(x_0 + h) - f(x_0 + 2h) = 3*f(x_0) + 2h*f'(x_0) - (2h^3/3)*f'''(xi)
%
%     Isolando f'(x_0), obtém-se a expressão da derivada:
%       f'(x_0) = \frac{-3*f(x_0) + 4*f(x_0 + h) - f(x_0 + 2h)}{2*h} + \frac{h^2}{3}*f'''(xi)
%
%     CORREÇÃO CRÍTICA DE PRECEDÊNCIA:
%       No código legado original constava:
%         dx(i) = (-3*u(i) + 4*u(i+1) - u(i+2)) / 2 * h;  <-- ERRO GRAVE!
%       Devido à associatividade da esquerda para a direita, a expressão acima
%       equivale matematicamente a:
%         [ (-3*u(i) + 4*u(i+1) - u(i+2)) / 2 ] * h
%       Isto multiplicava o resultado por h em vez de dividir por (2*h), gerando
%       um desvio catastrófico de fator h^2.
%       A expressão corrigida exige parênteses obrigatórios no denominador:
%         (-3*y_0 + 4*y_1 - y_2) / (2 * tamanho_passo_h)
%
%     Termo do Erro Global de Truncamento:
%       E_{prog3} = + \frac{h^2}{3} * f'''(xi),  com xi \in (x_0, x_0 + 2h).
%       Ordem de convergência: O(h^2) (segunda ordem).
%       Aplicação: Essencial em nós de fronteira esquerda de malhas discretas,
%       onde não existem nós anteriores para aplicar diferenças centrais.
%
% -------------------------------------------------------------------------
% 2. O QUE ENTREGA:
% -------------------------------------------------------------------------
%   - vetor_derivadas_numericas_dx: Vetor contendo a estimativa numérica da
%     derivada calculada em cada ponto x_k para k = 1, ..., N-2.
%   - vetor_erros_absolutos: Vetor com o erro absoluto ponto a ponto em relação
%     à derivada analítica exata de referência.
%   - Tabela comparativa impressa no console com amostragem dos primeiros nós.
%   - Painel gráfico bidimensional com dois subplots:
%       Subplot 1: Comparação visual entre derivada analítica e numérica.
%       Subplot 2: Perfil do erro absoluto evidenciando precisão O(h^2).
%
% -------------------------------------------------------------------------
% 3. O QUE PRECISA TER:
% -------------------------------------------------------------------------
%   - tamanho_passo_h: Escalar real estritamente positivo (h > 0).
%   - vetor_amostras_funcao_y: Vetor de dados amostrais com no mínimo 3 pontos.
%   - (Opcional demonstrativo) funcao_teste_f e sua derivada analítica
%     funcao_derivada_analitica_df para validação automática da precisão.
%
% -------------------------------------------------------------------------
% 4. EXEMPLO DE USO E RESULTADO ESPERADO:
% -------------------------------------------------------------------------
%   Função teste: f(x) = sin(2*x) => f'(x) = 2*cos(2*x)
%   Intervalo de amostragem: [0, pi] com h = 0.05
%   Resultado esperado:
%     Para x = 0.00: Derivada analítica = 2.000000, Derivada numérica = 2.006643
%     Erro absoluto máximo no domínio: ~ 6.66e-03 (conforme O(h^2) = 0.05^2/3 * 8).
% =========================================================================

% -------------------------------------------------------------------------
% ETAPA 0: LIMPEZA DO AMBIENTE DE TRABALHO
% -------------------------------------------------------------------------
% Limpa a tela do console de comandos
clc;
% Limpa todo o espaço de trabalho da memória
clear all;
% Fecha janelas de figuras anteriormente abertas
close all;

% -------------------------------------------------------------------------
% ETAPA 1: DEFINIÇÃO DOS PARÂMETROS E GERAÇÃO DOS DADOS DEMONSTRATIVOS
% -------------------------------------------------------------------------
% Espaçamento uniforme entre nós consecutivos da malha espacial (h)
tamanho_passo_h = 0.05;

% Limites inferior e superior para a malha demonstrativa de teste
limite_inferior_dominio = 0.0;
limite_superior_dominio = pi;

% Definição da função analítica de teste: f(x) = sin(2*x)
funcao_teste_f = @(coordenada_x) sin(2.0 .* coordenada_x);

% Definição da derivada analítica exata correspondente: f'(x) = 2*cos(2*x)
funcao_derivada_analitica_df = @(coordenada_x) 2.0 .* cos(2.0 .* coordenada_x);

% Construção do vetor com a malha discreta de coordenadas independentes x
vetor_coordenadas_x = limite_inferior_dominio : tamanho_passo_h : limite_superior_dominio;

% Avaliação das amostras dependentes y: y = f(x)
vetor_amostras_funcao_y = funcao_teste_f(vetor_coordenadas_x);

% -------------------------------------------------------------------------
% ETAPA 2: VALIDAÇÃO PRÉVIA DE CONSISTÊNCIA NUMÉRICA
% -------------------------------------------------------------------------
% Quantidade total de pontos amostrais disponíveis no vetor de dados
quantidade_total_pontos = length(vetor_amostras_funcao_y);

% Verifica se o passo de discretização h é estritamente positivo
if tamanho_passo_h <= 0
    error('Erro de validação: O tamanho do passo h deve ser estritamente positivo.');
end

% A fórmula de 3 pontos progressivos requer no mínimo 3 pontos para calcular ao menos uma derivada
if quantidade_total_pontos < 3
    error('Erro de validação: O vetor de amostras deve conter no mínimo 3 pontos (N >= 3).');
end

% -------------------------------------------------------------------------
% ETAPA 3: PROCESSAMENTO NUMÉRICO (DIFERENÇA FINITA PROGRESSIVA DE 3 PONTOS)
% -------------------------------------------------------------------------
% A fórmula progressiva de três pontos utiliza f(x_k), f(x_{k+1}) e f(x_{k+2}).
% Portanto, é avaliável do primeiro ponto (k = 1) até o antepenúltimo ponto (k = N - 2).
quantidade_pontos_calculaveis = quantidade_total_pontos - 2;

% Pré-aloca o vetor resultante das derivadas numéricas para máxima eficiência computacional
vetor_derivadas_numericas_dx = zeros(1, quantidade_pontos_calculaveis);

% Loop iterativo aplicando a quadratura diferencial em cada nó da malha
for indice_ponto_atual = 1 : quantidade_pontos_calculaveis
    % Amostra no ponto de avaliação: f(x_k)
    amostra_ponto_base = vetor_amostras_funcao_y(indice_ponto_atual);

    % Amostra no primeiro ponto avançado: f(x_k + h)
    amostra_primeiro_avanco = vetor_amostras_funcao_y(indice_ponto_atual + 1);

    % Amostra no segundo ponto avançado: f(x_k + 2h)
    amostra_segundo_avanco = vetor_amostras_funcao_y(indice_ponto_atual + 2);

    % Aplicação da fórmula progressiva de 3 pontos com parênteses obrigatórios no denominador:
    % dy/dx \approx (-3*y_k + 4*y_{k+1} - y_{k+2}) / (2 * h)
    vetor_derivadas_numericas_dx(indice_ponto_atual) = ...
        (-3.0 * amostra_ponto_base + 4.0 * amostra_primeiro_avanco - amostra_segundo_avanco) ...
        / (2.0 * tamanho_passo_h);
end

% -------------------------------------------------------------------------
% ETAPA 4: ANÁLISE DE PRECISÃO E ERROS CONTRA A SOLUÇÃO ANALÍTICA
% -------------------------------------------------------------------------
% Extrai as coordenadas correspondentes aos pontos onde a derivada foi estimada
vetor_coordenadas_calculadas_x = vetor_coordenadas_x(1 : quantidade_pontos_calculaveis);

% Calcula os valores exatos de referência através da derivada analítica
vetor_derivadas_analiticas_exatas = funcao_derivada_analitica_df(vetor_coordenadas_calculadas_x);

% Calcula o vetor de erros absolutos ponto a ponto: |f'_numérica - f'_exata|
vetor_erros_absolutos = abs(vetor_derivadas_numericas_dx - vetor_derivadas_analiticas_exatas);

% Determina o erro absoluto máximo e o erro médio em todo o domínio amostrado
erro_absoluto_maximo = max(vetor_erros_absolutos);
erro_absoluto_medio = mean(vetor_erros_absolutos);

% -------------------------------------------------------------------------
% ETAPA 5: EXIBIÇÃO DIDÁTICA DOS RESULTADOS NO CONSOLE
% -------------------------------------------------------------------------
fprintf('=========================================================================\n');
fprintf('     RELATÓRIO: DIFERENCIAÇÃO NUMÉRICA - PROGRESSIVA DE 3 PONTOS O(h^2)  \n');
fprintf('=========================================================================\n');
fprintf(' Função analisada                  : f(x) = sin(2*x)\n');
fprintf(' Derivada exata de referência      : f''(x) = 2*cos(2*x)\n');
fprintf(' Tamanho do passo uniforme (h)     : %16.8f\n', tamanho_passo_h);
fprintf(' Total de amostras no domínio (N)  : %16d\n', quantidade_total_pontos);
fprintf(' Pontos diferenciados (N - 2)      : %16d\n', quantidade_pontos_calculaveis);
fprintf(' Erro Absoluto Máximo              : %16.8e\n', erro_absoluto_maximo);
fprintf(' Erro Absoluto Médio               : %16.8e\n', erro_absoluto_medio);
fprintf('-------------------------------------------------------------------------\n');
fprintf(' AMOSTRAGEM COMPARATIVA (Primeiros nós da partição):\n');
fprintf('  Nó   |     x_k     |   f(x_k)   | f''_numérica |  f''_analítica | Erro Absoluto \n');
fprintf('-------+-------------+------------+--------------+--------------+---------------\n');
quantidade_amostras_exibicao = min(10, quantidade_pontos_calculaveis);
% Itera imprimindo cada linha correspondente aos primeiros nós
for indice_exibicao = 1 : quantidade_amostras_exibicao
    % Imprime os dados numéricos do nó atual na tabela comparativa
    fprintf('  %3d  | %11.6f | %10.6f | %12.6f | %12.6f | %13.6e\n', ...
            indice_exibicao, ...
            vetor_coordenadas_calculadas_x(indice_exibicao), ...
            vetor_amostras_funcao_y(indice_exibicao), ...
            vetor_derivadas_numericas_dx(indice_exibicao), ...
            vetor_derivadas_analiticas_exatas(indice_exibicao), ...
            vetor_erros_absolutos(indice_exibicao));
end
fprintf('=========================================================================\n\n');

% -------------------------------------------------------------------------
% ETAPA 6: VISUALIZAÇÃO GRÁFICA COMPARATIVA DIDÁTICA
% -------------------------------------------------------------------------
identificador_figura_grafico = figure('Name', 'Diferenciacao Numerica - Progressiva de 3 Pontos', 'NumberTitle', 'off');

% Subplot 1: Confronto da Derivada Analítica com a Derivada Numérica Calculada
subplot(2, 1, 1);
% Traça a derivada analítica como linha contínua suave
vetor_curva_suave_x = linspace(limite_inferior_dominio, limite_superior_dominio, 500);
plot(vetor_curva_suave_x, funcao_derivada_analitica_df(vetor_curva_suave_x), ...
     'b-', 'LineWidth', 1.8, 'DisplayName', 'Derivada Analítica f''(x) = 2 \cdot cos(2x)');
hold on;
% Plota as estimativas da derivada numérica como marcadores circulares
plot(vetor_coordenadas_calculadas_x, vetor_derivadas_numericas_dx, ...
     'ro', 'MarkerFaceColor', [1.0, 0.4, 0.4], 'MarkerSize', 4, ...
     'DisplayName', sprintf('Diferença Progressiva 3 Pontos (h = %.3f)', tamanho_passo_h));
grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6);
xlabel('Coordenada x', 'FontWeight', 'bold');
ylabel('Derivada f''(x)', 'FontWeight', 'bold');
title('Validação da Derivada Numérica vs. Analítica', 'FontSize', 11, 'FontWeight', 'bold');
legend('Location', 'northeast');
hold off;

% Subplot 2: Perfil do Erro Absoluto de Truncamento ao Longo do Domínio
subplot(2, 1, 2);
stem(vetor_coordenadas_calculadas_x, vetor_erros_absolutos, ...
     'Color', [0.8, 0.2, 0.2], 'LineWidth', 1.2, 'MarkerSize', 3, ...
     'DisplayName', 'Erro Absoluto |f''_num - f''_exata|');
hold on;
% Linha de referência indicando a ordem teórica do erro O(h^2)
linha_limite_erro_teorico = (tamanho_passo_h^2 / 3.0) * 8.0; % max |f'''(x)| = 8
plot([limite_inferior_dominio, limite_superior_dominio], ...
     [linha_limite_erro_teorico, linha_limite_erro_teorico], ...
     'k--', 'LineWidth', 1.2, 'DisplayName', sprintf('Cota Teórica Superior O(h^2) \\approx %.2e', linha_limite_erro_teorico));
grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6);
xlabel('Coordenada x', 'FontWeight', 'bold');
ylabel('Erro Absoluto', 'FontWeight', 'bold');
title('Distribuição Espacial do Erro de Truncamento O(h^2)', 'FontSize', 11, 'FontWeight', 'bold');
legend('Location', 'northeast');
hold off;
