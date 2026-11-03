% ==============================================================================
% UTILITÁRIOS MATEMÁTICOS: VISUALIZAÇÃO GRÁFICA DIDÁTICA DE FUNÇÕES ANALÍTICAS
% ==============================================================================
%
% 1. O QUE CALCULA:
%    Este script atua como uma rotina didática de visualização gráfica de curvas matemáticas e
%    análise de comportamento gráfico de funções univariadas f(x), unificando e
%    expandindo os estudos de cabos suspensos (catenárias) e de polinômios
%    ortogonais de aproximação numérica.
%
%    Casos de Estudo Implementados:
%    a) Curva Catenária Transcendental (Equação de Cabo de Transmissão Elétrica):
%       f_catenaria(x) = y_0 - c * cosh( x / c )
%       com parâmetros: y_0 = 211.50 m, c = 30.38 m.
%       Modela a catenária invertida ou flecha de sustentação sob tração uniforme.
%
%    b) Polinômio de Grau 5 e Derivada Analítica (Quadratura de Gauss-Legendre):
%       f_polinomial(x) = x^5 - (10/9)*x^3 + (5/21)*x
%       f_derivada(x)   = 5*x^4 - (10/3)*x^2 + (5/21)
%       Apresenta múltiplas raízes no intervalo [-1, +1], fundamentais para
%       integração numérica por pontos ótimos de Legendre-Gauss.
%
% 2. O QUE ENTREGA:
%    - Figura 1: Curva da Catenária no domínio [-30, +30] com ponto de vértice central
%      destacado, grade milimétrica, eixos cartesianos e anotações geométricas.
%    - Figura 2: Polinômio de grau 5 e sua derivada no domínio canônico [-1, +1],
%      com localização visual das raízes reais marcadas sobre o eixo nulo.
%    - Resumo quantitativo das propriedades críticas no console.
%
% 3. O QUE PRECISA TER (PRÉ-REQUISITOS E PARÂMETROS):
%    - Intervalos de visualização [x_min, x_max].
%    - Resolução da malha de amostragem (número de amostras suficiente para suavidade).
%
% 4. EXEMPLO DE USO E RESULTADO ESPERADO:
%    - Executar diretamente no terminal:
%      >> plotar_funcao
%    - Resultado esperado:
%      Geração de 2 janelas gráficas didáticas de alta legibilidade e relatório
%      numérico contendo os valores extremos das funções avaliadas.
% ==============================================================================

% Limpeza do ambiente de execução
clc;        % Limpa a tela do console de comandos
clear all;  % Remove todas as variáveis ativas do espaço de trabalho (workspace)
close all;  % Fecha todas as janelas de figuras abertas anteriormente

% ==============================================================================
% ESTUDO 1: CURVA CATENÁRIA (LINHA DE TRANSMISSÃO SOB TENSÃO MECÂNICA)
% ==============================================================================

% Parâmetros físicos e geométricos da catenária
parametro_altura_referencia_metros = 211.50; % Deslocamento vertical de ancoragem (y_0)
parametro_curvatura_escala_metros  = 30.38;  % Parâmetro característico da catenária (T / w)

% Definição da função anônima matemática da catenária
funcao_catenaria = @(x) parametro_altura_referencia_metros - ...
                         parametro_curvatura_escala_metros * cosh(x / parametro_curvatura_escala_metros);

% Limites do domínio de visualização espacial para a catenária
limite_inferior_catenaria_metros = -30.0;
limite_superior_catenaria_metros = +30.0;
numero_amostras_catenaria        = 300; % Malha densa para curva esteticamente suave

% Validação dos parâmetros geométricos
if parametro_curvatura_escala_metros <= 0.0
    error('O parâmetro de curvatura da catenária c deve ser estritamente positivo para evitar divisão por zero.');
end
if limite_superior_catenaria_metros <= limite_inferior_catenaria_metros
    error('O limite superior de visualização deve ser maior que o limite inferior.');
end

% Construção do vetor espacial x
vetor_espacial_catenaria_metros = linspace(limite_inferior_catenaria_metros, ...
                                           limite_superior_catenaria_metros, ...
                                           numero_amostras_catenaria);

% Avaliação numérica dos valores verticais da catenária
vetor_altura_catenaria_metros = funcao_catenaria(vetor_espacial_catenaria_metros);

% Determinação do vértice (ponto de máximo em x = 0 da função invertida)
posicao_vertice_x = 0.0;
posicao_vertice_y = funcao_catenaria(posicao_vertice_x);

% ------------------------------------------------------------------------------
% Construção Gráfica da Figura 1 (Catenária)
% ------------------------------------------------------------------------------
figura_catenaria = figure(1);
set(figura_catenaria, 'Name', 'Analise Grafica - Catenaria de Cabo Suspenso');
hold on;

% Traçado da curva matemática contínua principal
plot(vetor_espacial_catenaria_metros, vetor_altura_catenaria_metros, ...
     'Color', [0.0, 0.45, 0.75], 'LineWidth', 2.5);

% Marcador de destaque no vértice central
plot(posicao_vertice_x, posicao_vertice_y, ...
     'ro', 'MarkerSize', 9, 'MarkerFaceColor', [1.0, 0.2, 0.2]);

% Linha horizontal do nível do vértice como guia de flecha
plot([limite_inferior_catenaria_metros, limite_superior_catenaria_metros], ...
     [posicao_vertice_y, posicao_vertice_y], 'r:', 'LineWidth', 1.2);

grid on;
title('Perfil Geométrico de Curva Catenária: f(x) = 211.50 - 30.38 \cdot cosh(x / 30.38)', ...
      'FontSize', 12, 'FontWeight', 'bold');
xlabel('Posição Horizontal x (metros)', 'FontSize', 11);
ylabel('Altura Vertical f(x) (metros)', 'FontSize', 11);

identificador_legenda_figura_1 = legend('Catenária do Cabo de Transmissão', ...
                                        sprintf('Vértice Central (x=0, y=%.2f m)', posicao_vertice_y), ...
                                        'Nível de Referência da Flecha', ...
                                        'Location', 'south');
set(identificador_legenda_figura_1, 'FontSize', 10);
hold off;

% ==============================================================================
% ESTUDO 2: POLINÔMIO DE GRAU 5 E SUA DERIVADA (BASE DE GAUSS-LEGENDRE)
% ==============================================================================

% Definição da função analítica polinomial de grau 5
% f(x) = x^5 - (10/9)*x^3 + (5/21)*x
funcao_polinomial_grau_5 = @(x) x.^5 - (10.0 / 9.0) .* (x.^3) + (5.0 / 21.0) .* x;

% Derivada exata de primeira ordem: f'(x) = 5*x^4 - (10/3)*x^2 + (5/21)
funcao_derivada_primeira = @(x) 5.0 .* (x.^4) - (10.0 / 3.0) .* (x.^2) + (5.0 / 21.0);

% Intervalo canônico normalizado de análise
limite_inferior_polinomio = -1.0;
limite_superior_polinomio = +1.0;
numero_amostras_polinomio = 400;

vetor_coordenadas_x = linspace(limite_inferior_polinomio, limite_superior_polinomio, numero_amostras_polinomio);

% Avaliação do polinômio e de sua derivada na malha discreta
vetor_valores_polinomio = funcao_polinomial_grau_5(vetor_coordenadas_x);
vetor_valores_derivada  = funcao_derivada_primeira(vetor_coordenadas_x);

% Obtenção das raízes exatas do polinômio f(x) = 0 via roots
% Coeficientes: 1*x^5 + 0*x^4 - (10/9)*x^3 + 0*x^2 + (5/21)*x + 0
coeficientes_polinomio_grau_5 = [1.0, 0.0, -(10.0 / 9.0), 0.0, (5.0 / 21.0), 0.0];
vetor_raizes_reais_polinomio  = roots(coeficientes_polinomio_grau_5);

% Filtragem de raízes estritamente reais contidas no intervalo [-1, +1]
raizes_reais_validas = sort(real(vetor_raizes_reais_polinomio(abs(imag(vetor_raizes_reais_polinomio)) < 1e-12)));

% ------------------------------------------------------------------------------
% Construção Gráfica da Figura 2 (Polinômio e Derivada com Raízes)
% ------------------------------------------------------------------------------
figura_polinomio = figure(2);
set(figura_polinomio, 'Name', 'Analise Grafica - Polinomio de Legendre e Derivada');
hold on;

% Eixo cartesiano y = 0 em linha tracejada
plot([limite_inferior_polinomio, limite_superior_polinomio], [0, 0], ...
     'k--', 'LineWidth', 1.2);

% Traçado da função polinomial f(x)
plot(vetor_coordenadas_x, vetor_valores_polinomio, ...
     'Color', [0.1, 0.6, 0.2], 'LineWidth', 2.2);

% Traçado da função derivada f'(x)
plot(vetor_coordenadas_x, vetor_valores_derivada, ...
     'Color', [0.8, 0.4, 0.0], 'LineWidth', 1.8, 'LineStyle', '-.');

% Destaque visual das raízes reais de f(x) sobre o eixo nulo
plot(raizes_reais_validas, zeros(size(raizes_reais_validas)), ...
     'ro', 'MarkerSize', 8, 'MarkerFaceColor', [1.0, 0.0, 0.0]);

grid on;
title('Polinômio Ortogonal de Grau 5, Sua Derivada e Raízes no Intervalo [-1, +1]', ...
      'FontSize', 12, 'FontWeight', 'bold');
xlabel('Coordenada x (adimensional)', 'FontSize', 11);
ylabel('Amplitude da Função e Derivada', 'FontSize', 11);

identificador_legenda_figura_2 = legend('Eixo Nulo y = 0', ...
                                        'Polinômio f(x) = x^5 - \frac{10}{9}x^3 + \frac{5}{21}x', ...
                                        'Derivada f''(x) = 5x^4 - \frac{10}{3}x^2 + \frac{5}{21}', ...
                                        'Raízes Reais f(x) = 0', ...
                                        'Location', 'north');
set(identificador_legenda_figura_2, 'FontSize', 10);
hold off;

% ==============================================================================
% APRESENTAÇÃO DOS RESULTADOS NO CONSOLE
% ==============================================================================

fprintf('\n==============================================================================\n');
fprintf('     UTILITÁRIOS: RELATÓRIO DE ANÁLISE GRÁFICA DE FUNÇÕES ANALÍTICAS         \n');
fprintf('==============================================================================\n');
fprintf('  CASO 1: Curva Catenária Invertida\n');
fprintf('    - Ponto de Altura Máxima no Vértice (x = 0):   %10.4f m\n', posicao_vertice_y);
fprintf('    - Altura nos Extremos de Ancoragem (x = ±30): %10.4f m\n', funcao_catenaria(30.0));
fprintf('    - Flecha Total da Curva (y_max - y_min):      %10.4f m\n', ...
        posicao_vertice_y - funcao_catenaria(30.0));
fprintf('------------------------------------------------------------------------------\n');
fprintf('  CASO 2: Polinômio de Grau 5 e Raízes Calculadas via Roots:\n');
for indice_raiz = 1 : length(raizes_reais_validas)
    fprintf('    - Raiz Real %d:  x = %+10.6f  |  Resíduo f(x) = %11.4e\n', ...
            indice_raiz, ...
            raizes_reais_validas(indice_raiz), ...
            funcao_polinomial_grau_5(raizes_reais_validas(indice_raiz)));
end
fprintf('==============================================================================\n\n');
