% =========================================================================
% UNIVERSIDADE DE BRASÍLIA / CÁLCULO NUMÉRICO / COMPUTAÇÃO CIENTÍFICA
% SCRIPT: trapezio_composto.m
% OBJETIVO: Integração numérica pela Regra do Trapézio Composta (Newton-Cotes)
% AUTOR ORIGINAL: Davi Monteiro Alves
% REFATORAÇÃO: Engenharia Matemática e Especialista em Computação Científica
% =========================================================================
%
% -------------------------------------------------------------------------
% 1. O QUE CALCULA:
% -------------------------------------------------------------------------
%   Aproximação numérica da integral definida de uma função real f(x)
%   no intervalo fechado [limite_inferior_integracao, limite_superior_integracao]
%   utilizando a Regra do Trapézio Composta (Fórmula de Newton-Cotes fechada
%   de grau 1).
%
%   Formulação Matemática:
%     O intervalo [a, b] é subdividido em 'n' subintervalos de largura uniforme:
%       h = (b - a) / n
%     Gerando os nós de partição:
%       x_k = a + k * h,  com k = 0, 1, 2, ..., n
%
%     Em cada subintervalo [x_k, x_{k+1}], a função f(x) é aproximada pelo
%     polinômio interpolador de Lagrange de grau 1 (segmento de reta).
%     Integrando cada trapézio e somando as contribuições, obtém-se:
%       \int_{a}^{b} f(x) dx \approx h * [ (f(a) + f(b)) / 2 + \sum_{k=1}^{n-1} f(x_k) ]
%
%     Termo do Erro Global de Truncamento:
%       E_{trap} = - [ (b - a) / 12 ] * h^2 * f''(xi),  onde xi \in (a, b).
%       Ordem de convergência: O(h^2) (convergência de segunda ordem).
%
% -------------------------------------------------------------------------
% 2. O QUE ENTREGA:
% -------------------------------------------------------------------------
%   - valor_aproximado_integral: Escalar com o valor numérico estimado da integral.
%   - Relatório textual formatado no console com a discriminação dos parâmetros,
%     largura do passo h e valor final obtido.
%   - Gráfico didático bidimensional com a função contínua f(x), área sob a curva
%     e destaque para os limites de integração.
%
% -------------------------------------------------------------------------
% 3. O QUE PRECISA TER:
% -------------------------------------------------------------------------
%   - funcao_integrando: Handle de função matemática f(x), contínua em [a, b].
%   - limite_inferior_integracao: Escalar real definindo o início da integração (a).
%   - limite_superior_integracao: Escalar real definindo o término da integração (b), com b > a.
%   - numero_subintervalos: Número inteiro positivo (n >= 1) de subdivisões.
%
% -------------------------------------------------------------------------
% 4. EXEMPLO DE USO E RESULTADO ESPERADO:
% -------------------------------------------------------------------------
%   Função teste (comprimento de arco de catenária escalonada):
%     f(x) = sqrt(1 + 0.5 * sinh(x / 30)^2)
%     limite_inferior_integracao = -91.0
%     limite_superior_integracao = 91.0
%     numero_subintervalos      = 4000
%   Resultado esperado:
%     tamanho_passo_h           = 0.0455000000
%     valor_aproximado_integral = 467.6650834971
% =========================================================================

% -------------------------------------------------------------------------
% ETAPA 0: LIMPEZA DO AMBIENTE DE TRABALHO
% -------------------------------------------------------------------------
% Limpa a janela de comandos para facilitar a visualização dos resultados
clc;
% Remove todas as variáveis do espaço de trabalho (workspace)
clear all;
% Fecha todas as janelas gráficas que estejam abertas
close all;

% -------------------------------------------------------------------------
% ETAPA 1: DEFINIÇÃO DOS PARÂMETROS DE ENTRADA
% -------------------------------------------------------------------------
% Função integrando contínua f(x), implementada de forma vetorizada com operadores ponto (./ e .^)
funcao_integrando = @(coordenada_x) sqrt(1.0 + (0.5 .* (sinh(coordenada_x ./ 30.0)).^2));

% Limite inferior do intervalo de integração definida (a)
limite_inferior_integracao = -91.0;

% Limite superior do intervalo de integração definida (b)
limite_superior_integracao = 91.0;

% Quantidade total de subintervalos na partição uniforme (n)
numero_subintervalos = 4000;

% -------------------------------------------------------------------------
% ETAPA 2: VALIDAÇÃO PRÉVIA DE CONSISTÊNCIA MATEMÁTICA
% -------------------------------------------------------------------------
% Verifica se o limite inferior é estritamente menor que o superior
if limite_inferior_integracao >= limite_superior_integracao
    error('Erro de validação: O limite inferior deve ser estritamente menor que o limite superior.');
end

% Verifica se o número de subintervalos é um número inteiro válido e positivo
if (numero_subintervalos < 1) || (floor(numero_subintervalos) ~= numero_subintervalos)
    error('Erro de validação: O número de subintervalos deve ser um número inteiro estritamente positivo (>= 1).');
end

% -------------------------------------------------------------------------
% ETAPA 3: PROCESSAMENTO NUMÉRICO (QUADRATURA DO TRAPÉZIO COMPOSTA)
% -------------------------------------------------------------------------
% Calcula a largura constante de cada subintervalo: h = (b - a) / n
tamanho_passo_h = (limite_superior_integracao - limite_inferior_integracao) / numero_subintervalos;

% Constrói o vetor contendo todos os nós da partição equidistante: x_0, x_1, ..., x_n
vetor_pontos_particao = linspace(limite_inferior_integracao, limite_superior_integracao, numero_subintervalos + 1);

% Avalia a função integrando em cada nó da partição espacial
vetor_amostras_funcao = funcao_integrando(vetor_pontos_particao);

% Extrai a imagem da função no primeiro nó da partição: f(x_0) = f(a)
termo_fronteira_inicial = vetor_amostras_funcao(1);

% Extrai a imagem da função no último nó da partição: f(x_n) = f(b)
termo_fronteira_final = vetor_amostras_funcao(numero_subintervalos + 1);

% Soma os valores da função em todos os nós estritamente intermediários: \sum_{k=1}^{n-1} f(x_k)
soma_termos_intermediarios = sum(vetor_amostras_funcao(2 : numero_subintervalos));

% Aplica a fórmula da Quadratura Composta do Trapézio:
% I \approx h * [ (f(a) + f(b)) / 2 + \sum_{k=1}^{n-1} f(x_k) ]
valor_aproximado_integral = tamanho_passo_h * (((termo_fronteira_inicial + termo_fronteira_final) / 2.0) + soma_termos_intermediarios);

% -------------------------------------------------------------------------
% ETAPA 4: EXIBIÇÃO DIDÁTICA DOS RESULTADOS NO CONSOLE
% -------------------------------------------------------------------------
fprintf('=========================================================================\n');
fprintf('        RELATÓRIO: INTEGRAÇÃO NUMÉRICA - REGRA DO TRAPÉZIO COMPOSTA      \n');
fprintf('=========================================================================\n');
fprintf(' Limite inferior de integração (a) : %16.8f\n', limite_inferior_integracao);
fprintf(' Limite superior de integração (b) : %16.8f\n', limite_superior_integracao);
fprintf(' Número de subintervalos (n)       : %16d\n', numero_subintervalos);
fprintf(' Tamanho do passo uniforme (h)     : %16.10f\n', tamanho_passo_h);
fprintf('-------------------------------------------------------------------------\n');
fprintf(' VALOR APROXIMADO DA INTEGRAL (I)  : %16.10f\n', valor_aproximado_integral);
fprintf(' Ordem do erro de truncamento      : O(h^2) (Global)\n');
fprintf('=========================================================================\n\n');

% -------------------------------------------------------------------------
% ETAPA 5: VISUALIZAÇÃO GRÁFICA DIDÁTICA DA FUNÇÃO INTEGRANDO
% -------------------------------------------------------------------------
% Cria malha fina de coordenadas para traçar o gráfico suave da curva
vetor_coordenadas_curva_suave = linspace(limite_inferior_integracao, limite_superior_integracao, 1000);
vetor_valores_curva_suave = funcao_integrando(vetor_coordenadas_curva_suave);

% Inicializa a janela gráfica com identificador descritivo
identificador_figura_grafico = figure('Name', 'Integracao Numerica - Regra do Trapezio Composta', 'NumberTitle', 'off');

% Traça o preenchimento da área sob a curva f(x)
area(vetor_coordenadas_curva_suave, vetor_valores_curva_suave, ...
     'FaceColor', [0.85, 0.92, 1.0], 'EdgeColor', 'none', 'DisplayName', 'Área aproximada sob f(x)');
hold on;

% Traça a linha contínua representativa da função integrando f(x)
plot(vetor_coordenadas_curva_suave, vetor_valores_curva_suave, ...
     'Color', [0.0, 0.25, 0.65], 'LineWidth', 2.0, 'DisplayName', 'f(x) = \surd(1 + 0.5 \cdot sinh^2(x/30))');

% Ativa grade de fundo com estilo pontilhado para facilitar leitura de coordenadas
grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6);

% Configura os rótulos dos eixos coordenados com formatação matemática
xlabel('Coordenada x', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Imagem f(x)', 'FontSize', 11, 'FontWeight', 'bold');

% Configura o título com informações sintetizadas do cálculo
title(sprintf('Regra do Trapézio Composta (n = %d, h = %.4f) \\rightarrow I \\approx %.6f', ...
      numero_subintervalos, tamanho_passo_h, valor_aproximado_integral), ...
      'FontSize', 12, 'FontWeight', 'bold');

% Insere legenda didática
legend('Location', 'north');
hold off;
