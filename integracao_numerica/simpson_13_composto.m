% =========================================================================
% UNIVERSIDADE DE BRASÍLIA / CÁLCULO NUMÉRICO / COMPUTAÇÃO CIENTÍFICA
% SCRIPT: simpson_13_composto.m
% OBJETIVO: Integração numérica pela Regra de Simpson 1/3 Composta (Newton-Cotes)
% AUTOR ORIGINAL: Davi Monteiro Alves
% REFATORAÇÃO: Engenharia Matemática e Especialista em Computação Científica
% =========================================================================
%
% -------------------------------------------------------------------------
% 1. O QUE CALCULA:
% -------------------------------------------------------------------------
%   Aproximação numérica de alta ordem para a integral definida de uma função
%   real f(x) no intervalo [limite_inferior_integracao, limite_superior_integracao]
%   utilizando a Regra de Simpson 1/3 Composta (Fórmula de Newton-Cotes fechada
%   de grau 2).
%
%   Formulação Matemática:
%     O intervalo [a, b] é dividido em 'n' subintervalos de largura uniforme h:
%       h = (b - a) / n
%     IMPORTANTE: A Regra de Simpson 1/3 exige obrigatoriamente que 'n' seja um
%     número PAR (múltiplo de 2), pois cada parábola interpoladora aproxima a
%     função f(x) sobre pares de subintervalos contíguos [x_{2j-2}, x_{2j}].
%
%     A quadratura combina as imagens ponderadas pelos pesos de Newton-Cotes:
%       \int_{a}^{b} f(x) dx \approx \frac{h}{3} \left[ f(x_0) + 4 \sum_{k=1,3,5}^{n-1} f(x_k) + 2 \sum_{k=2,4,6}^{n-2} f(x_k) + f(x_n) \right]
%
%     Termo do Erro Global de Truncamento:
%       E_{simp} = - [ (b - a) / 180 ] * h^4 * f^{(4)}(xi),  onde xi \in (a, b).
%       Ordem de convergência: O(h^4) (convergência de quarta ordem).
%       Exatidão algébrica: Exato para polinômios de grau até 3.
%
% -------------------------------------------------------------------------
% 2. O QUE ENTREGA:
% -------------------------------------------------------------------------
%   - valor_aproximado_integral: Escalar com o valor numérico estimado da integral.
%   - valor_integral_convertido_pes: Valor convertido para pés (fator 3.28 m -> ft)
%     referente ao cálculo clássico de comprimento de catenária em engenharia.
%   - Relatório textual formatado no console detalhando a contribuição dos nós
%     de índices ímpares (peso 4) e de índices pares (peso 2).
%   - Gráfico didático bidimensional com a função contínua f(x), área sob a curva
%     e destaque para os limites de integração.
%
% -------------------------------------------------------------------------
% 3. O QUE PRECISA TER:
% -------------------------------------------------------------------------
%   - funcao_integrando: Handle de função matemática f(x), contínua em [a, b].
%   - limite_inferior_integracao: Escalar real definindo o início da integração (a).
%   - limite_superior_integracao: Escalar real definindo o término da integração (b), com b > a.
%   - numero_subintervalos: Número inteiro POSITIVO e OBRIGATORIAMENTE PAR (n >= 2).
%
% -------------------------------------------------------------------------
% 4. EXEMPLO DE USO E RESULTADO ESPERADO:
% -------------------------------------------------------------------------
%   Função teste (comprimento de arco de linha de transmissão/catenária):
%     f(x) = sqrt(1 + 0.476454 * sinh(x / 30.38)^2)
%     limite_inferior_integracao = -91.21
%     limite_superior_integracao = 91.21
%     numero_subintervalos      = 20 (PAR)
%   Resultado esperado:
%     tamanho_passo_h                 = 9.1210000000
%     valor_aproximado_integral       = 451.4159383748
%     valor_integral_convertido_pes   = 1480.6442778693 (usando fator 3.28)
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
% Função integrando contínua f(x), vetorizada com operadores elementares (.^ e ./)
funcao_integrando = @(coordenada_x) sqrt(1.0 + (0.476454 .* (sinh(coordenada_x ./ 30.38)).^2));

% Limite inferior do intervalo de integração definida (a)
limite_inferior_integracao = -91.21;

% Limite superior do intervalo de integração definida (b)
limite_superior_integracao = 91.21;

% Quantidade de subintervalos na partição (n) - DEVE SER PAR
numero_subintervalos = 20;

% Fator multiplicativo de conversão física (metros para pés no contexto da aplicação)
fator_conversao_unidades = 3.28;

% -------------------------------------------------------------------------
% ETAPA 2: VALIDAÇÃO PRÉVIA DE CONSISTÊNCIA MATEMÁTICA
% -------------------------------------------------------------------------
% Verifica se o limite inferior é estritamente menor que o limite superior
if limite_inferior_integracao >= limite_superior_integracao
    error('Erro de validação: O limite inferior deve ser estritamente menor que o superior.');
end

% Verifica se o número de subintervalos é inteiro positivo e no mínimo 2
if (numero_subintervalos < 2) || (floor(numero_subintervalos) ~= numero_subintervalos)
    error('Erro de validação: O número de subintervalos deve ser um inteiro positivo maior ou igual a 2.');
end

% Validação fundamental da Regra de Simpson 1/3 Composta: o número de subintervalos DEVE ser PAR
if mod(numero_subintervalos, 2) ~= 0
    error('Erro de validação: Na Regra de Simpson 1/3 Composta, o número de subintervalos (n) DEVE SER PAR (múltiplo de 2). Valor fornecido: %d', numero_subintervalos);
end

% -------------------------------------------------------------------------
% ETAPA 3: PROCESSAMENTO NUMÉRICO (QUADRATURA DE SIMPSON 1/3 COMPOSTA)
% -------------------------------------------------------------------------
% Calcula a largura constante de cada subintervalo: h = (b - a) / n
tamanho_passo_h = (limite_superior_integracao - limite_inferior_integracao) / numero_subintervalos;

% Avalia a contribuição dos nós extremos de fronteira: f(x_0) e f(x_n)
termo_fronteira_inicial = funcao_integrando(limite_inferior_integracao);
termo_fronteira_final = funcao_integrando(limite_superior_integracao);

% Inicializa a soma acumulada com os nós de contorno (pesos unitários)
soma_total_quadratura = termo_fronteira_inicial + termo_fronteira_final;

% Acumula os termos correspondentes aos nós internos de índices ÍMPARES (peso 4)
% Índices matemáticos: k = 1, 3, 5, ..., n - 1
soma_termos_indices_impares = 0.0;
for indice_no_impar = 1 : 2 : (numero_subintervalos - 1)
    % Calcula a coordenada espacial x_k do nó ímpar atual
    coordenada_no_atual = limite_inferior_integracao + (indice_no_impar * tamanho_passo_h);
    % Soma a avaliação funcional f(x_k) ao acumulador dos nós ímpares
    soma_termos_indices_impares = soma_termos_indices_impares + funcao_integrando(coordenada_no_atual);
end

% Acumula os termos correspondentes aos nós internos de índices PARES (peso 2)
% Índices matemáticos: k = 2, 4, 6, ..., n - 2
soma_termos_indices_pares = 0.0;
for indice_no_par = 2 : 2 : (numero_subintervalos - 2)
    % Calcula a coordenada espacial x_k do nó par atual
    coordenada_no_atual = limite_inferior_integracao + (indice_no_par * tamanho_passo_h);
    % Soma a avaliação funcional f(x_k) ao acumulador dos nós pares
    soma_termos_indices_pares = soma_termos_indices_pares + funcao_integrando(coordenada_no_atual);
end

% Incorpora as somas parciais ponderadas à quadratura total
soma_total_quadratura = soma_total_quadratura + (4.0 * soma_termos_indices_impares) + (2.0 * soma_termos_indices_pares);

% Multiplica pelo fator h / 3 para obter o valor aproximado final da integral
valor_aproximado_integral = (tamanho_passo_h / 3.0) * soma_total_quadratura;

% Converte o valor calculado para a unidade derivada (ex: pés) caso aplicável
valor_integral_convertido_pes = valor_aproximado_integral * fator_conversao_unidades;

% -------------------------------------------------------------------------
% ETAPA 4: EXIBIÇÃO DIDÁTICA DOS RESULTADOS NO CONSOLE
% -------------------------------------------------------------------------
fprintf('=========================================================================\n');
fprintf('       RELATÓRIO: INTEGRAÇÃO NUMÉRICA - REGRA DE SIMPSON 1/3 COMPOSTA   \n');
fprintf('=========================================================================\n');
fprintf(' Limite inferior de integração (a)       : %16.8f\n', limite_inferior_integracao);
fprintf(' Limite superior de integração (b)       : %16.8f\n', limite_superior_integracao);
fprintf(' Número de subintervalos (n - PAR)       : %16d\n', numero_subintervalos);
fprintf(' Tamanho do passo uniforme (h)           : %16.10f\n', tamanho_passo_h);
fprintf('-------------------------------------------------------------------------\n');
fprintf(' Soma dos termos de fronteira f(a)+f(b)  : %16.10f\n', termo_fronteira_inicial + termo_fronteira_final);
fprintf(' Contribuição dos nós ímpares (4 * soma) : %16.10f\n', 4.0 * soma_termos_indices_impares);
fprintf(' Contribuição dos nós pares   (2 * soma) : %16.10f\n', 2.0 * soma_termos_indices_pares);
fprintf('-------------------------------------------------------------------------\n');
fprintf(' VALOR DA INTEGRAL NO SI (metros)        : %16.10f\n', valor_aproximado_integral);
fprintf(' VALOR CONVERTIDO (fator %.2f -> pés)     : %16.10f\n', fator_conversao_unidades, valor_integral_convertido_pes);
fprintf(' Ordem do erro de truncamento            : O(h^4) (Global)\n');
fprintf('=========================================================================\n\n');

% -------------------------------------------------------------------------
% ETAPA 5: VISUALIZAÇÃO GRÁFICA DIDÁTICA DA FUNÇÃO INTEGRANDO
% -------------------------------------------------------------------------
% Gera vetor de coordenadas para a curva suave contínua
vetor_coordenadas_curva_suave = linspace(limite_inferior_integracao, limite_superior_integracao, 1000);
vetor_valores_curva_suave = funcao_integrando(vetor_coordenadas_curva_suave);

% Cria janela de visualização gráfica
identificador_figura_grafico = figure('Name', 'Integracao Numerica - Regra de Simpson 1/3 Composta', 'NumberTitle', 'off');

% Desenha o preenchimento da área sob a curva
area(vetor_coordenadas_curva_suave, vetor_valores_curva_suave, ...
     'FaceColor', [0.88, 0.95, 0.88], 'EdgeColor', 'none', 'DisplayName', 'Área aproximada sob f(x)');
hold on;

% Plota a curva contínua da função
plot(vetor_coordenadas_curva_suave, vetor_valores_curva_suave, ...
     'Color', [0.0, 0.5, 0.2], 'LineWidth', 2.0, 'DisplayName', 'f(x) = \surd(1 + 0.476454 \cdot sinh^2(x/30.38))');

% Plota os nós discretos da partição de Simpson
vetor_coordenadas_nos_particao = linspace(limite_inferior_integracao, limite_superior_integracao, numero_subintervalos + 1);
vetor_valores_nos_particao = funcao_integrando(vetor_coordenadas_nos_particao);
plot(vetor_coordenadas_nos_particao, vetor_valores_nos_particao, ...
     'ro', 'MarkerFaceColor', [0.85, 0.1, 0.1], 'MarkerSize', 5, 'DisplayName', 'Nós da Quadratura (n=20)');

% Ativa grade de fundo com estilo pontilhado
grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6);

% Configura os rótulos dos eixos
xlabel('Coordenada x', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Imagem f(x)', 'FontSize', 11, 'FontWeight', 'bold');

% Configura o título do gráfico com os resultados chave
title(sprintf('Regra de Simpson 1/3 Composta (n = %d, h = %.4f) \\rightarrow I \\approx %.4f (%.2f ft)', ...
      numero_subintervalos, tamanho_passo_h, valor_aproximado_integral, valor_integral_convertido_pes), ...
      'FontSize', 12, 'FontWeight', 'bold');

% Adiciona legenda didática
legend('Location', 'north');
hold off;
