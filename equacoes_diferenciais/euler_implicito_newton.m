% =========================================================================
% UNIVERSIDADE DE BRASÍLIA / CÁLCULO NUMÉRICO / COMPUTAÇÃO CIENTÍFICA
% SCRIPT: euler_implicito_newton.m
% OBJETIVO: Resolução de EDO via Euler Implícito (Backward Euler) e Newton-Raphson
% AUTOR ORIGINAL: Davi Monteiro Alves
% REFATORAÇÃO: Engenharia Matemática e Especialista em Computação Científica
% =========================================================================
%
% 1. O QUE CALCULA:
%    Resolve numericamente Problemas de Valor Inicial (PVI) para EDOs de 1a ordem
%    atraves do Metodo de Euler Implicito (Backward Euler):
%              dy/dx = f(x, y),   com condicao inicial y(x0) = y0
%
%    Fundamentacao teorica e estabilidade:
%    Diferente do Euler Explicito, a formulacao implicita avalia a derivada no
%    instante avancado x_{i+1}:
%              y_{i+1} = y_i + h * f(x_{i+1}, y_{i+1})
%    Esta formulacao e incondicionalmente estavel (A-estavel), sendo a escolha
%    padrao para equacoes diferenciais rigidas (stiff) e processos dissipativos.
%
%    Como a incognita y_{i+1} aparece dentro do argumento da funcao nao linear f,
%    define-se para cada passo i uma funcao de residuo algebrico G(u) = 0, onde u = y_{i+1}:
%              G(u) = u - y_i - h * f(x_{i+1}, u) = 0
%
%    A cada passo de integracao, o valor de u e encontrado iterativamente pelo
%    Metodo de Newton-Raphson:
%              u^{(k+1)} = u^{(k)} - [ G(u^{(k)}) / G'(u^{(k)}) ]
%    onde a derivada analitica do residuo em relacao a incognita u e:
%              G'(u) = 1 - h * [df/dy](x_{i+1}, u)
%
% 2. O QUE ENTREGA:
%    - Vetor da malha independente x discretizada.
%    - Vetor da solucao numerica y(x) com estabilidade incondicional.
%    - Vetor da solucao analitica exata para comparacao de acuracia.
%    - Vetor com a contagem de iteracoes necessarias de Newton em cada passo temporal.
%    - Erros absoluto e relativo percentual a cada no da malha.
%    - Tabela detalhada impressa no console com dados da convergencia de Newton.
%    - Graficos da trajetoria da solucao, curvas de erro e esforco iterativo de Newton.
%
% 3. O QUE PRECISA TER (PARAMETROS DE ENTRADA):
%    - funcao_derivada_f: funcao anonima f(x, y) representando dy/dx.
%    - funcao_derivada_parcial_df_dy: derivada parcial analitica df/dy(x, y).
%    - coordenada_inicial_x0: limite inferior do dominio.
%    - coordenada_final_xf: limite superior do dominio.
%    - condicao_inicial_y0: estado inicial da variavel dependente y(x0).
%    - tamanho_passo_h: espacamento da malha de discretizacao.
%    - tolerancia_newton_convergencia: criterio de parada para o residuo algebrico.
%    - numero_maximo_iteracoes_newton: teto de seguranca para evitar lacos infinitos.
%
% 4. EXEMPLO DE USO E RESULTADO ESPERADO:
%    PVI Nao Linear Rigido: dy/dx = 5*(1 - y^2),   y(0) = 0,   x em [0, 1], h = 0.05
%    Solucao analitica exata: y(x) = tanh(5*x)
%    Em x = 0.05:
%       Equacao do residuo: G(u) = u - 0 - 0.05*[5*(1 - u^2)] = 0.25*u^2 + u - 0.25 = 0
%       Convergencia de Newton atinge residuo < 1e-10 em apenas 3 iteracoes.
%       y(0.05) calculada = 0.236068 (Exato: tanh(0.25) = 0.244919).
% =========================================================================

% -------------------------------------------------------------------------
% ETAPA 0: LIMPEZA DO AMBIENTE DE EXECUCAO
% -------------------------------------------------------------------------
% Limpa o terminal de comandos
clc;
% Esvazia o espaco de variaveis da memoria
clear all;
% Fecha janelas de figuras graficas anteriores
close all;

% -------------------------------------------------------------------------
% ETAPA 1: DEFINICAO DOS PARAMETROS DO PROBLEMA E CONFIGURACOES DE NEWTON
% -------------------------------------------------------------------------
% Equacao diferencial nao linear dy/dx = f(x, y)
funcao_derivada_f = @(coordenada_x, variavel_y) (5.0 * (1.0 - variavel_y .^ 2));

% Derivada parcial analitica de f(x, y) em relacao a variavel dependente y: df/dy
funcao_derivada_parcial_df_dy = @(coordenada_x, variavel_y) (-10.0 * variavel_y);

% Solucao analitica exata de referencia: y(x) = tanh(5*x)
funcao_solucao_analitica_exata = @(coordenada_x) tanh(5.0 * coordenada_x);

% Limite inferior do intervalo de integracao
coordenada_inicial_x0 = 0.0;

% Limite superior do intervalo de integracao
coordenada_final_xf = 1.0;

% Condicao inicial conhecida y(x0) = y0
condicao_inicial_y0 = 0.0;

% Incremento uniforme da variavel independente (tamanho do passo temporal h)
tamanho_passo_h = 0.05;

% Tolerancia maxima permitida para convergencia do Metodo de Newton-Raphson
tolerancia_newton_convergencia = 1.0e-11;

% Numero maximo de iteracoes admitidas no resolvedor algebrico de Newton por passo
numero_maximo_iteracoes_newton = 100;

% -------------------------------------------------------------------------
% ETAPA 2: VALIDACAO PREVIA DE CONSISTENCIA MATEMATICA
% -------------------------------------------------------------------------
% Verifica se o limite final da variavel independente e estritamente superior ao inicial
if coordenada_final_xf <= coordenada_inicial_x0
    error('Erro de validacao: A coordenada final xf deve ser estritamente maior que a inicial x0.');
end

% Verifica se o tamanho do passo temporal h e estritamente positivo
if tamanho_passo_h <= 0.0
    error('Erro de validacao: O tamanho do passo h deve ser estritamente positivo (h > 0).');
end

% Verifica se a tolerancia de parada de Newton e estritamente positiva
if tolerancia_newton_convergencia <= 0.0
    error('Erro de validacao: A tolerancia de convergencia de Newton deve ser estritamente positiva.');
end

% Verifica se o teto de iteracoes de Newton e um inteiro valido maior ou igual a 1
if (numero_maximo_iteracoes_newton < 1) || (floor(numero_maximo_iteracoes_newton) ~= numero_maximo_iteracoes_newton)
    error('Erro de validacao: O numero maximo de iteracoes de Newton deve ser um inteiro positivo (>= 1).');
end

% Calcula e valida o numero total de subintervalos na malha discretizada
numero_total_passos_discretizacao = round((coordenada_final_xf - coordenada_inicial_x0) / tamanho_passo_h);
if numero_total_passos_discretizacao < 1
    error('Erro de validacao: O dominio e passo configurados geram menos de um passo de integracao.');
end

% -------------------------------------------------------------------------
% ETAPA 3: ALOCACAO DE ESTRUTURAS DE DADOS E CONDICOES INICIAIS
% -------------------------------------------------------------------------
% Vetor dos nos discretos da variavel independente
vetor_malha_independente_x = zeros(numero_total_passos_discretizacao + 1, 1);

% Vetor da solucao numerica de Euler Implicito
vetor_solucao_aproximada_y = zeros(numero_total_passos_discretizacao + 1, 1);

% Vetor com os valores da solucao teorica exata
vetor_solucao_analitica_exata_y = zeros(numero_total_passos_discretizacao + 1, 1);

% Vetor para armazenamento do erro absoluto pontual
vetor_erro_absoluto_pontual = zeros(numero_total_passos_discretizacao + 1, 1);

% Vetor para armazenamento do erro relativo percentual
vetor_erro_relativo_percentual = zeros(numero_total_passos_discretizacao + 1, 1);

% Vetor com a contagem de iteracoes gastas pelo Newton a cada passo
vetor_iteracoes_newton_por_passo = zeros(numero_total_passos_discretizacao, 1);

% Inicializacao do primeiro elemento com a condicao inicial prescrita
vetor_malha_independente_x(1) = coordenada_inicial_x0;
vetor_solucao_aproximada_y(1) = condicao_inicial_y0;
vetor_solucao_analitica_exata_y(1) = funcao_solucao_analitica_exata(coordenada_inicial_x0);
vetor_erro_absoluto_pontual(1) = 0.0;
vetor_erro_relativo_percentual(1) = 0.0;

% -------------------------------------------------------------------------
% ETAPA 4: LACO DE INTEGRACAO - EULER IMPLICITO COM NEWTON-RAPHSON POR PASSO
% -------------------------------------------------------------------------
for indice_passo_integracao = 1:numero_total_passos_discretizacao

    % Coordenada independente do ponto corrente x_i
    coordenada_x_corrente = vetor_malha_independente_x(indice_passo_integracao);

    % Valor da solucao no ponto corrente y_i
    valor_y_corrente = vetor_solucao_aproximada_y(indice_passo_integracao);

    % Avanca a coordenada independente de forma exata e sem acumulo de erro de arredondamento
    coordenada_x_avancada = coordenada_inicial_x0 + (indice_passo_integracao * tamanho_passo_h);
    vetor_malha_independente_x(indice_passo_integracao + 1) = coordenada_x_avancada;

    % Estimativa inicial para o algoritmo de Newton (chute preditor de Euler Explicito)
    % u^{(0)} = y_i + h * f(x_i, y_i)
    estimativa_incognita_u = valor_y_corrente + tamanho_passo_h * funcao_derivada_f(coordenada_x_corrente, valor_y_corrente);

    % Variavel booleana de controle de convergencia
    sinalizador_convergencia_atingida = false;

    % Sub-laco iterativo do Metodo de Newton-Raphson para solucao da equacao G(u) = 0
    for indice_iteracao_newton = 1:numero_maximo_iteracoes_newton

        % 1. Avalia o residuo da equacao algebrica nao linear de Euler Implicito:
        %    G(u) = u - y_i - h * f(x_{i+1}, u)
        valor_residuo_funcao_g = estimativa_incognita_u - valor_y_corrente - tamanho_passo_h * funcao_derivada_f(coordenada_x_avancada, estimativa_incognita_u);

        % 2. Avalia a derivada analitica do residuo G em relacao a u:
        %    G'(u) = 1 - h * (df/dy)(x_{i+1}, u)
        derivada_residuo_g_linha = 1.0 - tamanho_passo_h * funcao_derivada_parcial_df_dy(coordenada_x_avancada, estimativa_incognita_u);

        % Verificacao de seguranca numerica para evitar divisao por zero
        if abs(derivada_residuo_g_linha) < 1.0e-15
            error('Metodo de Newton falhou: Derivada proxima de zero (|G''(u)| < 1e-15) no passo %d, iteracao %d.', ...
                  indice_passo_integracao, indice_iteracao_newton);
        end

        % 3. Calcula o incremento corretivo de Newton: Delta_u = G(u) / G'(u)
        incremento_passo_newton = valor_residuo_funcao_g / derivada_residuo_g_linha;

        % 4. Atualiza a estimativa da incognita: u_{k+1} = u_k - Delta_u
        estimativa_incognita_u = estimativa_incognita_u - incremento_passo_newton;

        % Protecao contra divergencia numerica (valores NaN ou Infinitos)
        if isnan(estimativa_incognita_u) || isinf(estimativa_incognita_u)
            error('Divergencia numerica no Metodo de Newton: valor NaN/Inf detectado no passo %d, iteracao %d.', ...
                  indice_passo_integracao, indice_iteracao_newton);
        end

        % Avalia o residuo a posteriori apos a atualizacao
        valor_residuo_atualizado = estimativa_incognita_u - valor_y_corrente - tamanho_passo_h * funcao_derivada_f(coordenada_x_avancada, estimativa_incognita_u);

        % 5. Criterio de convergencia estrito (erro relativo da correcao ou norma do residuo)
        criterio_correcao_relativa = abs(incremento_passo_newton) / max(1.0, abs(estimativa_incognita_u));
        criterio_residuo_algebrico = abs(valor_residuo_atualizado);

        if (criterio_correcao_relativa < tolerancia_newton_convergencia) || (criterio_residuo_algebrico < tolerancia_newton_convergencia)
            sinalizador_convergencia_atingida = true;
            vetor_iteracoes_newton_por_passo(indice_passo_integracao) = indice_iteracao_newton;
            break;
        end

    % Fim do sub-laco de Newton
    end

    % Protecao contra estouro de iteracoes sem convergencia
    if ~sinalizador_convergencia_atingida
        warning('Newton-Raphson atingiu o teto maximo de %d iteracoes no passo %d (x = %.4f) sem atingir a tolerancia especificada (%.2e).', ...
                numero_maximo_iteracoes_newton, indice_passo_integracao, coordenada_x_avancada, tolerancia_newton_convergencia);
        vetor_iteracoes_newton_por_passo(indice_passo_integracao) = numero_maximo_iteracoes_newton;
    end

    % Atribui o valor convergido de u como a solucao no proximo passo: y_{i+1} = u
    vetor_solucao_aproximada_y(indice_passo_integracao + 1) = estimativa_incognita_u;

    % Avaliacao teorica da solucao analitica exata
    valor_solucao_exata_avancada = funcao_solucao_analitica_exata(coordenada_x_avancada);
    vetor_solucao_analitica_exata_y(indice_passo_integracao + 1) = valor_solucao_exata_avancada;

    % Calculo do erro absoluto pontual
    erro_absoluto_passo = abs(valor_solucao_exata_avancada - estimativa_incognita_u);
    vetor_erro_absoluto_pontual(indice_passo_integracao + 1) = erro_absoluto_passo;

    % Calculo do erro relativo percentual
    if abs(valor_solucao_exata_avancada) > 1.0e-12
        vetor_erro_relativo_percentual(indice_passo_integracao + 1) = (erro_absoluto_passo / abs(valor_solucao_exata_avancada)) * 100.0;
    else
        vetor_erro_relativo_percentual(indice_passo_integracao + 1) = erro_absoluto_passo * 100.0;
    end

% Fim do laco de integracao temporal
end

% -------------------------------------------------------------------------
% ETAPA 5: RELATORIO TABULAR DE RESULTADOS E CONVERGENCIA
% -------------------------------------------------------------------------
fprintf('=========================================================================================================\n');
fprintf('     EULER IMPLICITO (BACKWARD EULER) COM NEWTON-RAPHSON PARA EDO NAO LINEAR RIGIDA                     \n');
fprintf('=========================================================================================================\n');
fprintf(' Parametros: x0 = %.4f | xf = %.4f | y0 = %.4f | Passo h = %.4f | Tolerancia Newton = %.2e\n', ...
        coordenada_inicial_x0, coordenada_final_xf, condicao_inicial_y0, tamanho_passo_h, tolerancia_newton_convergencia);
fprintf('---------------------------------------------------------------------------------------------------------\n');
fprintf('  Passo  |      x      |  y (Implicito) |   y (Exato)   | Iteracoes Newton | Erro Absoluto  |  Erro (%%) \n');
fprintf('---------------------------------------------------------------------------------------------------------\n');

% Exibe o ponto inicial
fprintf('  %4d   |  %9.4f  |  %12.6f  |  %12.6f  |     --------     |  %12.6e  |  %8.4f %%\n', ...
        0, vetor_malha_independente_x(1), vetor_solucao_aproximada_y(1), ...
        vetor_solucao_analitica_exata_y(1), vetor_erro_absoluto_pontual(1), vetor_erro_relativo_percentual(1));

% Exibe os passos calculados e o esforco iterativo de Newton
for indice_exibicao_passo = 1:numero_total_passos_discretizacao
    fprintf('  %4d   |  %9.4f  |  %12.6f  |  %12.6f  |       %3d        |  %12.6e  |  %8.4f %%\n', ...
            indice_exibicao_passo, ...
            vetor_malha_independente_x(indice_exibicao_passo + 1), ...
            vetor_solucao_aproximada_y(indice_exibicao_passo + 1), ...
            vetor_solucao_analitica_exata_y(indice_exibicao_passo + 1), ...
            vetor_iteracoes_newton_por_passo(indice_exibicao_passo), ...
            vetor_erro_absoluto_pontual(indice_exibicao_passo + 1), ...
            vetor_erro_relativo_percentual(indice_exibicao_passo + 1));
end
fprintf('=========================================================================================================\n');

% -------------------------------------------------------------------------
% ETAPA 6: GERACAO DE PAINEL GRAFICO DE DIAGNOSTICO E CONVERGENCIA
% -------------------------------------------------------------------------
figura_euler_implicito = figure('Name', 'Euler Implicito com Newton-Raphson', 'NumberTitle', 'off');

% Painel 1: Trajetoria da solucao comparada com a analitica
subplot(3, 1, 1);
plot(vetor_malha_independente_x, vetor_solucao_analitica_exata_y, 'k-', 'LineWidth', 2.0);
hold on;
plot(vetor_malha_independente_x, vetor_solucao_aproximada_y, 'r--o', 'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerFaceColor', 'r');
grid on;
title('Solucao de EDO Nao Linear: Euler Implicito (A-Estavel) vs. Solucao Exata', 'FontSize', 11, 'FontWeight', 'bold');
xlabel('Variavel Independente (x)', 'FontSize', 10);
ylabel('Variavel Dependente y(x)', 'FontSize', 10);
legend('Solucao Analitica tanh(5x)', 'Euler Implicito + Newton', 'Location', 'SouthEast');
hold off;

% Painel 2: Perfil do Erro Absoluto
subplot(3, 1, 2);
plot(vetor_malha_independente_x, vetor_erro_absoluto_pontual, 'b-d', 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', 'b');
grid on;
title('Evolucao do Erro Absoluto de Truncamento Global', 'FontSize', 11, 'FontWeight', 'bold');
xlabel('Variavel Independente (x)', 'FontSize', 10);
ylabel('Erro Absoluto |y_{exato} - y_{num}|', 'FontSize', 10);
legend('Erro Absoluto Pontual', 'Location', 'NorthEast');

% Painel 3: Numero de iteracoes de Newton por passo de integracao
subplot(3, 1, 3);
stem(vetor_malha_independente_x(2:end), vetor_iteracoes_newton_por_passo, 'fill', 'LineWidth', 1.5, 'Color', [0 0.5 0]);
grid on;
title('Esforco Computacional: Iteracoes do Metodo de Newton por Passo Temporal', 'FontSize', 11, 'FontWeight', 'bold');
xlabel('Variavel Independente (x)', 'FontSize', 10);
ylabel('Iteracoes de Newton', 'FontSize', 10);
ylim([0, max(vetor_iteracoes_newton_por_passo) + 2]);

% =========================================================================
% FIM DO SCRIPT euler_implicito_newton.m
% =========================================================================
