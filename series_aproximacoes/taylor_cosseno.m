% =========================================================================
% UNIVERSIDADE DE BRASÍLIA / CÁLCULO NUMÉRICO / COMPUTAÇÃO CIENTÍFICA
% SCRIPT: taylor_cosseno.m
% OBJETIVO: Expansão em Série de Taylor/Maclaurin para Função Cosseno cos(x)
% AUTOR ORIGINAL: Davi Monteiro Alves
% REFATORAÇÃO: Engenharia Matemática e Especialista em Computação Científica
% =========================================================================
%
% 1. O QUE CALCULA:
%    Aproxima o valor do cosseno de um angulo x (em radianos) atraves do
%    somatorio de sua Serie de Taylor truncada em torno da origem x0 = 0
%    (formalmente denominada Serie de Maclaurin):
%              cos(x) = \sum_{n=0}^{\infty} \frac{(-1)^n}{(2n)!} x^{2n}
%                     = 1 - \frac{x^2}{2!} + \frac{x^4}{4!} - \frac{x^6}{6!} + \dots
%
%    Fundamentacao matematica e correcao de alternancia:
%    - Para a funcao f(x) = cos(x), as derivadas impares avaliadas na origem
%      sao nulas: f'(0) = -sin(0) = 0, f'''(0) = sin(0) = 0.
%    - As derivadas de ordem par f^{(2n)}(0) alternam entre +1 e -1:
%      n = 0 -> cos(0) = +1
%      n = 1 -> -cos(0) = -1
%      n = 2 -> +cos(0) = +1
%      n = 3 -> -cos(0) = -1
%    - A alternancia de sinais e rigorosamente implementada atraves da paridade:
%      se mod(n, 2) == 0 -> termo positivo (+1)
%      se mod(n, 2) == 1 -> termo negativo (-1)
%      ou diretamente via (-1)^n.
%    - Pelo Teorema de Taylor com resto de Lagrange, o erro de truncamento
%      apos N termos e limitado por: |R_{2N+1}(x)| <= |x|^{2N+2} / (2N+2)!
%
% 2. O QUE ENTREGA:
%    - Exibicao numerica em formato de alta precisao (format long).
%    - Tabela passo a passo a cada termo adicionado:
%      * Indice do termo n e grau do monomio correspondente (2n).
%      * Sinal do termo (+ ou -).
%      * Magnitude do termo individual |x^{2n} / (2n)!|.
%      * Somatorio acumulado parcial da serie.
%      * Valor analitico de referencia da funcao embutida cos(x).
%      * Erro absoluto pontual |cos(x) - soma|.
%      * Erro relativo percentual (|cos(x) - soma| / |cos(x)|) * 100%.
%    - Grafico da evolucao do erro relativo percentual em escala logaritmica.
%
% 3. O QUE PRECISA TER (PARAMETROS DE ENTRADA):
%    - angulo_entrada_x_radianos: angulo a ser avaliado (em radianos).
%    - numero_termos_expansao_n: quantidade de termos considerados na serie.
%
% 4. EXEMPLO DE USO E RESULTADO ESPERADO:
%    Para x = pi/3 radianos (60 graus) e 6 termos (n de 0 a 5, grau ate 10):
%    Valor exato embutido: cos(pi/3) = 0.500000000000000
%    Termo n=0 (grau  0): +1.000000000000000 -> Soma = 1.000000000000000 (Erro = 100.000 %)
%    Termo n=1 (grau  2): -0.548311355616075 -> Soma = 0.451688644383925 (Erro = 9.662271 %)
%    Termo n=2 (grau  4): +0.050107557116256 -> Soma = 0.501796201500181 (Erro = 0.359240 %)
%    Termo n=3 (grau  6): -0.001831636171268 -> Soma = 0.499964565328913 (Erro = 0.007087 %)
%    Termo n=4 (grau  8): +0.000035868104002 -> Soma = 0.500000433432915 (Erro = 8.668658e-05 %)
%    Termo n=5 (grau 10): -0.000000437041972 -> Soma = 0.499999996390943 (Erro = 7.218114e-07 %)
% =========================================================================

% -------------------------------------------------------------------------
% ETAPA 0: LIMPEZA DO AMBIENTE E CONFIGURACAO DE PRECISAO
% -------------------------------------------------------------------------
% Limpa o terminal de comandos
clc;
% Esvazia o espaco de memoria
clear all;
% Fecha janelas de figuras graficas preexistentes
close all;
% Configura a exibicao numerica para ponto flutuante de precisao dupla (15 digitos)
format long;

% -------------------------------------------------------------------------
% ETAPA 1: DEFINICAO DOS PARAMETROS DE ENTRADA DA SERIE
% -------------------------------------------------------------------------
% Angulo de entrada a ser avaliado em radianos (exemplo: pi / 3 = 60 graus)
angulo_entrada_x_radianos = pi / 3.0;

% Quantidade de termos a incluir na expansao polinomial truncada (n de 0 ate N-1)
numero_termos_expansao_n = 8;

% -------------------------------------------------------------------------
% ETAPA 2: VALIDACAO PREVIA DE CONSISTENCIA MATEMATICA
% -------------------------------------------------------------------------
% Verifica se o numero de termos e um inteiro positivo valido
if (numero_termos_expansao_n < 1) || (floor(numero_termos_expansao_n) ~= numero_termos_expansao_n)
    error('Erro de validacao: O numero de termos deve ser um inteiro positivo maior ou igual a 1.');
end

% Protecao contra estouro de ponto flutuante em fatorial: (2n)! ultrapassa limite de double se grau > 170 (n > 85)
if numero_termos_expansao_n > 85
    error('Erro de validacao: O numero maximo de termos suportado em precisao dupla sem overflow fatorial e 85 (grau 170).');
end

% Verifica se o angulo fornecido e um valor numerico real e finito
if isnan(angulo_entrada_x_radianos) || isinf(angulo_entrada_x_radianos) || ~isreal(angulo_entrada_x_radianos)
    error('Erro de validacao: O angulo de entrada deve ser um escalar real e finito.');
end

% -------------------------------------------------------------------------
% ETAPA 3: PRE-ALOCACAO DE VETORES E PREPARACAO MATEMATICA
% -------------------------------------------------------------------------
% Avalia o valor de referencia fornecido pela funcao trigonometrica embutida
valor_referencia_cos_embutido = cos(angulo_entrada_x_radianos);

% Pre-alocacao de memoria para vetores de analise estatistica de convergencia
vetor_indices_n = zeros(numero_termos_expansao_n, 1);
vetor_graus_polinomio = zeros(numero_termos_expansao_n, 1);
vetor_termos_individuais = zeros(numero_termos_expansao_n, 1);
vetor_somas_acumuladas = zeros(numero_termos_expansao_n, 1);
vetor_erros_absolutos = zeros(numero_termos_expansao_n, 1);
vetor_erros_relativos_percentuais = zeros(numero_termos_expansao_n, 1);

% Inicializa a variavel acumuladora da serie com zero
somatorio_acumulado_serie = 0.0;

% -------------------------------------------------------------------------
% ETAPA 4: COMPUTACAO TERMO A TERMO DA SERIE DE MACLAURIN
% -------------------------------------------------------------------------
for indice_termo_n = 0:(numero_termos_expansao_n - 1)

    % Posicao no vetor indexado em 1 (para compatibilidade Octave/MATLAB)
    posicao_vetorial = indice_termo_n + 1;

    % O grau da potencia de x para o termo n e sempre par: 2*n
    grau_polinomio_termo = 2 * indice_termo_n;

    % Calcula a potencia par do angulo no numerador: x^(2*n)
    potencia_numerador_x = angulo_entrada_x_radianos ^ grau_polinomio_termo;

    % Calcula o fatorial do grau no denominador: (2*n)!
    fatorial_denominador = factorial(grau_polinomio_termo);

    % Calcula a magnitude escalar absoluta do termo individual
    magnitude_termo_individual = potencia_numerador_x / fatorial_denominador;

    % CORRECAO CRITICA DO SINAL:
    % Utiliza o operador de modulo mod(n, 2) para alternancia exata:
    % Se n for par (0, 2, 4, ...), o sinal e positivo (+1)
    % Se n for impar (1, 3, 5, ...), o sinal e negativo (-1)
    if mod(indice_termo_n, 2) == 0
        fator_sinal_alternado = +1.0;
    else
        fator_sinal_alternado = -1.0;
    end

    % Termo assinado completo da serie: (-1)^n * x^(2n) / (2n)!
    termo_assinado_completo = fator_sinal_alternado * magnitude_termo_individual;

    % Adiciona o termo atual ao somatorio cumulativo da serie
    somatorio_acumulado_serie = somatorio_acumulado_serie + termo_assinado_completo;

    % Calcula o erro absoluto pontual em relacao ao valor de referencia
    erro_absoluto_atual = abs(valor_referencia_cos_embutido - somatorio_acumulado_serie);

    % Calcula o erro relativo percentual
    if abs(valor_referencia_cos_embutido) > 1.0e-15
        erro_relativo_percentual_atual = (erro_absoluto_atual / abs(valor_referencia_cos_embutido)) * 100.0;
    else
        erro_relativo_percentual_atual = erro_absoluto_atual * 100.0;
    end

    % Registra os dados nos vetores historicos para tabulacao e graficos
    vetor_indices_n(posicao_vetorial) = indice_termo_n;
    vetor_graus_polinomio(posicao_vetorial) = grau_polinomio_termo;
    vetor_termos_individuais(posicao_vetorial) = termo_assinado_completo;
    vetor_somas_acumuladas(posicao_vetorial) = somatorio_acumulado_serie;
    vetor_erros_absolutos(posicao_vetorial) = erro_absoluto_atual;
    vetor_erros_relativos_percentuais(posicao_vetorial) = erro_relativo_percentual_atual;

% Fim do laco de somatorio da serie
end

% -------------------------------------------------------------------------
% ETAPA 5: EXIBICAO DOS RESULTADOS TABELADOS COM FORMAT LONG
% -------------------------------------------------------------------------
fprintf('=========================================================================================================================\n');
fprintf('                           EXPANSAO EM SERIE DE TAYLOR / MACLAURIN PARA cos(x)                                          \n');
fprintf('=========================================================================================================================\n');
fprintf(' Angulo avaliado (x): %.15f radianos (equivalente a %.6f graus)\n', ...
        angulo_entrada_x_radianos, angulo_entrada_x_radianos * (180.0 / pi));
fprintf(' Valor teorico exato [cos(x) embutido]: %.15f\n', valor_referencia_cos_embutido);
fprintf(' Quantidade de termos calculados: %d (polinomio de grau maximo %d)\n', ...
        numero_termos_expansao_n, vetor_graus_polinomio(end));
fprintf('-------------------------------------------------------------------------------------------------------------------------\n');
fprintf('  n  | Grau | Sinal |      Valor do Termo      |       Soma Parcial       |    Erro Absoluto    |   Erro Relativo (%%)  \n');
fprintf('-------------------------------------------------------------------------------------------------------------------------\n');

for indice_tabela = 1:numero_termos_expansao_n
    if vetor_termos_individuais(indice_tabela) >= 0
        caractere_sinal = '+';
    else
        caractere_sinal = '-';
    end

    fprintf(' %2d  |  %2d  |   %c   | %24.15e |  %.15f | %19.12e | %19.12e %%\n', ...
            vetor_indices_n(indice_tabela), ...
            vetor_graus_polinomio(indice_tabela), ...
            caractere_sinal, ...
            abs(vetor_termos_individuais(indice_tabela)), ...
            vetor_somas_acumuladas(indice_tabela), ...
            vetor_erros_absolutos(indice_tabela), ...
            vetor_erros_relativos_percentuais(indice_tabela));
end
fprintf('=========================================================================================================================\n');
fprintf(' RESUMO FINAL:\n');
fprintf(' Aproximacao obtida com %d termos: %.15f\n', numero_termos_expansao_n, somatorio_acumulado_serie);
fprintf(' Valor de referencia da funcao cos(x): %.15f\n', valor_referencia_cos_embutido);
fprintf(' Erro Absoluto Final: %.15e\n', vetor_erros_absolutos(end));
fprintf(' Erro Relativo Percentual Final: %.15e %%\n', vetor_erros_relativos_percentuais(end));
fprintf('=========================================================================================================================\n');

% -------------------------------------------------------------------------
% ETAPA 6: VISUALIZACAO GRAFICA DA TAXA DE CONVERGENCIA DA SERIE
% -------------------------------------------------------------------------
figura_taylor_cosseno = figure('Name', 'Convergencia da Serie de Taylor para cos(x)', 'NumberTitle', 'off');

% Subplot superior: Somas parciais se aproximando do valor de referencia
subplot(2, 1, 1);
plot(vetor_graus_polinomio, vetor_somas_acumuladas, 'b-o', 'LineWidth', 1.8, 'MarkerSize', 6, 'MarkerFaceColor', 'b');
hold on;
% Linha horizontal com o valor exato
plot([vetor_graus_polinomio(1), vetor_graus_polinomio(end)], ...
     [valor_referencia_cos_embutido, valor_referencia_cos_embutido], ...
     'r--', 'LineWidth', 1.5);
grid on;
title(sprintf('Convergencia das Somas Parciais para cos(%.4f rad)', angulo_entrada_x_radianos), 'FontSize', 11, 'FontWeight', 'bold');
xlabel('Grau do Polinomio de Taylor (2n)', 'FontSize', 10);
ylabel('Valor da Soma Parcial', 'FontSize', 10);
legend('Somas Parciais de Maclaurin', 'Valor Exato cos(x)', 'Location', 'SouthEast');
hold off;

% Subplot inferior: Erro Relativo Percentual em escala semilogaritmica
subplot(2, 1, 2);
% Substitui erros identicamente nulos por eps para permitir plot em escala logaritmica
vetor_erros_ajustados_log = max(vetor_erros_relativos_percentuais, eps);
semilogy(vetor_graus_polinomio, vetor_erros_ajustados_log, 'r-s', 'LineWidth', 1.8, 'MarkerSize', 6, 'MarkerFaceColor', 'r');
grid on;
title('Decaimento Super-Exponencial do Erro Relativo Percentual', 'FontSize', 11, 'FontWeight', 'bold');
xlabel('Grau do Polinomio de Taylor (2n)', 'FontSize', 10);
ylabel('Erro Relativo Percentual (%) [Escala Log]', 'FontSize', 10);
legend('Erro Relativo (%)', 'Location', 'NorthEast');

% =========================================================================
% FIM DO SCRIPT taylor_cosseno.m
% =========================================================================
