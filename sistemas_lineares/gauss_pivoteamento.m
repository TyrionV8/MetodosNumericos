% =========================================================================
% METODOS NUMERICOS - SISTEMAS DE EQUACOES LINEARES
% ELIMINACAO GAUSSIANA COM PIVOTEAMENTO PARCIAL E RETROSSUBSTITUICAO
% =========================================================================
%
% 1. O QUE CALCULA:
%    Resolve sistemas lineares determinados da forma matricial:
%                           A * x = b
%    onde:
%      - matriz_coeficientes_A: matriz quadrada (ordem n x n) de coeficientes.
%      - vetor_termos_independentes_b: vetor de termos constantes (n x 1).
%      - vetor_solucao_final_x: vetor com os valores das incognitas (n x 1).
%
%    Formulacao Matematica e Teoremas Envolvidos:
%    a) Escalonamento por Operacoes Elementares:
%       O Teorema da Equivalencia de Sistemas Lineares garante que multiplicar
%       uma linha por uma constante nao nula ou somar a uma linha o multiplo
%       de outra preserva o conjunto solucao do sistema original.
%       A operacao elementar de eliminacao utilizada e:
%             L_i <- L_i - m_ik * L_k
%       onde m_ik = a_ik / a_kk e o fator multiplicador que zera o elemento
%       abaixo do elemento pivo a_kk na coluna k.
%
%    b) Estabilidade Numerica e Pivoteamento Parcial:
%       Em aritmetica de ponto flutuante (padrao IEEE 754), se o pivo |a_kk|
%       for nulo, ocorre divisao por zero imediata. Se for muito proximo de
%       zero, o multiplicador m_ik torna-se muito grande, gerando erros
%       catastroficos de arredondamento e cancelamento subtrativo.
%       O pivoteamento parcial resolve isso inspecionando a coluna k de k ate n:
%             |a_lk| = max_{k <= p <= n} |a_pk|
%       e permutando as linhas k e l (L_k <-> L_l). Isso forca que todos os
%       multiplicadores satisfacam |m_ik| <= 1, garantindo estabilidade.
%
%    c) Retrossubstituicao Regressiva:
%       Uma vez obtida a matriz triangular superior U, resolve-se o sistema
%       de tras para frente (da ultima incognita x_n ate a primeira x_1):
%             x_n = c_n / u_nn
%             x_i = ( c_i - sum_{j=i+1}^n u_ij * x_j ) / u_ii
%
% 2. O QUE ENTREGA:
%    - Exibicao da matriz aumentada inicial [A | b].
%    - Notificacoes didaticas no console a cada troca de linhas efetuada.
%    - Matriz aumentada escalonada em forma triangular superior [U | c].
%    - Vetor solucao final x contendo o valor exato de cada variavel.
%    - Analise residual: vetor de residuo r = A*x - b e norma euclidiana ||r||_2.
%
% 3. O QUE PRECISA TER (Pre-requisitos):
%    - Matriz de coeficientes A de ordem n x n nao-singular (det(A) ~= 0).
%    - Vetor de termos independentes b com tamanho compativel n.
%    - Ambiente GNU Octave ou MATLAB.
%
% 4. EXEMPLO DE USO E RESULTADO ESPERADO:
%    Para o sistema demonstrativo fornecido:
%      [ 4  3  2 ] [ x1 ]   [ 960 ]
%      [ 1  3  1 ] [ x2 ] = [ 510 ]
%      [ 2  1  3 ] [ x3 ]   [ 610 ]
%
%    Resultado esperado:
%      x1 = 120.000000
%      x2 = 100.000000
%      x3 =  90.000000
%      Norma Euclidiana do Residuo = 0.000000e+00
% =========================================================================

% -------------------------------------------------------------------------
% ETAPA 0: LIMPEZA DO AMBIENTE DE EXECUCAO
% -------------------------------------------------------------------------
% Limpa a janela de comandos do console para melhor visualizacao dos resultados
clc;
% Elimina todas as variaveis do espaco de trabalho para evitar conflitos residuais
clear all;
% Fecha eventuais janelas graficas que estejam abertas
close all;

% -------------------------------------------------------------------------
% ETAPA 1: DEFINICAO DOS PARAMETROS DE ENTRADA (MODIFIQUE AQUI)
% -------------------------------------------------------------------------
% Matriz de coeficientes do sistema linear A * x = b (dimensao n x n)
matriz_coeficientes_A = [ ...
    4 , 3 , 2 ; ...
    1 , 3 , 1 ; ...
    2 , 1 , 3   ...
];

% Vetor de termos independentes (pode ser definido como vetor linha ou coluna)
vetor_termos_independentes_b = [960 ; 510 ; 610];

% Tolerancia numerica para deteccao de matriz singular ou pivo quase nulo
tolerancia_singularidade_pivo = 1.0e-12;

% -------------------------------------------------------------------------
% ETAPA 2: VALIDACAO DIMENSIONAL E PREPARACAO DA MATRIZ AUMENTADA
% -------------------------------------------------------------------------
% Obtencao das dimensoes da matriz de coeficientes
[quantidade_linhas_matriz_A, quantidade_colunas_matriz_A] = size(matriz_coeficientes_A);

% Verificacao se a matriz de coeficientes e rigorosamente quadrada
if quantidade_linhas_matriz_A ~= quantidade_colunas_matriz_A
    error('Erro dimensional: A matriz de coeficientes deve ser quadrada (n x n).');
end

% A ordem do sistema linear corresponde a quantidade de linhas da matriz quadrada
ordem_sistema_linear = quantidade_linhas_matriz_A;

% Adequacao do vetor de termos independentes para formato de vetor coluna (n x 1)
if isrow(vetor_termos_independentes_b)
    vetor_termos_independentes_b = vetor_termos_independentes_b(:);
end

% Verificacao da compatibilidade dimensional entre a matriz A e o vetor b
if length(vetor_termos_independentes_b) ~= ordem_sistema_linear
    error('Erro dimensional: O vetor de termos independentes b deve ter dimensao igual a ordem de A.');
end

% Construcao da matriz aumentada do sistema linear [A | b] para operacoes combinadas
matriz_aumentada_sistema = [matriz_coeficientes_A, vetor_termos_independentes_b];

% Determinacao da quantidade de colunas da matriz aumentada (igual a n + 1)
quantidade_colunas_matriz_aumentada = ordem_sistema_linear + 1;

% Exibicao explicativa da matriz aumentada no estado inicial
disp('=========================================================================');
disp('   METODO DE ELIMINACAO GAUSSIANA COM PIVOTEAMENTO PARCIAL');
disp('=========================================================================');
disp('Matriz Aumentada Inicial [A | b]:');
disp(matriz_aumentada_sistema);

% -------------------------------------------------------------------------
% ETAPA 3: ESCALONAMENTO DIRETO COM PIVOTEAMENTO PARCIAL REAL
% -------------------------------------------------------------------------
% Itera atraves das colunas de 1 ate (n - 1) para triangularizar a matriz
for coluna_pivo_atual = 1 : (ordem_sistema_linear - 1)

    % Extrai a subcoluna a partir da diagonal principal ate a ultima linha
    vetor_subcoluna_candidata_pivo = matriz_aumentada_sistema(coluna_pivo_atual : ordem_sistema_linear, coluna_pivo_atual);

    % Localiza o elemento com maior valor absoluto na coluna para ser o novo pivo
    [maior_valor_absoluto_pivo, posicao_relativa_maior_pivo] = max(abs(vetor_subcoluna_candidata_pivo));

    % Calcula o indice real da linha onde se encontra o maior pivo
    linha_maior_pivo_encontrada = (coluna_pivo_atual - 1) + posicao_relativa_maior_pivo;

    % Verifica se o maior pivo encontrado e nulo ou proximo de zero dentro da tolerancia
    if maior_valor_absoluto_pivo < tolerancia_singularidade_pivo
        error('Matriz singular ou quase-singular: Pivo nulo detectado na coluna %d. O sistema nao admite solucao unica.', coluna_pivo_atual);
    end

    % Aplica a permutacao de linhas se a linha do maior pivo nao for a linha pivô atual
    if linha_maior_pivo_encontrada ~= coluna_pivo_atual
        % Registra didaticamente a operacao elementar de troca de linhas realizada
        fprintf('-> Pivoteamento Parcial: Permutando Linha %d com Linha %d (Maior pivo em modulo = %.6f)\n', ...
            coluna_pivo_atual, linha_maior_pivo_encontrada, maior_valor_absoluto_pivo);

        % Salva temporariamente a linha pivo atual para troca
        linha_auxiliar_troca_temporaria = matriz_aumentada_sistema(coluna_pivo_atual, :);

        % Substitui a linha pivo atual pela linha que contem o maior pivo
        matriz_aumentada_sistema(coluna_pivo_atual, :) = matriz_aumentada_sistema(linha_maior_pivo_encontrada, :);

        % Atribui a linha salva a posicao da linha permutada
        matriz_aumentada_sistema(linha_maior_pivo_encontrada, :) = linha_auxiliar_troca_temporaria;
    end

    % Executa a eliminacao gaussiana em todas as linhas abaixo da linha pivo
    for linha_para_eliminar = (coluna_pivo_atual + 1) : ordem_sistema_linear

        % Calcula o fator multiplicador m_ik = a_ik / a_kk
        fator_multiplicador_eliminacao = matriz_aumentada_sistema(linha_para_eliminar, coluna_pivo_atual) / ...
                                         matriz_aumentada_sistema(coluna_pivo_atual, coluna_pivo_atual);

        % Aplica a operacao elementar: Linha_i <- Linha_i - m_ik * Linha_k
        matriz_aumentada_sistema(linha_para_eliminar, :) = matriz_aumentada_sistema(linha_para_eliminar, :) - ...
            fator_multiplicador_eliminacao * matriz_aumentada_sistema(coluna_pivo_atual, :);

        % Forca matematicamente o elemento eliminado a zero para eliminar ruidos de ponto flutuante
        matriz_aumentada_sistema(linha_para_eliminar, coluna_pivo_atual) = 0.0;
    end
end

% Verificacao de singularidade no ultimo elemento da diagonal principal (u_nn)
if abs(matriz_aumentada_sistema(ordem_sistema_linear, ordem_sistema_linear)) < tolerancia_singularidade_pivo
    error('Matriz singular detectada no ultimo elemento diagonal. O sistema nao admite solucao unica.');
end

% Exibicao didatica da matriz aumentada triangular superior resultante
disp(' ');
disp('Matriz Aumentada Escalonada [U | c] (Triangular Superior):');
disp(matriz_aumentada_sistema);

% -------------------------------------------------------------------------
% ETAPA 4: RETROSSUBSTITUICAO REGRESSIVA
% -------------------------------------------------------------------------
% Aloca espaco de memoria para o vetor solucao final
vetor_solucao_final_x = zeros(ordem_sistema_linear, 1);

% Calcula a ultima incognita x_n = c_n / u_nn diretamente da ultima equacao
vetor_solucao_final_x(ordem_sistema_linear) = ...
    matriz_aumentada_sistema(ordem_sistema_linear, quantidade_colunas_matriz_aumentada) / ...
    matriz_aumentada_sistema(ordem_sistema_linear, ordem_sistema_linear);

% Itera regressivamente da penultima linha (n - 1) ate a primeira linha (1)
for linha_regressiva_atual = (ordem_sistema_linear - 1) : -1 : 1

    % Inicializa o acumulador de soma para os produtos coeficientes * incognitas ja obtidas
    soma_produtos_retrossubstituicao = 0.0;

    % Soma a contribuicao de todas as incognitas ja determinadas a direita da diagonal
    for coluna_incognita_conhecida = (linha_regressiva_atual + 1) : ordem_sistema_linear
        soma_produtos_retrossubstituicao = soma_produtos_retrossubstituicao + ...
            matriz_aumentada_sistema(linha_regressiva_atual, coluna_incognita_conhecida) * ...
            vetor_solucao_final_x(coluna_incognita_conhecida);
    end

    % Verificacao de seguranca contra divisao por zero na retrossubstituicao
    if abs(matriz_aumentada_sistema(linha_regressiva_atual, linha_regressiva_atual)) < tolerancia_singularidade_pivo
        error('Divisao por zero na retrossubstituicao: elemento diagonal nulo na linha %d.', linha_regressiva_atual);
    end

    % Isola a incognita atual: x_i = (c_i - soma_produtos) / u_ii
    vetor_solucao_final_x(linha_regressiva_atual) = ...
        (matriz_aumentada_sistema(linha_regressiva_atual, quantidade_colunas_matriz_aumentada) - ...
         soma_produtos_retrossubstituicao) / ...
        matriz_aumentada_sistema(linha_regressiva_atual, linha_regressiva_atual);
end

% -------------------------------------------------------------------------
% ETAPA 5: EXIBICAO DOS RESULTADOS E AUDITORIA DE EXATIDAO RESIDUAL
% -------------------------------------------------------------------------
disp('-------------------------------------------------------------------------');
disp('Vetor Solucao Final Obtido (x):');
for indice_exibicao_solucao = 1 : ordem_sistema_linear
    fprintf('  x(%d) = %14.6f\n', indice_exibicao_solucao, vetor_solucao_final_x(indice_exibicao_solucao));
end

% Calculo do vetor residual r = A * x - b para atestar a precisao numerica
vetor_residuo_linear = matriz_coeficientes_A * vetor_solucao_final_x - vetor_termos_independentes_b;

% Calculo da norma euclidiana do residuo ||r||_2 = sqrt(sum(r_i^2))
norma_euclidiana_residuo = norm(vetor_residuo_linear, 2);

fprintf('\nNorma Euclidiana do Residuo ||A*x - b||_2: %e\n', norma_euclidiana_residuo);
disp('=========================================================================');
