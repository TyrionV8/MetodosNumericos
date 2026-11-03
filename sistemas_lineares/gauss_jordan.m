% =========================================================================
% METODOS NUMERICOS - SISTEMAS DE EQUACOES LINEARES
% METODO DE GAUSS-JORDAN COM PIVOTEAMENTO PARCIAL E NORMALIZACAO UNITARIA
% =========================================================================
%
% 1. O QUE CALCULA:
%    Resolve sistemas de equacoes lineares algebricas da forma matricial:
%                           A * x = b
%    onde:
%      - matriz_coeficientes_A: matriz quadrada (ordem n x n) de coeficientes.
%      - vetor_termos_independentes_b: vetor de termos constantes (n x 1).
%      - vetor_solucao_final_x: vetor contendo as solucoes exatas (n x 1).
%
%    Formulacao Matematica e Principio de Gauss-Jordan:
%    a) Diferenca em Relacao a Eliminacao de Gauss Convencional:
%       Enquanto o metodo de Gauss classico apenas transforma a matriz em
%       triangular superior exigindo retrossubstituicao regressiva posterior,
%       o metodo de Gauss-Jordan zera sistematicamente os coeficientes
%       ABAIXO e ACIMA da diagonal principal. Apos a normalizacao de cada
%       linha pelo seu pivo, a matriz de coeficientes e reduzida a matriz
%       identidade I_n, permitindo a leitura imediata da solucao no vetor
%       aumentado final:
%             [ A | b ]  --->  [ I_n | x ]
%
%    b) Correcao Critica do Codigo Legado (Pivoteamento Parcial Real):
%       O codigo de referencia anterior apresentava uma tentativa defeituosa
%       de contornar pivos nulos, trocando indevidamente a primeira linha e
%       indexando variaveis inexistentes. Nesta implementacao corrigida:
%       - Realiza-se o pivoteamento parcial rigoroso: na coluna k, busca-se
%         o elemento de maior modulo entre todas as linhas candidatas de k ate n:
%               |a_lk| = max_{k <= p <= n} |a_pk|
%       - Permuta-se integralmente a linha k com a linha l na matriz aumentada.
%       - Isso elimina o risco de divisao por zero e impede o crescimento
%         descontrolado dos multiplicadores, garantindo estabilidade numerica.
%
%    c) Fases do Procedimento de Diagonalizacao:
%       - Fase 1 (Eliminacao Descendente): zera os elementos abaixo da diagonal,
%         construindo a forma triangular superior.
%       - Fase 2 (Eliminacao Ascendente): zera os elementos acima da diagonal,
%         de baixo para cima, isolando os termos diagonais.
%       - Fase 3 (Normalizacao da Diagonal): divide cada linha pelo seu pivo
%         diagonal, transformando a diagonal principal em 1.
%       - Fase 4 (Entrega Direta): a ultima coluna aumentada e a solucao x.
%
% 2. O QUE ENTREGA:
%    - Exibicao da matriz aumentada inicial [A | b].
%    - Registro de cada operacao de troca de linhas realizada no pivoteamento.
%    - Matriz aumentada apos eliminacao descendente (triangular superior).
%    - Matriz aumentada completamente diagonalizada e normalizada [I | x].
%    - Vetor solucao final x com identificacao das variaveis fisicas (a, R, T).
%    - Norma euclidiana do residuo ||r||_2 = ||A*x - b||_2.
%
% 3. O QUE PRECISA TER (Pre-requisitos):
%    - Matriz de coeficientes A de dimensao n x n nao-singular (det(A) ~= 0).
%    - Vetor de termos independentes b com dimensao n compativel.
%    - Ambiente GNU Octave ou MATLAB.
%
% 4. EXEMPLO DE USO E RESULTADO ESPERADO:
%    Sistema representativo de calibracao experimental de parametros:
%      [ 55   0   1 ] [ a ]   [ 499.55 ]
%      [ 75   1  -1 ] [ R ] = [ 687.75 ]
%      [ 60  -1   0 ] [ T ]   [ 528.60 ]
%
%    Resultado esperado:
%      a =  9.031053 (aprox.  9.0311)
%      R = 13.263158 (aprox. 13.2632)
%      T =  2.842105 (aprox.  2.8421)
%      Norma Euclidiana do Residuo < 1.000000e-12
% =========================================================================

% -------------------------------------------------------------------------
% ETAPA 0: LIMPEZA DO AMBIENTE DE EXECUCAO
% -------------------------------------------------------------------------
% Limpa a tela do console para exibicao limpa e organizada
clc;
% Remove todas as variaveis da memoria para evitar interferencia de execucoes anteriores
clear all;
% Fecha todas as figuras graficas abertas
close all;

% -------------------------------------------------------------------------
% ETAPA 1: DEFINICAO DOS PARAMETROS DE ENTRADA (MODIFIQUE AQUI)
% -------------------------------------------------------------------------
% Matriz de coeficientes do sistema linear A (dimensao n x n)
matriz_coeficientes_A = [ ...
    55 ,  0 ,  1 ; ...
    75 ,  1 , -1 ; ...
    60 , -1 ,  0   ...
];

% Vetor de termos independentes b
vetor_termos_independentes_b = [499.55 ; 687.75 ; 528.60];

% Tolerancia numerica para verificacao de singularidade e deteccao de pivo nulo
tolerancia_singularidade_pivo = 1.0e-12;

% -------------------------------------------------------------------------
% ETAPA 2: VALIDACAO DIMENSIONAL E MONTAGEM DA MATRIZ AUMENTADA
% -------------------------------------------------------------------------
% Captura as dimensoes da matriz de coeficientes
[quantidade_linhas_matriz_A, quantidade_colunas_matriz_A] = size(matriz_coeficientes_A);

% Valida se a matriz de coeficientes e perfeitamente quadrada
if quantidade_linhas_matriz_A ~= quantidade_colunas_matriz_A
    error('Erro dimensional: A matriz de coeficientes A deve ser quadrada (n x n).');
end

% Define a ordem do sistema linear
ordem_sistema_linear = quantidade_linhas_matriz_A;

% Converte o vetor de termos independentes em vetor coluna, caso fornecido como linha
if isrow(vetor_termos_independentes_b)
    vetor_termos_independentes_b = vetor_termos_independentes_b(:);
end

% Valida se o numero de termos independentes coincide com o numero de equacoes
if length(vetor_termos_independentes_b) ~= ordem_sistema_linear
    error('Erro dimensional: O vetor b deve possuir a mesma quantidade de elementos que as linhas de A.');
end

% Concatena a matriz de coeficientes com o vetor b, gerando a matriz aumentada [A | b]
matriz_aumentada_sistema = [matriz_coeficientes_A, vetor_termos_independentes_b];

% Quantidade total de colunas da matriz aumentada (n + 1)
quantidade_colunas_matriz_aumentada = ordem_sistema_linear + 1;

% Exibe o estado inicial da matriz aumentada no terminal
disp('=========================================================================');
disp('   METODO DE GAUSS-JORDAN COM PIVOTEAMENTO PARCIAL E DIAGONALIZACAO');
disp('=========================================================================');
disp('Matriz Aumentada Inicial [A | b]:');
disp(matriz_aumentada_sistema);

% -------------------------------------------------------------------------
% ETAPA 3: FASE 1 - ELIMINACAO DESCENDENTE COM PIVOTEAMENTO PARCIAL REAL
% -------------------------------------------------------------------------
% Zera todos os elementos abaixo da diagonal principal (triangular superior)
for coluna_pivo_atual = 1 : (ordem_sistema_linear - 1)

    % Extrai os coeficientes candidatos a pivo na coluna atual da diagonal para baixo
    vetor_subcoluna_candidata_pivo = matriz_aumentada_sistema(coluna_pivo_atual : ordem_sistema_linear, coluna_pivo_atual);

    % Identifica o maior elemento em modulo para garantir estabilidade numerica
    [maior_valor_absoluto_pivo, posicao_relativa_maior_pivo] = max(abs(vetor_subcoluna_candidata_pivo));

    % Converte o indice relativo para o indice global da linha do maior pivo
    linha_maior_pivo_encontrada = (coluna_pivo_atual - 1) + posicao_relativa_maior_pivo;

    % Verifica se a coluna e identicamente nula, indicando singularidade
    if maior_valor_absoluto_pivo < tolerancia_singularidade_pivo
        error('Matriz singular: Pivo nulo detectado na coluna %d. Nao ha solucao unica.', coluna_pivo_atual);
    end

    % Se a linha de maior pivo nao for a linha diagonal atual, realiza a permutacao elementar
    if linha_maior_pivo_encontrada ~= coluna_pivo_atual
        fprintf('-> Pivoteamento Parcial: Permutando Linha %d com Linha %d (Pivo em modulo = %.6f)\n', ...
            coluna_pivo_atual, linha_maior_pivo_encontrada, maior_valor_absoluto_pivo);

        % Variavel auxiliar para realizar a troca segura das duas linhas inteiras
        linha_auxiliar_troca_temporaria = matriz_aumentada_sistema(coluna_pivo_atual, :);

        % Linha atual recebe a linha contendo o maior elemento
        matriz_aumentada_sistema(coluna_pivo_atual, :) = matriz_aumentada_sistema(linha_maior_pivo_encontrada, :);

        % Linha onde estava o maior pivo recebe a linha anterior
        matriz_aumentada_sistema(linha_maior_pivo_encontrada, :) = linha_auxiliar_troca_temporaria;
    end

    % Elimina os coeficientes de todas as linhas abaixo da linha pivo
    for linha_para_eliminar = (coluna_pivo_atual + 1) : ordem_sistema_linear

        % Calcula o multiplicador de eliminacao m_ik = a_ik / a_kk
        fator_multiplicador_eliminacao = matriz_aumentada_sistema(linha_para_eliminar, coluna_pivo_atual) / ...
                                         matriz_aumentada_sistema(coluna_pivo_atual, coluna_pivo_atual);

        % Operacao elementar: Linha_i <- Linha_i - m_ik * Linha_k
        matriz_aumentada_sistema(linha_para_eliminar, :) = matriz_aumentada_sistema(linha_para_eliminar, :) - ...
            fator_multiplicador_eliminacao * matriz_aumentada_sistema(coluna_pivo_atual, :);

        % Forca o zero analitico para evitar residuos residuais de ponto flutuante
        matriz_aumentada_sistema(linha_para_eliminar, coluna_pivo_atual) = 0.0;
    end
end

% Verificacao de singularidade no ultimo pivo diagonal
if abs(matriz_aumentada_sistema(ordem_sistema_linear, ordem_sistema_linear)) < tolerancia_singularidade_pivo
    error('Matriz singular: Ultimo elemento diagonal e nulo. O sistema nao admite solucao unica.');
end

disp(' ');
disp('Matriz Apos Eliminacao Descendente (Triangular Superior):');
disp(matriz_aumentada_sistema);

% -------------------------------------------------------------------------
% ETAPA 4: FASE 2 - ELIMINACAO ASCENDENTE (ZERAR TERMOS ACIMA DA DIAGONAL)
% -------------------------------------------------------------------------
% Itera de baixo para cima, da ultima coluna ate a segunda
for linha_ascendente_pivo = ordem_sistema_linear : -1 : 2

    % Itera em todas as linhas situadas acima da linha pivo
    for linha_acima_para_zerar = (linha_ascendente_pivo - 1) : -1 : 1

        % Verificacao defensiva de seguranca contra pivo nulo na eliminacao ascendente
        if abs(matriz_aumentada_sistema(linha_ascendente_pivo, linha_ascendente_pivo)) < tolerancia_singularidade_pivo
            error('Divisao por zero na eliminacao ascendente: pivo nulo na linha %d.', linha_ascendente_pivo);
        end

        % Calcula o multiplicador ascendente m_ik = a_ik / a_kk
        fator_multiplicador_ascendente = matriz_aumentada_sistema(linha_acima_para_zerar, linha_ascendente_pivo) / ...
                                         matriz_aumentada_sistema(linha_ascendente_pivo, linha_ascendente_pivo);

        % Operacao elementar: Linha_acima <- Linha_acima - m_ik * Linha_pivo
        matriz_aumentada_sistema(linha_acima_para_zerar, :) = matriz_aumentada_sistema(linha_acima_para_zerar, :) - ...
            fator_multiplicador_ascendente * matriz_aumentada_sistema(linha_ascendente_pivo, :);

        % Forca o elemento acima da diagonal a zero exato
        matriz_aumentada_sistema(linha_acima_para_zerar, linha_ascendente_pivo) = 0.0;
    end
end

disp(' ');
disp('Matriz Apos Eliminacao Ascendente (Diagonal Pura):');
disp(matriz_aumentada_sistema);

% -------------------------------------------------------------------------
% ETAPA 5: FASE 3 - NORMALIZACAO DA DIAGONAL PRINCIPAL (FORMA REDUZIDA [I | x])
% -------------------------------------------------------------------------
% Divide cada linha pelo seu elemento diagonal principal para obter pivos unitarios
for linha_normalizacao_atual = 1 : ordem_sistema_linear

    % Guarda o valor do pivo diagonal antes da divisao da linha
    elemento_pivo_diagonal_normalizacao = matriz_aumentada_sistema(linha_normalizacao_atual, linha_normalizacao_atual);

    % Verificacao defensiva de seguranca contra pivo nulo na normalizacao
    if abs(elemento_pivo_diagonal_normalizacao) < tolerancia_singularidade_pivo
        error('Divisao por zero na normalizacao: pivo nulo na linha %d.', linha_normalizacao_atual);
    end

    % Aplica a operacao elementar: Linha_k <- Linha_k / a_kk
    matriz_aumentada_sistema(linha_normalizacao_atual, :) = ...
        matriz_aumentada_sistema(linha_normalizacao_atual, :) / elemento_pivo_diagonal_normalizacao;

    % Garante que a diagonal principal seja exatamente 1.0
    matriz_aumentada_sistema(linha_normalizacao_atual, linha_normalizacao_atual) = 1.0;
end

disp(' ');
disp('Matriz Aumentada Final Diagonalizada e Normalizada [I_n | x]:');
disp(matriz_aumentada_sistema);

% -------------------------------------------------------------------------
% ETAPA 6: ENTREGA DIRETA DA SOLUCAO E ANALISE DE EXATIDAO RESIDUAL
% -------------------------------------------------------------------------
% No metodo de Gauss-Jordan, a ultima coluna da matriz aumentada ja e a solucao direta
vetor_solucao_final_x = matriz_aumentada_sistema(:, quantidade_colunas_matriz_aumentada);

disp('-------------------------------------------------------------------------');
disp('Vetor Solucao Final Obtido (x):');
for indice_exibicao_solucao = 1 : ordem_sistema_linear
    fprintf('  x(%d) = %14.6f\n', indice_exibicao_solucao, vetor_solucao_final_x(indice_exibicao_solucao));
end

% Interpretacao fisica dos coeficientes do problema demonstrativo de calibracao
if ordem_sistema_linear == 3
    coeficiente_calibracao_a = vetor_solucao_final_x(1);
    resistencia_eletrica_R = vetor_solucao_final_x(2);
    temperatura_referencia_T = vetor_solucao_final_x(3);
    fprintf('\nInterpretacao dos Parametros Fisicos do Sistema:\n');
    fprintf('  Coeficiente de calibracao a = %12.6f\n', coeficiente_calibracao_a);
    fprintf('  Resistencia eletrica R      = %12.6f\n', resistencia_eletrica_R);
    fprintf('  Temperatura de referencia T = %12.6f\n', temperatura_referencia_T);
end

% Calculo do vetor residual para confirmacao analitica da exatidao
vetor_residuo_linear = matriz_coeficientes_A * vetor_solucao_final_x - vetor_termos_independentes_b;

% Norma euclidiana do residuo ||r||_2
norma_euclidiana_residuo = norm(vetor_residuo_linear, 2);

fprintf('\nNorma Euclidiana do Residuo ||A*x - b||_2: %e\n', norma_euclidiana_residuo);
disp('=========================================================================');
