% =========================================================================
% METODOS NUMERICOS - CALCULO DE RAIZES DE EQUACOES NAO LINEARES
% METODO DE NEWTON-RAPHSON (METODO DAS TANGENTES)
% =========================================================================
% 
% 1. O QUE CALCULA:
%    Determina a raiz real de uma equacao nao linear f(x) = 0 a partir de uma
%    estimativa inicial x0, utilizando a derivada analitica f'(x).
%    O metodo baseia-se na expansao em Serie de Taylor de primeira ordem da
%    funcao f(x) em torno do ponto corrente x_k:
%        f(x) ~= f(x_k) + f'(x_k) * (x - x_k) = 0
%    Isolando x, obtem-se a equacao de recorrencia de Newton-Raphson:
%        x_(k+1) = x_k - f(x_k) / f'(x_k)
%    Geometricamente, x_(k+1) e o ponto onde a reta tangente a curva y = f(x)
%    em (x_k, f(x_k)) intercepta o eixo das abscissas (y = 0).
%    Sob condicoes favoraveis (raiz simples e estimativa inicial proxima),
%    a convergencia do metodo e quadratica.
%
% 2. O QUE ENTREGA:
%    - Tabela passo a passo com: Iteracao, x_k, f(x_k), f'(x_k), x_(k+1), Erro Absoluto e Erro Relativo.
%    - Raiz aproximada final ('raiz_aproximada_final').
%    - Numero total de iteracoes executadas.
%    - Mensagens detalhadas de convergencia ou indicacao de possivel divergencia/falha.
%
% 3. CORRECOES APLICADAS EM RELACAO A VERSAO LEGADA:
%    [BUG CRITICO 1 CORRIGIDO]: A variavel da solucao 'xz' era criada apenas
%    dentro da clausula if quando havia convergencia. Caso o numero maximo de
%    ciclos fosse atingido sem convergir, 'xz' nao existia e o comando fprintf
%    gerava erro de execucao ('xz' undefined). Nesta versao, 'raiz_aproximada_final'
%    e inicializada com NaN no topo.
%    [BUG CRITICO 2 CORRIGIDO]: O erro relativo na raiz do repositorio usava a
%    formula 'abs(xi-x0)/x0' sem modulo no denominador, gerando erros negativos
%    quando x0 < 0. Corrigido para 'abs((proxima - atual) / atual)' com protecao
%    contra divisao por zero caso 'atual' seja nulo.
%    [ROBUSTEZ ADICIONADA]: Verificacao explicita de derivada nula (|f'(x)| < eps)
%    para evitar divisao por zero e parada brusca quando a tangente e horizontal.
%
% 4. PARAMETROS E PRE-REQUISITOS:
%    - Funcao f(x) diferenciavel no intervalo que contem a raiz.
%    - Derivada analitica f'(x) continua e nao nula na raiz.
%    - Estimativa inicial suficientemente proxima da raiz desejada.
%    - Tolerancia de convergencia e numero maximo de iteracoes permitidas.
%
% 5. EXEMPLO DE USO E RESULTADO ESPERADO:
%    Funcao: f(x) = x^3 + 3.8*x^2 - 8.6*x - 24.4
%    Derivada: f'(x) = 3*x^2 + 7.6*x - 8.6
%    Estimativa inicial x0 = 4.0
%    Tolerancia: 1e-6
%    Resultado esperado: raiz aproximada x ~= 2.706938 (obtida em 5 iteracoes).
% =========================================================================

% Limpeza inicial do ambiente de trabalho do Octave / MATLAB
clc;        % Limpa a janela de comandos
clear all;  % Limpa todas as variaveis da memoria do workspace
close all;  % Fecha todas as janelas graficas abertas

%% =========================================================================
%% BLOCO DE PARAMETROS DE ENTRADA DO USUARIO
%% =========================================================================

% Definicao da funcao matematica objetivo f(x)
funcao_objetivo = @(x)(x.^3 + 3.8*x.^2 - 8.6*x - 24.4);

% Definicao analitica da primeira derivada da funcao objetivo f'(x)
derivada_funcao_objetivo = @(x)(3.0*x.^2 + 7.6*x - 8.6);

% Estimativa inicial para o ponto de partida do metodo iterativo
estimativa_inicial_x = 4.0;

% Tolerancia para o criterio de convergencia (erro relativo e residuo)
tolerancia_convergencia = 1e-6;

% Numero maximo de iteracoes admitidas antes de abortar por seguranca
numero_maximo_iteracoes = 50;

%% =========================================================================
%% BLOCO DE PROCESSAMENTO E CALCULO NUMERICO
%% =========================================================================

% Inicializacao da variavel que armazena a raiz com NaN para seguranca
raiz_aproximada_final = NaN;

% Variavel booleana de estado para registrar se houve convergencia
convergencia_atingida = false;

% Contador para o total de iteracoes executadas
total_iteracoes_executadas = 0;

% Inicializacao da estimativa de busca com o chute inicial fornecido
estimativa_atual_x = estimativa_inicial_x;

% Inicializacao defensiva das variaveis para evitar erro em caso de saida prematura
proxima_estimativa_x = estimativa_inicial_x;
erro_absoluto_passo = Inf;
erro_relativo_passo = Inf;

% Impressao do cabecalho explicativo da tabela de iteracoes
fprintf('\n=========================================================================================================================\n');
fprintf('                           METODO DE NEWTON-RAPHSON - CONVERGENCIA PASSO A PASSO                                         \n');
fprintf('=========================================================================================================================\n');
fprintf('  Iter |       x_atual     |       f(x)        |       f''(x)      |      x_proximo    |   Erro Absoluto   |   Erro Relativo\n');
fprintf('-------+-------------------+-------------------+-------------------+-------------------+-------------------+----------------\n');

% Laco principal de repeticao do metodo de Newton-Raphson
for iteracao_atual = 1:numero_maximo_iteracoes
    total_iteracoes_executadas = iteracao_atual;

    % Avaliacao da funcao objetivo no ponto atual
    valor_funcao_atual = funcao_objetivo(estimativa_atual_x);

    % Avaliacao da derivada da funcao objetivo no ponto atual
    valor_derivada_atual = derivada_funcao_objetivo(estimativa_atual_x);

    % Verificacao de seguranca contra derivada nula (tangente horizontal)
    if abs(valor_derivada_atual) < eps
        fprintf('\n[ERRO NUMERICO]: A derivada f''(x) e nula no ponto x = %.8f.\n', estimativa_atual_x);
        fprintf('A tangente e horizontal e nao intercepta o eixo x. O metodo nao pode prosseguir.\n');
        break;
    end

    % Calculo do passo de correcao de Newton: Delta x = f(x) / f'(x)
    passo_correcao_newton = valor_funcao_atual / valor_derivada_atual;

    % Calculo da proxima estimativa da raiz: x_(k+1) = x_k - Delta x
    proxima_estimativa_x = estimativa_atual_x - passo_correcao_newton;

    % Calculo do deslocamento absoluto entre passos sucessivos
    erro_absoluto_passo = abs(proxima_estimativa_x - estimativa_atual_x);

    % Calculo robusto do erro relativo com protecao contra raiz nula ou negativa
    if abs(estimativa_atual_x) > eps
        erro_relativo_passo = abs(erro_absoluto_passo / estimativa_atual_x);
    else
        erro_relativo_passo = erro_absoluto_passo;
    end

    % Avaliacao da funcao objetivo na nova aproximacao gerada
    valor_funcao_proxima = funcao_objetivo(proxima_estimativa_x);

    % Exibicao dos resultados formatados da iteracao corrente
    fprintf('  %4d |   %15.8f |   %15.8e |   %15.8f |   %15.8f |   %17.8e |   %14.8e\n', ...
            iteracao_atual, ...
            estimativa_atual_x, ...
            valor_funcao_atual, ...
            valor_derivada_atual, ...
            proxima_estimativa_x, ...
            erro_absoluto_passo, ...
            erro_relativo_passo);

    % Teste de parada: verifica residuo da funcao OU erro relativo OU erro absoluto
    if (abs(valor_funcao_proxima) < tolerancia_convergencia) || ...
       (erro_relativo_passo < tolerancia_convergencia) || ...
       (erro_absoluto_passo < tolerancia_convergencia)
        raiz_aproximada_final = proxima_estimativa_x;
        convergencia_atingida = true;
        break; % Interrompe o laco pois a tolerancia requerida foi satisfeita
    end

    % Atualiza a estimativa para o proximo ciclo de iteracao
    estimativa_atual_x = proxima_estimativa_x;
end

fprintf('-------------------------------------------------------------------------------------------------------------------------\n');

% Exibicao do relatorio final de encerramento do script
if convergencia_atingida
    fprintf('\n>>> SUCESSO: O Metodo de Newton-Raphson convergiu satisfatoriamente!\n');
    fprintf('    Raiz aproximada encontrada : %.10f\n', raiz_aproximada_final);
    fprintf('    Residuo f(raiz)            : %.10e\n', funcao_objetivo(raiz_aproximada_final));
    fprintf('    Ultimo erro absoluto       : %.10e (Tolerancia: %.1e)\n', erro_absoluto_passo, tolerancia_convergencia);
    fprintf('    Ultimo erro relativo       : %.10e\n', erro_relativo_passo);
    fprintf('    Total de iteracoes         : %d\n\n', total_iteracoes_executadas);
else
    fprintf('\n>>> AVISO: O metodo nao convergiu dentro do limite de %d iteracoes.\n', numero_maximo_iteracoes);
    fprintf('    Ultima aproximacao obtida  : %.10f\n', proxima_estimativa_x);
    fprintf('    Erro absoluto alcancado    : %.10e (Tolerancia: %.1e)\n', erro_absoluto_passo, tolerancia_convergencia);
    fprintf('    Erro relativo alcancado    : %.10e\n\n', erro_relativo_passo);
end
