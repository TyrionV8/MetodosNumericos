% =========================================================================
% METODOS NUMERICOS - CALCULO DE RAIZES DE EQUACOES NAO LINEARES
% METODO DA SECANTE
% =========================================================================
% 
% 1. O QUE CALCULA:
%    Determina uma aproximacao numerica para a raiz real da equacao f(x) = 0
%    a partir de duas estimativas iniciais (x0 e x1), sem requerer a avaliacao
%    analitica da derivada f'(x) e sem a obrigatoriedade de que o intervalo
%    inicial contenha uma troca de sinal.
%    O metodo aproxima a derivada de Newton-Raphson pelo quociente das diferencas
%    finitas entre os dois pontos mais recentes:
%        f'(x_k) ~= (f(x_k) - f(x_(k-1))) / (x_k - x_(k-1))
%    Substituindo na formula de Newton, obtem-se a recorrencia da Secante:
%        x_(k+1) = x_k - f(x_k) * (x_k - x_(k-1)) / (f(x_k) - f(x_(k-1)))
%    Geometricamente, traca-se uma reta secante pelos pontos (x_(k-1), f(x_(k-1)))
%    e (x_k, f(x_k)), determinando-se a intersecao dessa reta com o eixo x.
%    A ordem de convergencia assintotica e de aproximadamente 1.618 (razao aurea).
%
% 2. O QUE ENTREGA:
%    - Tabela formatada passo a passo exibida no console com:
%      Iteracao, x_(k-1), x_k, x_(k+1), f(x_(k+1)), Erro Absoluto e Erro Relativo.
%    - Raiz aproximada final ('raiz_aproximada_final').
%    - Total de iteracoes executadas.
%    - Mensagens informativas sobre o encerramento do processo iterativo.
%
% 3. CORRECOES APLICADAS EM RELACAO A VERSAO LEGADA:
%    [BUG CRITICO 1 CORRIGIDO]: A variavel 'xz' so era criada caso a condicao
%    de parada fosse atingida antes do fim do laco. Se o numero maximo de ciclos
%    esgotasse sem convergir, 'xz' nao existia e a impressao final falhava.
%    Nesta versao, 'raiz_aproximada_final' e inicializada previamente com NaN.
%    [BUG CRITICO 2 CORRIGIDO]: Na versao legada da raiz do projeto, a expressao
%    de erro era 'abs(xi-b)/b' sem modulo no denominador, gerando erros negativos
%    ou incorretos caso 'b' fosse negativo. Corrigido com modulo total e protecao
%    para valores proximos de zero.
%    [ROBUSTEZ ADICIONADA]: Protecao explicita contra divisao por zero caso
%    |f(x_k) - f(x_(k-1))| seja menor que a precisao de maquina.
%
% 4. PARAMETROS E PRE-REQUISITOS:
%    - Funcao f(x) continua na vizinhanca da raiz.
%    - Duas estimativas iniciais distintas: x0 e x1.
%    - f(x0) ~= f(x1) para evitar divisao por zero imediata.
%    - Tolerancia de convergencia e numero maximo de iteracoes admitidas.
%
% 5. EXEMPLO DE USO E RESULTADO ESPERADO:
%    Funcao: f(x) = 2*x^3 - 11.7*x^2 + 17.7*x - 5
%    Estimativa 0: x0 = 0.0
%    Estimativa 1: x1 = 1.0
%    Tolerancia: 1e-3
%    Resultado esperado: raiz aproximada x ~= 0.365099 (em 7 iteracoes).
% =========================================================================

% Limpeza inicial do ambiente de trabalho do Octave / MATLAB
clc;        % Limpa o console de texto
clear all;  % Remove todas as variaveis da memoria do workspace
close all;  % Fecha todas as janelas graficas abertas

%% =========================================================================
%% BLOCO DE PARAMETROS DE ENTRADA DO USUARIO
%% =========================================================================

% Funcao matematica f(x) cuja raiz real deseja-se calcular
funcao_objetivo = @(x)(2.0*x.^3 - 11.7*x.^2 + 17.7*x - 5.0);

% Primeira estimativa inicial (ponto x0)
ponto_inicial_a = 0.0;

% Segunda estimativa inicial (ponto x1)
ponto_inicial_b = 1.0;

% Tolerancia para os criterios de parada (erro relativo e residuo)
tolerancia_convergencia = 1e-3;

% Limite maximo de iteracoes permitidas para evitar lacos infinitos
numero_maximo_iteracoes = 50;

%% =========================================================================
%% BLOCO DE PROCESSAMENTO E CALCULO NUMERICO
%% =========================================================================

% Inicializacao da variavel de saida da raiz com NaN
raiz_aproximada_final = NaN;

% Variavel booleana de estado para registrar se a convergencia foi obtida
convergencia_atingida = false;

% Contador do total de iteracoes executadas
total_iteracoes_executadas = 0;

% Atribuicao dos pontos de trabalho a partir das entradas
estimativa_anterior_a = ponto_inicial_a;
estimativa_anterior_b = ponto_inicial_b;

% Inicializacao defensiva de variaveis para prevencao de falha em saida prematura
proxima_estimativa_x = ponto_inicial_b;
erro_absoluto_passo = Inf;
erro_relativo_passo = Inf;

% Impressao do cabecalho explicativo da tabela de iteracoes
fprintf('\n=========================================================================================================================\n');
fprintf('                           METODO DA SECANTE - CONVERGENCIA PASSO A PASSO                                                \n');
fprintf('=========================================================================================================================\n');
fprintf('  Iter |      x_(k-1)      |        x_k        |      x_(k+1)      |     f(x_(k+1))    |   Erro Absoluto   |   Erro Relativo\n');
fprintf('-------+-------------------+-------------------+-------------------+-------------------+-------------------+----------------\n');

% Laco de iteracoes do metodo da secante
for iteracao_atual = 1:numero_maximo_iteracoes
    total_iteracoes_executadas = iteracao_atual;

    % Avaliacao da funcao objetivo nos dois pontos mais recentes
    valor_funcao_ponto_a = funcao_objetivo(estimativa_anterior_a);
    valor_funcao_ponto_b = funcao_objetivo(estimativa_anterior_b);

    % Calculo da diferenca entre os valores da funcao no denominador
    diferenca_valores_funcao = valor_funcao_ponto_b - valor_funcao_ponto_a;

    % Verificacao de seguranca contra divisao por zero (reta secante horizontal)
    if abs(diferenca_valores_funcao) < eps
        fprintf('\n[ERRO NUMERICO]: Divisao por zero detectada. f(x_k) e f(x_(k-1)) sao numericamente iguais.\n');
        fprintf('A reta secante e perfeitamente horizontal. O algoritmo foi interrompido.\n');
        break;
    end

    % Formula iterativa do metodo da secante: x_(k+1) = x_k - f(x_k)*(x_k - x_(k-1)) / (f(x_k) - f(x_(k-1)))
    proxima_estimativa_x = estimativa_anterior_b - ...
        valor_funcao_ponto_b * (estimativa_anterior_b - estimativa_anterior_a) / diferenca_valores_funcao;

    % Avaliacao da funcao objetivo na nova raiz estimada
    valor_funcao_proxima_estimativa = funcao_objetivo(proxima_estimativa_x);

    % Calculo do erro absoluto entre aproximacoes consecutivas
    erro_absoluto_passo = abs(proxima_estimativa_x - estimativa_anterior_b);

    % Calculo robusto do erro relativo com protecao contra raiz nula ou negativa
    if abs(estimativa_anterior_b) > eps
        erro_relativo_passo = abs(erro_absoluto_passo / estimativa_anterior_b);
    else
        erro_relativo_passo = erro_absoluto_passo;
    end

    % Exibicao dos dados formatados da iteracao corrente
    fprintf('  %4d |   %15.8f |   %15.8f |   %15.8f |   %15.8e |   %17.8e |   %14.8e\n', ...
            iteracao_atual, ...
            estimativa_anterior_a, ...
            estimativa_anterior_b, ...
            proxima_estimativa_x, ...
            valor_funcao_proxima_estimativa, ...
            erro_absoluto_passo, ...
            erro_relativo_passo);

    % Criterio de parada: residuo proximo de zero OU erro relativo inferior a tol OU erro absoluto inferior a tol
    if (abs(valor_funcao_proxima_estimativa) < tolerancia_convergencia) || ...
       (erro_relativo_passo < tolerancia_convergencia) || ...
       (erro_absoluto_passo < tolerancia_convergencia)
        raiz_aproximada_final = proxima_estimativa_x;
        convergencia_atingida = true;
        break; % Convergencia atingida com sucesso
    end

    % Atualizacao dos pontos para o proximo ciclo: x_(k-1) recebe x_k e x_k recebe x_(k+1)
    estimativa_anterior_a = estimativa_anterior_b;
    estimativa_anterior_b = proxima_estimativa_x;
end

fprintf('-------------------------------------------------------------------------------------------------------------------------\n');

% Exibicao do relatorio final de resultados
if convergencia_atingida
    fprintf('\n>>> SUCESSO: O Metodo da Secante convergiu satisfatoriamente!\n');
    fprintf('    Raiz aproximada encontrada : %.10f\n', raiz_aproximada_final);
    fprintf('    Residuo f(raiz)            : %.10e\n', funcao_objetivo(raiz_aproximada_final));
    fprintf('    Ultimo erro absoluto       : %.10e (Tolerancia: %.1e)\n', erro_absoluto_passo, tolerancia_convergencia);
    fprintf('    Ultimo erro relativo       : %.10e\n', erro_relativo_passo);
    fprintf('    Total de iteracoes         : %d\n\n', total_iteracoes_executadas);
else
    fprintf('\n>>> AVISO: O metodo nao convergiu dentro do limite de %d iteracoes.\n', numero_maximo_iteracoes);
    fprintf('    Ultima estimativa calculada: %.10f\n', proxima_estimativa_x);
    fprintf('    Erro absoluto alcancado    : %.10e (Tolerancia: %.1e)\n', erro_absoluto_passo, tolerancia_convergencia);
    fprintf('    Erro relativo alcancado    : %.10e\n\n', erro_relativo_passo);
end
