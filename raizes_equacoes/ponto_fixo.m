% =========================================================================
% METODOS NUMERICOS - CALCULO DE RAIZES DE EQUACOES NAO LINEARES
% METODO DO PONTO FIXO (ITERACAO LINEAR / PICARD-BANACH)
% =========================================================================
% 
% 1. O QUE CALCULA:
%    Determina uma raiz real da equacao f(x) = 0 reescrevendo o problema
%    algebricamente na forma equivalente de ponto fixo x = g(x).
%    O processo iterativo e regido pela relacao de recorrencia:
%        x_(k+1) = g(x_k)
%    Pelo Teorema do Ponto Fixo de Banach, a sequencia de aproximacoes converge
%    para o ponto fixo unico x* (onde x* = g(x*) e f(x*) = 0) se a funcao
%    de iteracao g(x) satisfizer a condicao de contracao de Lipschitz:
%        |g'(x)| <= L < 1
%    em uma vizinhanca fechada da raiz contendo a estimativa inicial x0.
%
% 2. O QUE ENTREGA:
%    - Tabela formatada passo a passo exibida no console com as colunas:
%      Iteracao, x_k, x_(k+1) = g(x_k), |f(x_(k+1))|, |x_(k+1) - x_k| e Erro Relativo.
%    - Raiz aproximada final ('raiz_aproximada_final').
%    - Numero total de iteracoes executadas.
%    - Diagnostico conclusivo sobre a convergencia ou deteccao de divergencia.
%
% 3. CORRECOES APLICADAS EM RELACAO A VERSAO LEGADA:
%    [BUG CRITICO 1 CORRIGIDO]: A solucao 'xz' era definida apenas quando o if
%    de convergencia era acionado. Caso o limite maximo de iteracoes fosse
%    atingido sem convergir, 'xz' nao existia e o fprintf final falhava.
%    Aqui, 'raiz_aproximada_final' e inicializada previamente com NaN.
%    [BUG CRITICO 2 CORRIGIDO]: Na versao legada havia uma linha solta
%    'abs(F(G(xi)))' sem ponto e virgula que poluia a saida com impressao
%    desformatada a cada ciclo. Essa chamada foi removida e incorporada
%    a tabela formatada.
%    [CRITERIOS DE PARADA ROBUSTOS]: Na versao legada verificava-se apenas
%    |f(g(x))| < tol. Foram adicionados:
%      a) Verificacao do deslocamento absoluto |x_(k+1) - x_k| <= tolerancia.
%      b) Verificacao do residuo da funcao objetivo |f(x_(k+1))| <= tolerancia.
%      c) Deteccao preventiva de divergencia (NaN, Inf ou valores > 1e12).
%
% 4. PARAMETROS E PRE-REQUISITOS:
%    - Funcao objetivo f(x) e funcao de iteracao equivalente g(x) tal que f(x)=0 <=> x=g(x).
%    - Estimativa inicial x0 no intervalo onde |g'(x)| < 1.
%    - Tolerancia de convergencia e numero maximo de iteracoes.
%
% 5. EXEMPLO DE USO E RESULTADO ESPERADO:
%    Funcao objetivo : f(x) = exp(x) - 4*x^2
%    Funcao iteracao : g(x) = 0.5 * exp(0.5 * x)
%    Estimativa inicial : x0 = 0.5
%    Tolerancia : 1e-3
%    Resultado esperado: raiz aproximada x ~= 0.714393 (em 6 iteracoes).
% =========================================================================

% Limpeza inicial do ambiente de trabalho do Octave / MATLAB
clc;        % Limpa o console de texto
clear all;  % Remove todas as variaveis do workspace
close all;  % Fecha todas as figuras graficas abertas

%% =========================================================================
%% BLOCO DE PARAMETROS DE ENTRADA DO USUARIO
%% =========================================================================

% Funcao matematica objetivo original f(x) = 0
funcao_objetivo = @(x)(exp(x) - 4.0 * x.^2);

% Funcao de iteracao de ponto fixo equivalente g(x) tal que x = g(x)
funcao_iteracao_ponto_fixo = @(x)(0.5 * exp(0.5 * x));

% Estimativa inicial para o inicio da sequencia iterativa
estimativa_inicial_x = 0.5;

% Tolerancia para os criterios de convergencia
tolerancia_convergencia = 1e-3;

% Limite maximo de iteracoes admitidas
numero_maximo_iteracoes = 50;

%% =========================================================================
%% BLOCO DE PROCESSAMENTO E CALCULO NUMERICO
%% =========================================================================

% Inicializacao da raiz aproximada com NaN para consistencia caso nao convirja
raiz_aproximada_final = NaN;

% Variavel booleana de estado para registrar se a convergencia foi alcancada
convergencia_atingida = false;

% Contador do total de iteracoes executadas
total_iteracoes_executadas = 0;

% Inicializacao do ponto de busca com a estimativa inicial
estimativa_atual_x = estimativa_inicial_x;

% Inicializacao defensiva das variaveis para protecao em caso de divergencia precoce
proxima_estimativa_x = estimativa_inicial_x;
erro_deslocamento_passo = Inf;
erro_relativo_passo = Inf;

% Impressao do cabecalho explicativo da tabela de iteracoes
fprintf('\n===================================================================================================\n');
fprintf('                 METODO DO PONTO FIXO - CONVERGENCIA PASSO A PASSO                                 \n');
fprintf('===================================================================================================\n');
fprintf('  Iter |       x_k         |    x_(k+1)=g(x_k) |   |f(x_(k+1))|    |  |x_(k+1) - x_k|  | Erro Relativo \n');
fprintf('-------+-------------------+-------------------+-------------------+-------------------+---------------\n');

% Laco principal de repeticao do metodo do ponto fixo
for iteracao_atual = 1:numero_maximo_iteracoes
    total_iteracoes_executadas = iteracao_atual;

    % Calculo da proxima estimativa pela funcao de iteracao: x_(k+1) = g(x_k)
    proxima_estimativa_x = funcao_iteracao_ponto_fixo(estimativa_atual_x);

    % Deteccao preventiva de divergencia numerica (overflow ou indeterminacao)
    if isnan(proxima_estimativa_x) || isinf(proxima_estimativa_x) || abs(proxima_estimativa_x) > 1e12
        fprintf('\n[ALERTA DE DIVERGENCIA]: A sequencia iterativa divergiu na iteracao %d (|x| > 1e12 ou NaN/Inf).\n', iteracao_atual);
        fprintf('Possivel causa: A condicao de contracao |g''(x)| < 1 nao e satisfeita nesta regiao.\n');
        break;
    end

    % Calculo do deslocamento absoluto entre passos sucessivos
    erro_deslocamento_passo = abs(proxima_estimativa_x - estimativa_atual_x);

    % Calculo do erro relativo com protecao contra divisao por zero
    if abs(proxima_estimativa_x) > eps
        erro_relativo_passo = abs(erro_deslocamento_passo / proxima_estimativa_x);
    else
        erro_relativo_passo = erro_deslocamento_passo;
    end

    % Avaliacao do residuo da funcao objetivo f(x) na nova estimativa
    valor_funcao_objetivo_atual = funcao_objetivo(proxima_estimativa_x);

    % Exibicao dos dados da iteracao formatados em colunas alinhadas
    fprintf('  %4d |   %15.8f |   %15.8f |   %15.8e |   %15.8e |   %11.8e\n', ...
            iteracao_atual, ...
            estimativa_atual_x, ...
            proxima_estimativa_x, ...
            abs(valor_funcao_objetivo_atual), ...
            erro_deslocamento_passo, ...
            erro_relativo_passo);

    % Criterios de parada robustos: satisfaz o residuo de f, a variacao absoluta de x ou erro relativo
    if (abs(valor_funcao_objetivo_atual) < tolerancia_convergencia) || ...
       (erro_deslocamento_passo < tolerancia_convergencia) || ...
       (erro_relativo_passo < tolerancia_convergencia)
        raiz_aproximada_final = proxima_estimativa_x;
        convergencia_atingida = true;
        break; % Convergencia alcancada com sucesso
    end

    % Atualiza a estimativa para o proximo ciclo iterativo
    estimativa_atual_x = proxima_estimativa_x;
end

fprintf('---------------------------------------------------------------------------------------------------\n');

% Exibicao dos resultados finais
if convergencia_atingida
    fprintf('\n>>> SUCESSO: O Metodo do Ponto Fixo convergiu com exito!\n');
    fprintf('    Raiz aproximada encontrada : %.10f\n', raiz_aproximada_final);
    fprintf('    Residuo |f(raiz)|          : %.10e\n', abs(funcao_objetivo(raiz_aproximada_final)));
    fprintf('    Ultimo deslocamento |dx|   : %.10e (Tolerancia: %.1e)\n', erro_deslocamento_passo, tolerancia_convergencia);
    fprintf('    Ultimo erro relativo       : %.10e\n', erro_relativo_passo);
    fprintf('    Total de iteracoes         : %d\n\n', total_iteracoes_executadas);
else
    fprintf('\n>>> AVISO: O metodo nao convergiu dentro do limite de %d iteracoes.\n', numero_maximo_iteracoes);
    if ~isnan(proxima_estimativa_x) && ~isinf(proxima_estimativa_x)
        fprintf('    Ultima estimativa obtida   : %.10f\n', proxima_estimativa_x);
        fprintf('    Ultimo deslocamento |dx|   : %.10e (Tolerancia: %.1e)\n', erro_deslocamento_passo, tolerancia_convergencia);
        fprintf('    Ultimo erro relativo       : %.10e\n\n', erro_relativo_passo);
    else
        fprintf('    Divergencia numerica confirmada: valor de x ultrapassou os limites finitos de calculo.\n\n');
    end
end
