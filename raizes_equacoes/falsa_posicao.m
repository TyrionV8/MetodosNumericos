% =========================================================================
% METODOS NUMERICOS - CALCULO DE RAIZES DE EQUACOES NAO LINEARES
% METODO DA FALSA POSICAO (REGULA FALSI)
% =========================================================================
% 
% 1. O QUE CALCULA:
%    Determina a raiz real da equacao nao linear f(x) = 0 contida em um
%    intervalo [limite_inferior, limite_superior] com troca de sinal
%    (Teorema de Bolzano: f(a) * f(b) < 0).
%    Diferente da bissecao (que divide o intervalo geometricamente ao meio),
%    o metodo da falsa posicao traca uma reta secante unindo os pontos
%    (a, f(a)) e (b, f(b)). A raiz aproximada e a intersecao dessa reta
%    com o eixo das abscissas (y = 0):
%        x = (a * f(b) - b * f(a)) / (f(b) - f(a))
%          = b - (f(b) * (b - a)) / (f(b) - f(a))
%
% 2. O QUE ENTREGA:
%    - Tabela formatada passo a passo exibida no console com:
%      Iteracao, Limite Inferior, Limite Superior, Raiz Estimada, f(Raiz), Erro Absoluto (dx) e Erro Relativo.
%    - Raiz aproximada final ('raiz_aproximada_final').
%    - Numero de iteracoes necessarias para a convergencia.
%    - Mensagens informativas de convergencia e verificacoes de consistencia.
%
% 3. CORRECOES APLICADAS EM RELACAO A VERSAO LEGADA:
%    [BUG CRITICO 1 CORRIGIDO]: Na versao anterior, f(a) e f(b) eram avaliados
%    apenas uma vez antes do laco for. A cada iteracao, quando os limites do
%    intervalo mudavam, a formula da reta continuava usando os f(a) e f(b) antigos.
%    Nesta versao, os valores da funcao nos limites sao recalculados e atualizados
%    a cada passo.
%    [BUG CRITICO 2 CORRIGIDO]: Na versao anterior, o erro era calculado como (b-a)/2,
%    o que e valido apenas na bissecao. Na falsa posicao, um dos extremos pode ficar
%    fixo (fenomeno do extremo estagnado em curvas convexas/concavas), fazendo com que
%    (b - a) nao convirja para zero. O erro correto e o modulo do deslocamento entre
%    duas aproximacoes sucessivas da raiz: erro = |x_atual - x_anterior|.
%
% 4. PARAMETROS E PRE-REQUISITOS:
%    - Funcao f(x) continua no intervalo fechado.
%    - Teorema de Bolzano satisfeito: f(limite_inferior) * f(limite_superior) < 0.
%    - Tolerancia de convergencia e numero maximo de iteracoes.
%
% 5. EXEMPLO DE USO E RESULTADO ESPERADO:
%    Funcao: f(x) = x^2 - exp(-x)
%    Intervalo: [0.0, 1.0]
%    Tolerancia: 1e-4
%    Resultado esperado: raiz aproximada x ~= 0.703467 (convergencia em cerca de 5-6 iteracoes).
% =========================================================================

% Limpeza inicial do ambiente de trabalho do Octave / MATLAB
clc;        % Limpa o console de texto
clear all;  % Limpa todas as variaveis da memoria
close all;  % Fecha todas as figuras graficas abertas

%% =========================================================================
%% BLOCO DE PARAMETROS DE ENTRADA DO USUARIO
%% =========================================================================

% Funcao matematica continua f(x) cuja raiz real deseja-se calcular
funcao_objetivo = @(x)(x.^2 - exp(-x));

% Limite inferior do intervalo inicial que contem a raiz
limite_inferior_intervalo = 0.0;

% Limite superior do intervalo inicial que contem a raiz
limite_superior_intervalo = 1.0;

% Tolerancia para o criterio de parada (baseada no deslocamento da raiz)
tolerancia_convergencia = 1e-4;

% Numero maximo de iteracoes para prevencao de loops infinitos
numero_maximo_iteracoes = 50;

%% =========================================================================
%% BLOCO DE PROCESSAMENTO E CALCULO NUMERICO
%% =========================================================================

% Inicializacao da raiz encontrada com NaN para evitar quebrar caso nao convirja
raiz_aproximada_final = NaN;

% Indicador booleano de sucesso na convergencia
convergencia_atingida = false;

% Contador de iteracoes efetivamente realizadas
total_iteracoes_executadas = 0;

% Variavel para guardar a estimativa da iteracao anterior (iniciada como NaN)
raiz_estimada_anterior = NaN;

% Inicializacao defensiva das variaveis de estimativa e erro
raiz_estimada_reta_secante = NaN;
erro_deslocamento_raiz = Inf;
erro_relativo_deslocamento = Inf;

% Avaliacao da funcao nos limites iniciais do intervalo
valor_funcao_limite_inferior = funcao_objetivo(limite_inferior_intervalo);
valor_funcao_limite_superior = funcao_objetivo(limite_superior_intervalo);

% Verificacao se algum dos extremos do intervalo ja e raiz exata (dentro da precisao de maquina)
if abs(valor_funcao_limite_inferior) <= eps
    raiz_aproximada_final = limite_inferior_intervalo;
    fprintf('A raiz exata foi encontrada diretamente no limite inferior: x = %.6f\n', raiz_aproximada_final);
    return;
elseif abs(valor_funcao_limite_superior) <= eps
    raiz_aproximada_final = limite_superior_intervalo;
    fprintf('A raiz exata foi encontrada diretamente no limite superior: x = %.6f\n', raiz_aproximada_final);
    return;
end

% Verificacao do Teorema de Bolzano para troca de sinal
if (valor_funcao_limite_inferior * valor_funcao_limite_superior) > 0.0
    fprintf('\n[ERRO CRITICO]: O Teorema de Bolzano nao e satisfeito no intervalo [%.6f, %.6f].\n', ...
            limite_inferior_intervalo, limite_superior_intervalo);
    fprintf('Motivo: f(a) = %.6f e f(b) = %.6f possuem o mesmo sinal algebraico.\n', ...
            valor_funcao_limite_inferior, valor_funcao_limite_superior);
    fprintf('Nao e possivel garantir a convergencia para a raiz neste intervalo.\n\n');
    return;
end

% Impressao do cabecalho da tabela informativa
fprintf('\n============================================================================================================================\n');
fprintf('                           METODO DA FALSA POSICAO - CONVERGENCIA PASSO A PASSO                                             \n');
fprintf('============================================================================================================================\n');
fprintf('  Iter |    Lim. Inferior  |    Lim. Superior  |    Raiz Estimada  |      f(Raiz)      |   Erro Abs. (dx)  |   Erro Relativo   \n');
fprintf('-------+-------------------+-------------------+-------------------+-------------------+-------------------+-------------------\n');

% Laco iterativo do metodo da falsa posicao
for iteracao_atual = 1:numero_maximo_iteracoes
    total_iteracoes_executadas = iteracao_atual;

    % CORRECAO DE BUG 1: Atualiza/reavalia f(a) e f(b) com os limites atuais do intervalo
    valor_funcao_limite_inferior = funcao_objetivo(limite_inferior_intervalo);
    valor_funcao_limite_superior = funcao_objetivo(limite_superior_intervalo);

    % Prevencao de divisao por zero caso f(b) e f(a) se tornem numericamente identicos
    diferenca_valores_funcao = valor_funcao_limite_superior - valor_funcao_limite_inferior;
    if abs(diferenca_valores_funcao) < eps
        fprintf('\n[AVISO]: Diferenca entre f(b) e f(a) e numericamente nula. Interrupcao para evitar divisao por zero.\n');
        break;
    end

    % Formula da interpolacao linear (zero da reta secante que liga (a, f(a)) a (b, f(b)))
    raiz_estimada_reta_secante = (limite_inferior_intervalo * valor_funcao_limite_superior - ...
                                  limite_superior_intervalo * valor_funcao_limite_inferior) / ...
                                  diferenca_valores_funcao;

    % Avaliacao da funcao objetivo na nova raiz estimada
    valor_funcao_raiz_estimada = funcao_objetivo(raiz_estimada_reta_secante);

    % CORRECAO DE BUG 2: Calculo do erro baseado no deslocamento da raiz (|x_novo - x_anterior|)
    if iteracao_atual == 1
        % Na primeira iteracao, toma-se a largura total do intervalo inicial como referencia
        erro_deslocamento_raiz = abs(limite_superior_intervalo - limite_inferior_intervalo);
        if abs(raiz_estimada_reta_secante) > eps
            erro_relativo_deslocamento = abs(erro_deslocamento_raiz / raiz_estimada_reta_secante);
        else
            erro_relativo_deslocamento = erro_deslocamento_raiz;
        end
    else
        % A partir da segunda iteracao, o erro real e o deslocamento absoluto entre estimativas
        erro_deslocamento_raiz = abs(raiz_estimada_reta_secante - raiz_estimada_anterior);
        if abs(raiz_estimada_reta_secante) > eps
            erro_relativo_deslocamento = abs(erro_deslocamento_raiz / raiz_estimada_reta_secante);
        else
            erro_relativo_deslocamento = erro_deslocamento_raiz;
        end
    end

    % Exibicao dos dados da iteracao atual formatados em colunas
    fprintf('  %4d |   %15.8f |   %15.8f |   %15.8f |   %15.8e |   %17.8e |   %17.8e\n', ...
            iteracao_atual, ...
            limite_inferior_intervalo, ...
            limite_superior_intervalo, ...
            raiz_estimada_reta_secante, ...
            valor_funcao_raiz_estimada, ...
            erro_deslocamento_raiz, ...
            erro_relativo_deslocamento);

    % Teste de parada: verifica se o deslocamento absoluto, relativo ou valor da funcao atingiram a tolerancia
    if (iteracao_atual > 1 && (erro_deslocamento_raiz <= tolerancia_convergencia || erro_relativo_deslocamento <= tolerancia_convergencia)) || ...
       (abs(valor_funcao_raiz_estimada) <= tolerancia_convergencia)
        raiz_aproximada_final = raiz_estimada_reta_secante;
        convergencia_atingida = true;
        break; % Convergencia alcancada com sucesso
    end

    % Armazena a estimativa atual para calcular o deslocamento na proxima iteracao
    raiz_estimada_anterior = raiz_estimada_reta_secante;

    % Atualizacao dos limites do intervalo mantendo a troca de sinal (Bolzano)
    if (valor_funcao_limite_inferior * valor_funcao_raiz_estimada) < 0.0
        % A raiz esta entre o limite inferior e a estimativa atual da reta
        limite_superior_intervalo = raiz_estimada_reta_secante;
    else
        % A raiz esta entre a estimativa atual da reta e o limite superior
        limite_inferior_intervalo = raiz_estimada_reta_secante;
    end
end

fprintf('----------------------------------------------------------------------------------------------------------------------------\n');

% Exibicao dos resultados finais
if convergencia_atingida
    fprintf('\n>>> SUCESSO: O Metodo da Falsa Posicao convergiu com exito!\n');
    fprintf('    Raiz aproximada encontrada : %.10f\n', raiz_aproximada_final);
    fprintf('    Residuo f(raiz)            : %.10e\n', funcao_objetivo(raiz_aproximada_final));
    fprintf('    Ultimo erro de deslocamento: %.10e (Tolerancia: %.1e)\n', erro_deslocamento_raiz, tolerancia_convergencia);
    fprintf('    Ultimo erro relativo       : %.10e\n', erro_relativo_deslocamento);
    fprintf('    Total de iteracoes         : %d\n\n', total_iteracoes_executadas);
else
    fprintf('\n>>> AVISO: Numero maximo de iteracoes (%d) atingido sem convergencia satisfatoria.\n', numero_maximo_iteracoes);
    if ~isnan(raiz_estimada_reta_secante)
        fprintf('    Ultima raiz calculada      : %.10f\n', raiz_estimada_reta_secante);
        fprintf('    Erro remanescente          : %.10e (Tolerancia: %.1e)\n', erro_deslocamento_raiz, tolerancia_convergencia);
        fprintf('    Erro relativo remanescente : %.10e\n\n', erro_relativo_deslocamento);
    else
        fprintf('    Nenhuma aproximacao valida foi calculada antes da interrupcao.\n\n');
    end
end
