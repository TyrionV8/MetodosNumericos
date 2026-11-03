% =========================================================================
% METODOS NUMERICOS - CALCULO DE RAIZES DE EQUACOES NAO LINEARES
% METODO DA BISSECAO (DICOTOMIA INTERVALAR)
% =========================================================================
% 
% 1. O QUE CALCULA:
%    Determina uma aproximacao numerica para a raiz real da equacao nao linear
%    f(x) = 0 contida em um intervalo fechado [limite_inferior, limite_superior].
%    O metodo fundamenta-se no Teorema de Bolzano (Teorema do Valor Intermediario):
%    se uma funcao f(x) e continua em [a, b] e f(a) * f(b) < 0, entao existe
%    pelo menos uma raiz x* pertencente ao intervalo aberto (a, b).
%    A cada passo, o intervalo e dividido ao meio:
%        ponto_medio = (limite_inferior + limite_superior) / 2
%    O erro absoluto maximo cometido na k-esima iteracao e delimitado por:
%        erro = (limite_superior - limite_inferior) / 2 = (b0 - a0) / (2^k)
%
% 2. O QUE ENTREGA:
%    - Tabela passo a passo exibida no console com as colunas:
%      Iteracao, Limite Inferior, Limite Superior, Ponto Medio, f(Ponto Medio), Erro Absoluto e Erro Relativo.
%    - Valor final da raiz aproximada ('raiz_aproximada_final').
%    - Numero total de iteracoes executadas para atingir a tolerancia.
%    - Mensagens claras sobre convergencia ou falha por esgotamento de iteracoes.
%
% 3. O QUE PRECISA TER (PARAMETROS E PRE-REQUISITOS):
%    - Funcao f(x) continua no intervalo de interesse.
%    - Troca de sinal nos extremos do intervalo: f(a) * f(b) < 0.
%    - Tolerancia desejada para o criterio de parada (tolerancia_convergencia).
%    - Limite maximo de iteracoes admitidas (numero_maximo_iteracoes).
%
% 4. EXEMPLO DE USO E RESULTADO ESPERADO:
%    Funcao: f(x) = x^2 - exp(-x)
%    Intervalo: [0.0, 1.0]
%    Tolerancia: 1e-4
%    Iteracoes maximas: 100
%    Resultado esperado: raiz aproximada x ~= 0.703491 (em 13-14 iteracoes).
% =========================================================================

% Limpeza inicial do ambiente de trabalho do Octave / MATLAB
clc;        % Limpa a tela do console de comandos
clear all;  % Remove todas as variaveis da memoria do workspace
close all;  % Fecha todas as janelas graficas que estejam abertas

%% =========================================================================
%% BLOCO DE PARAMETROS DE ENTRADA DO USUARIO
%% =========================================================================

% Definicao da funcao matematica continua cuja raiz f(x) = 0 deseja-se encontrar
funcao_objetivo = @(x)(x.^2 - exp(-x));

% Limite inferior do intervalo inicial de busca
limite_inferior_intervalo = 0.0;

% Limite superior do intervalo inicial de busca
limite_superior_intervalo = 1.0;

% Tolerancia maxima admitida para o erro absoluto da raiz aproximada
tolerancia_convergencia = 1e-4;

% Numero maximo de iteracoes para evitar laco infinito caso haja divergencia
numero_maximo_iteracoes = 100;

%% =========================================================================
%% BLOCO DE PROCESSAMENTO E CALCULO NUMERICO
%% =========================================================================

% Inicializacao da variavel que armazenara o resultado final da raiz
raiz_aproximada_final = NaN;

% Variavel booleana que indica se o metodo atingiu o criterio de convergencia
convergencia_atingida = false;

% Contador para armazenar quantas iteracoes foram efetivamente executadas
total_iteracoes_executadas = 0;

% Inicializacao defensiva das estimativas e erros antes do laco iterativo
ponto_medio_estimado = (limite_inferior_intervalo + limite_superior_intervalo) / 2.0;
erro_estimado_absoluto = (limite_superior_intervalo - limite_inferior_intervalo) / 2.0;
erro_estimado_relativo = Inf;

% Avaliacao da funcao objetivo nos dois extremos do intervalo inicial
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

% Verificacao da condicao do Teorema de Bolzano (troca de sinal f(a) * f(b) < 0)
if (valor_funcao_limite_inferior * valor_funcao_limite_superior) > 0.0
    fprintf('\n[ERRO CRITICO]: O Teorema de Bolzano nao e satisfeito no intervalo [%.6f, %.6f].\n', ...
            limite_inferior_intervalo, limite_superior_intervalo);
    fprintf('Motivo: f(a) = %.6f e f(b) = %.6f possuem o mesmo sinal algebraico.\n', ...
            valor_funcao_limite_inferior, valor_funcao_limite_superior);
    fprintf('Nao e possivel garantir a existencia ou convergencia para uma raiz neste intervalo.\n\n');
    return;
end

% Impressao do cabecalho explicativo da tabela de convergencia iterativa
fprintf('\n=====================================================================================================================\n');
fprintf('                             METODO DA BISSECAO - CONVERGENCIA PASSO A PASSO                                         \n');
fprintf('=====================================================================================================================\n');
fprintf('  Iter |    Lim. Inferior  |    Lim. Superior  |     Ponto Medio   |    f(Ponto Medio) |     Erro Abs.   |   Erro Rel.   \n');
fprintf('-------+-------------------+-------------------+-------------------+-------------------+-----------------+---------------\n');

% Laco principal de iteracoes do metodo da bissecao
for iteracao_atual = 1:numero_maximo_iteracoes
    % Registra o numero da iteracao corrente
    total_iteracoes_executadas = iteracao_atual;

    % Calculo do ponto medio do intervalo corrente (estimativa da raiz)
    ponto_medio_estimado = (limite_inferior_intervalo + limite_superior_intervalo) / 2.0;

    % Calculo do limitante superior do erro absoluto nesta iteracao: erro = (b - a) / 2
    erro_estimado_absoluto = (limite_superior_intervalo - limite_inferior_intervalo) / 2.0;

    % Calculo do erro relativo aproximado com protecao contra ponto medio nulo
    if abs(ponto_medio_estimado) > eps
        erro_estimado_relativo = erro_estimado_absoluto / abs(ponto_medio_estimado);
    else
        erro_estimado_relativo = erro_estimado_absoluto;
    end

    % Avaliacao da funcao objetivo no ponto medio calculado
    valor_funcao_ponto_medio = funcao_objetivo(ponto_medio_estimado);

    % Exibicao dos dados da iteracao formatados em colunas alinhadas
    fprintf('  %4d |   %15.8f |   %15.8f |   %15.8f |   %15.8e |   %13.8e |   %11.8e\n', ...
            iteracao_atual, ...
            limite_inferior_intervalo, ...
            limite_superior_intervalo, ...
            ponto_medio_estimado, ...
            valor_funcao_ponto_medio, ...
            erro_estimado_absoluto, ...
            erro_estimado_relativo);

    % Criterio de parada: erro absoluto <= tol OU erro relativo <= tol OU residuo da funcao nulo/desprezivel
    if (erro_estimado_absoluto <= tolerancia_convergencia) || ...
       (erro_estimado_relativo <= tolerancia_convergencia) || ...
       (abs(valor_funcao_ponto_medio) <= tolerancia_convergencia)
        raiz_aproximada_final = ponto_medio_estimado;
        convergencia_atingida = true;
        break; % Interrompe o laco pois a precisao requerida foi alcancada
    end

    % Aplicacao do Teorema de Bolzano para selecionar o subintervalo que preserva a troca de sinal
    if (valor_funcao_limite_inferior * valor_funcao_ponto_medio) < 0.0
        % A troca de sinal ocorre na metade esquerda [limite_inferior, ponto_medio]
        limite_superior_intervalo = ponto_medio_estimado;
        valor_funcao_limite_superior = valor_funcao_ponto_medio;
    else
        % A troca de sinal ocorre na metade direita [ponto_medio, limite_superior]
        limite_inferior_intervalo = ponto_medio_estimado;
        valor_funcao_limite_inferior = valor_funcao_ponto_medio;
    end
end

fprintf('---------------------------------------------------------------------------------------------------------------------\n');

% Exibicao do resultado final conforme o status de convergencia
if convergencia_atingida
    fprintf('\n>>> SUCESSO: Convergencia alcancada com exito!\n');
    fprintf('    Raiz aproximada encontrada : %.10f\n', raiz_aproximada_final);
    fprintf('    Residuo f(raiz)            : %.10e\n', funcao_objetivo(raiz_aproximada_final));
    fprintf('    Erro absoluto final        : %.10e (Tolerancia: %.1e)\n', erro_estimado_absoluto, tolerancia_convergencia);
    fprintf('    Erro relativo final        : %.10e\n', erro_estimado_relativo);
    fprintf('    Total de iteracoes         : %d\n\n', total_iteracoes_executadas);
else
    fprintf('\n>>> AVISO: O numero maximo de iteracoes (%d) foi atingido sem satisfazer a tolerancia!\n', numero_maximo_iteracoes);
    fprintf('    Ultima estimativa obtida   : %.10f\n', ponto_medio_estimado);
    fprintf('    Erro absoluto remanescente : %.10e (Tolerancia: %.1e)\n', erro_estimado_absoluto, tolerancia_convergencia);
    fprintf('    Erro relativo remanescente : %.10e\n\n', erro_estimado_relativo);
end
