function estrutura_resultado_conversoes = conversao_si_imperial(massa_quilogramas_entrada, comprimento_metros_entrada)
% ==============================================================================
% UTILITÁRIOS: CONVERSÃO DE UNIDADES DO SISTEMA INTERNACIONAL (SI) PARA IMPERIAL
% ==============================================================================
%
% 1. O QUE CALCULA:
%    Esta rotina converte grandezas físicas fundamentais de massa e comprimento
%    expressas nas unidades padrão do Sistema Internacional (quilogramas e metros)
%    para as unidades do Sistema Imperial Britânico / Tradicional Americano:
%    - Massa: Quilogramas [kg] -> Libras-massa [lb] e Onças [oz].
%    - Comprimento: Metros [m] -> Polegadas [in], Pés [ft] e Jardas [yd].
%
%    Formulação e Fatores de Conversão Padronizados:
%    - Pelo acordo internacional de 1959 (International Yard and Pound Agreement):
%      1 libra (lb) = 0.45359237 kg  ==>  1 kg ≈ 2.20462262 lb (aprox. 2.205 lb)
%      1 onça (oz)  = 1/16 lb         ==>  1 kg ≈ 35.27396 oz
%      1 polegada (in) = 0.0254 m     ==>  1 m  ≈ 39.3700787 in (aprox. 39.37 in)
%      1 pé (ft) = 12 polegadas = 0.3048 m ==> 1 m ≈ 3.28084 ft
%      1 jarda (yd) = 3 pés = 0.9144 m     ==> 1 m ≈ 1.09361 yd
%
%    Esta implementação corrige o bug sintático crítico da versão legada da raiz,
%    onde vírgulas foram digitadas no lugar de pontos decimais ('2,205' e '39,37').
%
% 2. O QUE ENTREGA:
%    - Se executado sem argumentos de entrada (nargin == 0):
%      Executa demonstração didática com impressões detalhadas no console.
%    - Se chamado como função com 2 argumentos:
%      Retorna uma estrutura organizada (estrutura_resultado_conversoes) e vetor
%      [massa_libras, comprimento_polegadas] com todas as unidades derivadas.
%
% 3. O QUE PRECISA TER (PRÉ-REQUISITOS E PARÂMETROS):
%    - massa_quilogramas_entrada: Escalar positivo em kg (ex.: 75.0 kg).
%    - comprimento_metros_entrada: Escalar positivo em metros (ex.: 1.80 m).
%
% 4. EXEMPLO DE USO E RESULTADO ESPERADO:
%    - Exemplo de chamada como script:
%      >> conversao_si_imperial
%    - Exemplo de chamada como função:
%      >> res = conversao_si_imperial(80.0, 1.82);
%      res.massa_libras ≈ 176.37 lb
%      res.comprimento_polegadas ≈ 71.65 in (5 ft 11.65 in)
% ==============================================================================

% ==============================================================================
% DEFINIÇÃO DOS FATORES DE CONVERSÃO EXATOS E CONVENCIONAIS
% ==============================================================================

% Fator exato de conversão de massa: quilogramas para libras
fator_quilogramas_para_libras = 1.0 / 0.45359237; % Aprox. 2.20462262 lb/kg

% Fator de conversão de massa: quilogramas para onças (1 lb = 16 oz)
fator_quilogramas_para_oncas = fator_quilogramas_para_libras * 16.0;

% Fator exato de conversão de comprimento: metros para polegadas
fator_metros_para_polegadas = 1.0 / 0.0254; % Aprox. 39.3700787 in/m

% Fator de conversão de comprimento: metros para pés (1 ft = 12 in = 0.3048 m)
fator_metros_para_pes = 1.0 / 0.3048; % Aprox. 3.28084 ft/m

% Fator de conversão de comprimento: metros para jardas (1 yd = 3 ft = 0.9144 m)
fator_metros_para_jardas = 1.0 / 0.9144; % Aprox. 1.09361 yd/m

% ==============================================================================
% TRATAMENTO DE PARÂMETROS DE ENTRADA (SUPORTE A DEMONSTRAÇÃO AUTÔNOMA)
% ==============================================================================

modo_demonstracao_ativo = false;

if nargin == 0
    % Modo demonstração ativado: carrega valores antropométricos típicos de teste
    modo_demonstracao_ativo = true;
    massa_quilogramas_local = 70.0;  % Massa corporal de referência: 70 kg
    comprimento_metros_local = 1.75; % Altura de referência: 1.75 m
elseif nargin == 1
    massa_quilogramas_local = massa_quilogramas_entrada;
    comprimento_metros_local = 1.0; % Comprimento padrão unitário se omitido
else
    massa_quilogramas_local = massa_quilogramas_entrada;
    comprimento_metros_local = comprimento_metros_entrada;
end

% Validação de não negatividade das grandezas físicas
if massa_quilogramas_local < 0.0
    error('A massa fornecida não pode ser negativa.');
end
if comprimento_metros_local < 0.0
    error('O comprimento fornecido não pode ser negativo.');
end

% ==============================================================================
% PROCESSAMENTO DAS CONVERSÕES DE UNIDADES
% ==============================================================================

% Conversão de grandezas de massa
massa_libras_calculada   = massa_quilogramas_local * fator_quilogramas_para_libras;
massa_oncas_calculada    = massa_quilogramas_local * fator_quilogramas_para_oncas;

% Conversão de grandezas de comprimento
comprimento_polegadas_calculado = comprimento_metros_local * fator_metros_para_polegadas;
comprimento_pes_calculado       = comprimento_metros_local * fator_metros_para_pes;
comprimento_jardas_calculado    = comprimento_metros_local * fator_metros_para_jardas;

% Decomposição didática de altura em pés inteiros e polegadas remanescentes
pes_inteiros_altura          = floor(comprimento_pes_calculado);
polegadas_excedentes_altura = (comprimento_pes_calculado - pes_inteiros_altura) * 12.0;

% ==============================================================================
% EMPACOTAMENTO DA ESTRUTURA DE RETORNO
% ==============================================================================

estrutura_resultado_conversoes = struct();
estrutura_resultado_conversoes.massa_original_quilogramas       = massa_quilogramas_local;
estrutura_resultado_conversoes.massa_convertida_libras          = massa_libras_calculada;
estrutura_resultado_conversoes.massa_convertida_oncas           = massa_oncas_calculada;
estrutura_resultado_conversoes.comprimento_original_metros      = comprimento_metros_local;
estrutura_resultado_conversoes.comprimento_convertido_polegadas = comprimento_polegadas_calculado;
estrutura_resultado_conversoes.comprimento_convertido_pes       = comprimento_pes_calculado;
estrutura_resultado_conversoes.comprimento_convertido_jardas    = comprimento_jardas_calculado;

% Vetor compacto legado para retrocompatibilidade direta com [massa_lb, altura_in]
estrutura_resultado_conversoes.vetor_compacto_legado            = [massa_libras_calculada, comprimento_polegadas_calculado];

% ==============================================================================
% EXIBIÇÃO DE RELATÓRIO QUANDO EXECUTADO EM MODO DEMONSTRAÇÃO OU CONSOLE
% ==============================================================================

if modo_demonstracao_ativo || nargout == 0
    fprintf('\n==============================================================================\n');
    fprintf('       UTILITÁRIOS: CONVERSÃO DE UNIDADES (SISTEMA INTERNACIONAL -> IMPERIAL) \n');
    fprintf('==============================================================================\n');
    fprintf('  GRANDEZAS DE ENTRADA NO SI:\n');
    fprintf('    - Massa fornecida:                 %10.4f kg\n', massa_quilogramas_local);
    fprintf('    - Comprimento/Altura fornecido:    %10.4f m\n', comprimento_metros_local);
    fprintf('------------------------------------------------------------------------------\n');
    fprintf('  CONVERSÕES EQUIVALENTES NO SISTEMA IMPERIAL:\n');
    fprintf('    - Massa em Libras (lb):            %10.4f lb\n', massa_libras_calculada);
    fprintf('    - Massa em Onças (oz):             %10.4f oz\n', massa_oncas_calculada);
    fprintf('    - Comprimento em Polegadas (in):   %10.4f in\n', comprimento_polegadas_calculado);
    fprintf('    - Comprimento em Pés (ft):         %10.4f ft\n', comprimento_pes_calculado);
    fprintf('    - Altura mista (Pés e Polegadas):   %d ft %4.1f in\n', pes_inteiros_altura, polegadas_excedentes_altura);
    fprintf('    - Comprimento em Jardas (yd):      %10.4f yd\n', comprimento_jardas_calculado);
    fprintf('------------------------------------------------------------------------------\n');
    fprintf('  VETOR COMPACTO DE SAÍDA (LEGADO):   [%10.4f lb, %10.4f in]\n', ...
            massa_libras_calculada, comprimento_polegadas_calculado);
    fprintf('==============================================================================\n\n');
end

end
