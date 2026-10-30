# 🧮 Métodos Numéricos — Octave / MATLAB

Scripts em **GNU Octave** para resolver problemas de cálculo numérico: raízes de equações, sistemas lineares, integrais, derivadas, EDOs, além de alguns extras (sinais, química, farmacocinética).

> **Este README é um índice de consulta.** Use `Ctrl+F` com o nome do método, um sinônimo ou o que você quer fazer (ex.: *"zero da função"*, *"integral"*, *"pivô"*, *"RMS"*). Cada seção tem a fórmula, o arquivo, o que editar e o resultado esperado.

---

## 📑 Sumário

1. [⚡ Quero fazer X → abra o arquivo Y](#-quero-fazer-x--abra-o-arquivo-y)
2. [🗺️ Mapa do repositório](#️-mapa-do-repositório)
3. [▶️ Como rodar](#️-como-rodar)
4. [🎯 Raízes de equações — f(x) = 0](#-raízes-de-equações--fx--0)
5. [🧱 Sistemas lineares — Ax = b](#-sistemas-lineares--ax--b)
6. [∫ Integração numérica](#-integração-numérica)
7. [📐 Diferenciação numérica](#-diferenciação-numérica)
8. [📈 EDOs — equações diferenciais](#-edos--equações-diferenciais)
9. [🔁 Séries de Taylor](#-séries-de-taylor)
10. [〰️ Sinais: senoide, quadrada, triangular, RMS](#️-sinais-senoide-quadrada-triangular-rms)
11. [⚗️ Química / Física](#️-química--física)
12. [💊 Farmacocinética (pasta `pessoal/`)](#-farmacocinética-pasta-pessoal)
13. [📝 Exercícios de lista / aula e utilitários](#-exercícios-de-lista--aula-e-utilitários)
14. [📋 Cola rápida de Octave](#-cola-rápida-de-octave)
15. [🔄 Octave × MATLAB](#-octave--matlab)
16. [⚠️ Bugs conhecidos e pegadinhas](#️-bugs-conhecidos-e-pegadinhas)

---

## ⚡ Quero fazer X → abra o arquivo Y

| Quero... | Método | Arquivo |
|---|---|---|
| Achar **raiz / zero** de f(x) tendo um **intervalo [a, b]** com troca de sinal | Bisseção | [`Raizes Reais/Bissecao.m`](Raizes%20Reais/Bissecao.m) |
| Achar raiz com intervalo, convergindo mais rápido que a bisseção | Falsa Posição (Regula Falsi) | [`Raizes Reais/PosFalsa.m`](Raizes%20Reais/PosFalsa.m) |
| Achar raiz **sabendo a derivada** f'(x) | Newton-Raphson | [`Raizes Reais/Newton.m`](Raizes%20Reais/Newton.m) |
| Achar raiz **sem derivada** e sem precisar de troca de sinal | Secante | [`Raizes Reais/Secante.m`](Raizes%20Reais/Secante.m) |
| Achar raiz reescrevendo f(x)=0 como **x = g(x)** | Ponto Fixo (MPF) | [`Raizes Reais/PontoFixo.m`](Raizes%20Reais/PontoFixo.m) |
| Raízes de um **polinômio** | `roots` nativo | [`EletonsDesemp.m`](EletonsDesemp.m) (exemplo) |
| **Resolver sistema linear** Ax = b (completo, com resultado x) | Gauss + pivoteamento parcial + retrossubstituição | [`EscalonamentoGAUSScomSubst.m`](EscalonamentoGAUSScomSubst.m) ✅ |
| **Diagonalizar** a matriz e ler a solução direto | Gauss-Jordan | [`EscalonamentoGAUSSJORDAM2.m`](EscalonamentoGAUSSJORDAM2.m) ✅ |
| Só **escalonar** (triangular superior), sem resolver | Eliminação de Gauss | [`EscalonamentoGAUSS2.m`](EscalonamentoGAUSS2.m), [`teste2.m`](teste2.m) |
| Calcular **integral definida** | Regra do Trapézio composta | [`Integração/Trapez.m`](Integra%C3%A7%C3%A3o/Trapez.m) |
| Integral com mais precisão | Simpson 1/3 composto | [`Integração/Simpson13composto.m`](Integra%C3%A7%C3%A3o/Simpson13composto.m) |
| Calcular **comprimento de arco** de uma curva | ∫√(1+f'(x)²)dx | [`Integração/`](Integra%C3%A7%C3%A3o/) (ambos fazem isso) |
| **Derivada numérica** de dados tabelados | Diferença progressiva de 3 pontos | [`Diferenciacao/Progressiva3Pontos.m`](Diferenciacao/Progressiva3Pontos.m) |
| Resolver **EDO** y' = f(x, y) com valor inicial | Euler explícito | [`EDO/edoEULER.m`](EDO/edoEULER.m) |
| Resolver **sistema de 2 EDOs** acopladas | Runge-Kutta 2ª ordem (Heun) | [`EDO/Sys2EDOsRK2.m`](EDO/Sys2EDOsRK2.m) |
| Runge-Kutta 2 para uma EDO genérica (com `feval`) | RK2 | [`testeRK.m`](testeRK.m) |
| Aproximar **cos(x)** por série | Série de Taylor/Maclaurin | [`Coseno.m`](Coseno.m) |
| **Plotar** uma função rapidinho | `linspace` + `plot` | [`Plotar.m`](Plotar.m), [`QuestaoD3M3_EmAula.m`](QuestaoD3M3_EmAula.m) |
| Gerar **onda senoidal / quadrada / triangular** com média e RMS | `sin`, `square`, `sawtooth`, `rms` | [`Octave/func_generator.m`](Octave/func_generator.m) |
| Calcular **energia reticular** de cristal iônico | Equação de Born-Mayer | [`Eq_BornMayer.m`](Eq_BornMayer.m) |
| **Susceptibilidade magnética** a partir de nº de elétrons desemparelhados | μ = √(n(n+2)) | [`SuscepMag.m`](SuscepMag.m) |
| **Nº de elétrons desemparelhados** a partir da susceptibilidade | inverso do anterior | [`EletonsDesemp.m`](EletonsDesemp.m) |
| Simular **concentração de remédio no sangue** com doses repetidas | Modelo de meia-vida, passo de 1 h | [`pessoal/`](pessoal/) |
| Converter **kg → lb** e **m → polegadas** | Função simples | [`my_functions/SItoING.m`](my_functions/SItoING.m) |
| Matriz aleatória, **transposta e inversa** | `randi`, `'`, `inv` | [`L1Q2.m`](L1Q2.m) |

---

## 🗺️ Mapa do repositório

```
MetodosNumericos/
│
├── Raizes Reais/            ← 🎯 raízes de f(x)=0 (versões mais novas)
│   ├── Bissecao.m              f(x) = x² − e^(−x)            em [0, 1]
│   ├── PosFalsa.m              f(x) = x² − e^(−x)            em [0, 1]
│   ├── Newton.m                f(x) = x³ + 3.8x² − 8.6x − 24.4,  x0 = 4
│   ├── Secante.m               f(x) = 2x³ − 11.7x² + 17.7x − 5,  x0=0, x1=1
│   └── PontoFixo.m             f(x) = eˣ − 4x²,  g(x) = 0.5·e^(x/2)
│
├── Integração/              ← ∫ integrais definidas
│   ├── Trapez.m                Trapézio composto (n = 4000)
│   └── Simpson13composto.m     Simpson 1/3 composto (n = 20)
│
├── Diferenciacao/
│   └── Progressiva3Pontos.m    derivada por diferença progressiva O(h²)
│
├── EDO/                     ← 📈 equações diferenciais ordinárias
│   ├── edoEULER.m              Euler explícito
│   ├── Sys2EDOsRK2.m           sistema de 2 EDOs por RK2
│   └── edoEULERimpNewton.m     (INCOMPLETO) ideia: Euler implícito + Newton
│
├── Octave/                  ← 〰️ geração de sinais
│   ├── func_generator.m        senoide + quadrada + triangular, média e RMS
│   ├── Gerador de função.m     senoide simples
│   ├── rmssaw.m                RMS de dente-de-serra (com bug, ver abaixo)
│   └── 01                      workspace salvo do Octave (dados de uma senoide)
│
├── my_functions/
│   └── SItoING.m               conversão SI → imperial (versão CORRETA)
│
├── pessoal/                 ← 💊 farmacocinética
│   ├── ConcentPlasm_SemAbs.m   doses instantâneas (bolus)
│   ├── ConcentPlasm_ComAbs.m   com taxa de absorção
│   └── intXConc.m              varredura: intervalo/meia-vida × concentração
│
│   ── Arquivos soltos na raiz ──────────────────────────────────────
├── Bissecao.m  PosFalsa.m  Newton.m  Secante.m  PontoFixo.m
│                               versões antigas dos métodos de raiz
│                               (quase todos com f(x) = 8 − 4.5(x − sin x))
├── EscalonamentoGAUSScomSubst.m  ✅ Gauss completo (melhor versão)
├── EscalonamentoGAUSS.m          Gauss + retrossubstituição (com bug)
├── EscalonamentoGAUSS2.m         Gauss na matriz aumentada [A|b]
├── EscalonamentoGAUSSJORDAM.m    apesar do nome, só faz Gauss
├── EscalonamentoGAUSSJORDAM2.m   ✅ Gauss-Jordan de verdade
├── teste.m  teste2.m             rascunhos de Gauss
├── testeRK.m                     RK2 (rascunho de função)
├── Coseno.m                      série de Taylor do cosseno
├── Plotar.m  QuestaoD3M3_EmAula.m  plots rápidos
├── Eq_BornMayer.m  SuscepMag.m  EletonsDesemp.m   química
├── l1q1a.m  l1q1b.m  L1Q2.m      exercícios da Lista 1
└── SItoING.m                     conversão SI → imperial (versão com bug)
```

> 💡 **Duplicatas:** os métodos de raiz existem na raiz **e** em `Raizes Reais/`. O algoritmo é o mesmo; muda a função de exemplo. As versões de `Raizes Reais/` são as mais novas (e o Newton/Secante de lá têm o critério de parada corrigido com `abs(...)`).

---

## ▶️ Como rodar

Os scripts usam **sintaxe do Octave** (`endfor`, `endif`, `printf`, `++`, `+=`, `##`, constante `e`). Use o [GNU Octave](https://octave.org) — no MATLAB vários vão dar erro (veja [Octave × MATLAB](#-octave--matlab)).

```bash
# Linux (Debian/Ubuntu)
sudo apt install octave
```

```bash
# Rodar um script direto do terminal
octave --no-gui "Raizes Reais/Bissecao.m"
```

Ou abra o arquivo na GUI do Octave e aperte **F5**.

**Fluxo padrão de uso** — quase todo script é assim:

1. Começa com `clc; clear all; close all;` (limpa tela, variáveis e gráficos).
2. Tem um bloco de **parâmetros no topo** → **é só isso que você edita** (função `F`, intervalo, chute inicial, tolerância, nº de iterações, matriz…).
3. Roda e lê o resultado no console / gráfico.

**Pacotes** (só para a pasta `Octave/`):

```octave
pkg install -forge signal   % uma vez só (precisa de internet)
pkg load signal             % todo script que usa square/sawtooth
```

> `pkg load nan` também aparece nos scripts, mas não é necessário — pode apagar a linha se der erro.

**Arquivo `Octave/01`** é um *workspace* salvo (variáveis `Senoide`, `t`, `Amplitude_Sinal`, `Frequencia_Sinal`, `Frequencia_Amostragem`, `Valor_Medio`…). Para recuperar:

```octave
load('Octave/01'); whos      % carrega e lista as variáveis
plot(t, Senoide)
```

---

## 🎯 Raízes de equações — f(x) = 0

*Palavras-chave: zero da função, raiz real, equação não linear, convergência, tolerância, iteração.*

### Qual método usar?

| Método | Precisa de | Converge | Sempre converge? | Arquivo |
|---|---|---|---|---|
| **Bisseção** | intervalo [a,b] com f(a)·f(b) < 0 | linear (lento) | ✅ sim | [`Raizes Reais/Bissecao.m`](Raizes%20Reais/Bissecao.m) |
| **Falsa Posição** | intervalo [a,b] com f(a)·f(b) < 0 | superlinear* | ✅ sim | [`Raizes Reais/PosFalsa.m`](Raizes%20Reais/PosFalsa.m) |
| **Newton-Raphson** | chute x₀ **e** derivada f'(x) | quadrática (muito rápido) | ❌ pode divergir | [`Raizes Reais/Newton.m`](Raizes%20Reais/Newton.m) |
| **Secante** | dois chutes x₀, x₁ | ≈1,618 (rápido) | ❌ pode divergir | [`Raizes Reais/Secante.m`](Raizes%20Reais/Secante.m) |
| **Ponto Fixo** | chute x₀ e uma g(x) com \|g'(x)\| < 1 perto da raiz | linear | ❌ só se \|g'\| < 1 | [`Raizes Reais/PontoFixo.m`](Raizes%20Reais/PontoFixo.m) |

\* pode ficar lenta quando um dos extremos "trava".

**Dica:** antes de tudo, **plote a função** ([`Plotar.m`](Plotar.m)) para ver onde ela cruza o zero e escolher o intervalo / chute.

### Bisseção

```
x_m = (a + b) / 2
se f(a)·f(x_m) < 0  →  b = x_m      senão  →  a = x_m
erro = (b − a) / 2
```

- **Edite:** `F`, `a`, `b`, `Tolerancia`, `iMax`.
- **Saída:** tabela `i | a | b | PMed | F(PMed) | Erro` a cada iteração + `Convergiu para: ...`.
- Se `F(a)*F(b) > 0` o script avisa que não há troca de sinal no intervalo.

### Falsa Posição (Regula Falsi)

```
x_r = (a·f(b) − b·f(a)) / (f(b) − f(a))      ← zero da reta (corda) entre (a,f(a)) e (b,f(b))
```

- Mesma estrutura da bisseção, trocando o ponto médio pelo zero da reta.
- **Edite:** `F`, `a`, `b`, `Tolerancia`, `iMax`.
- ⚠️ Ver [bug](#️-bugs-conhecidos-e-pegadinhas): `Fa` e `Fb` não são recalculados dentro do laço.

### Newton-Raphson

```
x_{i+1} = x_i − f(x_i) / f'(x_i)
para quando |f(x_i)| < tol   ou   |(x_{i+1} − x_i) / x_i| < tol
```

- **Edite:** `F`, `dF` (derivada **calculada à mão**), `x0`, `Tolerancia`, `iMax`.
- **Saída:** `Solução xz = ... obtida em N interações`.

### Secante

```
x_{i+1} = x_i − f(x_i)·(x_{i−1} − x_i) / (f(x_{i−1}) − f(x_i))
```

- Igual Newton, mas troca a derivada por uma inclinação aproximada com os dois últimos pontos. No código: `a` = x_{i−1}, `b` = x_i.
- **Edite:** `F`, `a`, `b` (dois chutes), `Tolerancia`, `iMax`.

### Ponto Fixo

```
reescreva f(x) = 0  como  x = g(x)
x_{i+1} = g(x_i)        converge se |g'(x)| < 1 perto da raiz
```

- Exemplo do repo: eˣ − 4x² = 0 → eˣ = 4x² → x = ½·e^(x/2) → `G = @(x)(0.5*exp(0.5*x))`.
- **Edite:** `F` (só para checar), `G`, `xi` (chute), `Tolerancia`, `iMax`.

### Resultados esperados (para conferir)

| Função | Onde | Raiz |
|---|---|---|
| 8 − 4.5(x − sin x) | arquivos da raiz | **x ≈ 2.430466** |
| x² − e^(−x) | `Raizes Reais/Bissecao.m`, `PosFalsa.m` | **x ≈ 0.703467** |
| x³ + 3.8x² − 8.6x − 24.4 | `Raizes Reais/Newton.m` (x0 = 4) | **x ≈ 2.706938** (outras: −2, −4.506938) |
| 2x³ − 11.7x² + 17.7x − 5 | `Raizes Reais/Secante.m` (0 e 1) | **x ≈ 0.365098** (outras: 1.921741, 3.563161) |
| eˣ − 4x² | `PontoFixo.m` | **x ≈ 0.714806** |

### Atalhos nativos do Octave

```octave
fzero(@(x) 8-4.5*(x-sin(x)), [2 3])   % raiz em um intervalo
fzero(@(x) x.^2-exp(-x), 0.5)         % raiz a partir de um chute
roots([1 3.8 -8.6 -24.4])             % TODAS as raízes de um polinômio
```

---

## 🧱 Sistemas lineares — Ax = b

*Palavras-chave: escalonamento, eliminação gaussiana, pivô, pivoteamento parcial, matriz aumentada, retrossubstituição, triangular superior, Gauss-Jordan, matriz diagonal.*

### Qual arquivo?

| Arquivo | O que faz | Resolve? |
|---|---|---|
| [`EscalonamentoGAUSScomSubst.m`](EscalonamentoGAUSScomSubst.m) ✅ | Gauss com **pivoteamento parcial** + **retrossubstituição** | ✅ imprime `x` |
| [`EscalonamentoGAUSSJORDAM2.m`](EscalonamentoGAUSSJORDAM2.m) ✅ | **Gauss-Jordan**: zera abaixo e acima da diagonal, normaliza | ✅ imprime a, R, T |
| [`EscalonamentoGAUSS.m`](EscalonamentoGAUSS.m) | Gauss + retrossubstituição | ⚠️ tem bug em `x(n)` |
| [`EscalonamentoGAUSS2.m`](EscalonamentoGAUSS2.m) | Gauss na matriz aumentada `M = [A \| b]` | só escalona |
| [`EscalonamentoGAUSSJORDAM.m`](EscalonamentoGAUSSJORDAM.m) | igual ao GAUSS2 (o nome engana) | só escalona |
| [`teste2.m`](teste2.m) | Gauss + compara com `rref(A)` | só escalona |
| [`teste.m`](teste.m) | primeira versão, pivô com laço manual | só escalona |

### O algoritmo (Gauss com pivoteamento parcial)

```
para cada coluna k = 1 … n−1:
    1. PIVÔ: acha a linha l (de k até n) com maior |A(l,k)|  →  troca linha l ↔ linha k  (em A e em b)
    2. ELIMINAÇÃO: para cada linha i abaixo de k:
           m = A(i,k) / A(k,k)
           linha_i = linha_i − m · linha_k        (e  b(i) = b(i) − m·b(k))

RETROSSUBSTITUIÇÃO (de baixo para cima):
    x(n) = b(n) / A(n,n)
    x(i) = ( b(i) − Σ_{j>i} A(i,j)·x(j) ) / A(i,i)
```

**Gauss-Jordan** continua: depois de zerar abaixo da diagonal, zera **acima** também (de baixo para cima) e divide cada linha pelo seu pivô. A última coluna de `M` vira a solução.

### Como usar

Edite no topo:

```octave
A = [4 , 3 , 2 ; 1 , 3 , 1 ; 2 , 1 , 3];   % matriz dos coeficientes (; separa linhas)
b = [960 510 610];                          % termos independentes (vetor LINHA)
n = 3;                                      % ordem da matriz
```

### Resultados esperados

| Sistema | Arquivo | Solução |
|---|---|---|
| A = [4 3 2; 1 3 1; 2 1 3], b = [960 510 610] | `EscalonamentoGAUSScomSubst.m` | **x = [120  100  90]** |
| A = [55 0 1; 75 1 −1; 60 −1 0], b = [499.55 687.75 528.60] | `EscalonamentoGAUSSJORDAM2.m` | **a ≈ 9.0311, R ≈ 13.2632, T ≈ 2.8421** |

### Atalhos nativos

```octave
x = A \ b'            % resolve Ax = b (forma recomendada)
rref([A b'])          % forma escalonada reduzida (Gauss-Jordan pronto)
[L, U, P] = lu(A)     % decomposição LU com pivoteamento
inv(A), det(A)        % inversa e determinante
```

---

## ∫ Integração numérica

*Palavras-chave: integral definida, quadratura, área sob a curva, trapézio, Simpson, comprimento de arco, catenária.*

| Método | Fórmula | Erro | Restrição | Arquivo |
|---|---|---|---|---|
| **Trapézio composto** | h·[f(x₀)/2 + f(x₁) + … + f(xₙ₋₁) + f(xₙ)/2] | O(h²) | qualquer n | [`Integração/Trapez.m`](Integra%C3%A7%C3%A3o/Trapez.m) |
| **Simpson 1/3 composto** | (h/3)·[f(x₀) + 4·Σf(ímpares) + 2·Σf(pares) + f(xₙ)] | O(h⁴) | **n par** | [`Integração/Simpson13composto.m`](Integra%C3%A7%C3%A3o/Simpson13composto.m) |

com `h = (b − a) / n`.

- **Edite:** `F`, `LowerLimit`, `UpperLimit`, `n` (nº de subintervalos).
- **Exemplo do repo:** os dois calculam o **comprimento de arco de uma catenária** — L = ∫√(1 + f'(x)²) dx — com f(x) = 211.5 − 20.97·cosh(x/30.38) (curva plotada em [`Plotar.m`](Plotar.m), formato de arco tipo Gateway Arch). O Simpson multiplica por `3.28` no final para converter **metros → pés**.

### Atalhos nativos

```octave
integral(@(x) sqrt(1+0.476454*sinh(x/30.38).^2), -91.21, 91.21)   % função
trapz(x, y)                                                       % dados tabelados
quad(f, a, b)                                                     % alternativa
```

---

## 📐 Diferenciação numérica

*Palavras-chave: derivada numérica, diferenças finitas, progressiva, forward difference, dados experimentais.*

Arquivo: [`Diferenciacao/Progressiva3Pontos.m`](Diferenciacao/Progressiva3Pontos.m)

```
f'(xᵢ) ≈ ( −3·f(xᵢ) + 4·f(xᵢ₊₁) − f(xᵢ₊₂) ) / (2h)        erro O(h²)
```

- Pensado para **dados tabelados**: espera que já existam no workspace os vetores `u` (valores da função) e `y`; passo `h = 0.002`.
- ⚠️ O código escreve `/2*h` (multiplica por h) — o correto é `/(2*h)`. Ver [bugs](#️-bugs-conhecidos-e-pegadinhas).

Outras fórmulas úteis para lembrar:

| Nome | Fórmula | Erro |
|---|---|---|
| Progressiva 2 pontos | (f(x+h) − f(x)) / h | O(h) |
| Regressiva 2 pontos | (f(x) − f(x−h)) / h | O(h) |
| **Central** 2 pontos | (f(x+h) − f(x−h)) / (2h) | O(h²) |
| Regressiva 3 pontos | (3f(x) − 4f(x−h) + f(x−2h)) / (2h) | O(h²) |
| 2ª derivada central | (f(x+h) − 2f(x) + f(x−h)) / h² | O(h²) |

Nativo: `diff(y) ./ diff(x)` ou `gradient(y, h)`.

---

## 📈 EDOs — equações diferenciais

*Palavras-chave: equação diferencial ordinária, problema de valor inicial, PVI, passo de integração, Euler, Runge-Kutta, RK2, Heun, sistema de EDOs.*

| Método | Fórmula | Arquivo |
|---|---|---|
| **Euler explícito** | yᵢ₊₁ = yᵢ + h·f(xᵢ, yᵢ) | [`EDO/edoEULER.m`](EDO/edoEULER.m) |
| **RK2 (Heun) para sistema** de 2 EDOs | k₁ = f(xᵢ, yᵢ); k₂ = f(xᵢ₊₁, yᵢ + h·k₁); yᵢ₊₁ = yᵢ + (h/2)(k₁ + k₂) | [`EDO/Sys2EDOsRK2.m`](EDO/Sys2EDOsRK2.m) |
| **RK2** genérico com `feval` | idem | [`testeRK.m`](testeRK.m) |
| Euler implícito + Newton | — | [`EDO/edoEULERimpNewton.m`](EDO/edoEULERimpNewton.m) ⚠️ **incompleto** (hoje é só Euler explícito) |

### `edoEULER.m`

- **Edite:** `F = @(X,Y)(...)` (a EDO y' = f(x,y)), `X(1)` (x inicial), `B` (x final), `H` (passo), `Y(1)` (condição inicial).
- **Saída:** gráfico de Y(x).
- Exemplo: y' = −2x³ + 12x² − 20x + 8.5, y(0) = 1, de 0 a 4 com h = 0.1 (exemplo clássico do Chapra; solução exata y = −0.5x⁴ + 4x³ − 10x² + 8.5x + 1).

### `Sys2EDOsRK2.m`

- Sistema: y' = EDO1(x, y, z) e z' = EDO2(x, y, z).
- **Edite:** `EDO1`, `EDO2`, `x(1)`, `xn` (x final), `h`, `y(1)`, `z(1)`.
- **Saída:** gráfico de y(x) e z(x) juntos.

### `testeRK.m`

Rascunho de função (o cabeçalho `function` está comentado). Antes de rodar, defina no workspace:

```octave
fun = @(t,y) -2*t*y;   alpha = 1;   a = 0;   b = 2;   N = 20;
testeRK
plot(T, Y)
```

### Atalho nativo

```octave
[x, y] = ode45(@(x,y) -2*x.^3 + 12*x.^2 - 20*x + 8.5, [0 4], 1);
plot(x, y)
```

---

## 🔁 Séries de Taylor

*Palavras-chave: Maclaurin, aproximação de função, cosseno, fatorial, polinômio de Taylor.*

Arquivo: [`Coseno.m`](Coseno.m)

```
cos(x) = Σ (−1)ⁱ · x^(2i) / (2i)!     = 1 − x²/2! + x⁴/4! − x⁶/6! + …
```

- Calcula cos(π/3) com 6 termos (`i = 0:5`) e `format long`. Valor esperado: **0.5**.
- ⚠️ A condição de sinal está errada (`i/2 == 0`) — use `mod(i,2) == 0`. Ver [bugs](#️-bugs-conhecidos-e-pegadinhas).

---

## 〰️ Sinais: senoide, quadrada, triangular, RMS

*Palavras-chave: gerador de função, forma de onda, valor eficaz, valor médio, offset, sawtooth, square, dente de serra, frequência de amostragem.*

Pasta [`Octave/`](Octave/) — precisa de `pkg load signal`.

| Arquivo | O que faz |
|---|---|
| [`func_generator.m`](Octave/func_generator.m) | Gera **senoide** (1 kHz), **quadrada** (15 kHz) e **triangular** (4 kHz), cada uma com offset (valor médio); plota a onda + linha do valor médio + linha do **RMS**, e no fim todas sobrepostas. Imprime as médias. |
| [`Gerador de função.m`](Octave/Gerador%20de%20fun%C3%A7%C3%A3o.m) | Senoide simples. |
| [`rmssaw.m`](Octave/rmssaw.m) | Tentativa de RMS de dente-de-serra (com bug). |
| [`01`](Octave/01) | Workspace salvo com uma senoide (`load('Octave/01')`). |

Receita das ondas:

```octave
pkg load signal
t = linspace(0, 1e-3, 8000);                        % tempo
sen = Vmed + A*sin(2*pi*f*t);                       % senoide
qua = Vmed + A*square(2*pi*f*t);                    % quadrada
tri = Vmed + A*sawtooth(2*pi*f*t, 0.5);             % triangular (0.5 = simétrica)
ser = Vmed + A*sawtooth(2*pi*f*t);                  % dente de serra
rms(sen), mean(sen)                                 % valor eficaz e médio
```

Valores teóricos (sem offset): RMS senoide = A/√2 · RMS quadrada = A · RMS triangular/serra = A/√3. Com offset: RMS = √(Vmed² + RMS_AC²).

---

## ⚗️ Química / Física

| Arquivo | O que calcula | Fórmula |
|---|---|---|
| [`Eq_BornMayer.m`](Eq_BornMayer.m) | **Energia reticular** de cristal iônico (exemplo: NaCl — raios 102 pm + 181 pm, Madelung 1.74756) | U = −(N_A·\|z₊z₋\|·M·e²)/(4πε₀·r₀) · (1 − ρ/r₀) → ≈ **−753 kJ/mol** |
| [`SuscepMag.m`](SuscepMag.m) | Momento magnético *spin-only* e **susceptibilidade molar** a partir de **n** elétrons desemparelhados | μ = √(n(n+2)) ; χ_M = (μ/2.84)² / T |
| [`EletonsDesemp.m`](EletonsDesemp.m) | O inverso: **n elétrons desemparelhados** a partir de χ_M | n² + 2n − 8.0656·χ_M·T = 0 → `roots` |

- **Edite:** cargas `Za`, `Zb`, cte de Madelung `m`, raios `Rio` (pm), `r` (ρ) no Born-Mayer; `n` ou `xM` e `TempK` nos outros.

---

## 💊 Farmacocinética (pasta `pessoal/`)

*Palavras-chave: concentração plasmática, meia-vida, dose, intervalo de administração, absorção, média móvel, estado de equilíbrio.*

Simulação hora a hora (passo de 1 h) da concentração de um fármaco com doses repetidas.

| Arquivo | Modelo |
|---|---|
| [`ConcentPlasm_SemAbs.m`](pessoal/ConcentPlasm_SemAbs.m) | Dose entra **instantaneamente** a cada `Intervalo_Adm_h`. Plota concentração + **média móvel**. |
| [`ConcentPlasm_ComAbs.m`](pessoal/ConcentPlasm_ComAbs.m) | Dose vai para um "depósito" e é **absorvida aos poucos** (`TaxaDeAbs_porhora`). ⚠️ usa `horasmedia` sem definir → erro na penúltima linha. |
| [`intXConc.m`](pessoal/intXConc.m) | Varre 10 intervalos de administração e plota **(intervalo / meia-vida) × (concentração final / dose)**. |

- **Edite:** `Dose`, `Intervalo_Adm_h`, `Meia_Vida_h`, tempo total, janela da média móvel.
- Modelo de eliminação usado: C ← C − C/(2·T½) por hora (aproximação; o decaimento exponencial exato seria C ← C·2^(−1/T½), ≈ C − 0.693·C/T½).

---

## 📝 Exercícios de lista / aula e utilitários

| Arquivo | O que é |
|---|---|
| [`l1q1a.m`](l1q1a.m), [`l1q1b.m`](l1q1b.m) | Lista 1, Questão 1 (a/b): calcular expressões com A=22, B=183, C=108, D=16, E=6 (`power`, `sqrt`). |
| [`L1Q2.m`](L1Q2.m) | Lista 1, Questão 2: matriz 6×6 aleatória (`randi(9,6)`), transposta e inversa. |
| [`QuestaoD3M3_EmAula.m`](QuestaoD3M3_EmAula.m) | Plot de x⁵ − (10/9)x³ + (5/21)x em [−1, 1] (polinômio de Legendre P₅ escalado). |
| [`Plotar.m`](Plotar.m) | Template para plotar qualquer `F` rapidinho (tem algumas funções comentadas para trocar). |
| [`my_functions/SItoING.m`](my_functions/SItoING.m) | `SItoING(massa_kg, altura_m)` → `[libras, polegadas]`. Exemplo de **como escrever uma função** em Octave. |

---

## 📋 Cola rápida de Octave

```octave
%% ---------- Básico ----------
clc; clear all; close all;      % limpa console, variáveis e figuras
format long                     % mostra mais casas decimais (format short volta)
x = 3;                          % ; no fim = não imprime
x                               % sem ; = imprime

%% ---------- Função anônima ----------
F = @(x) x.^2 - exp(-x);        % use .^ .* ./ para funcionar com vetores
F(0.5)                          % avalia
G = @(x,y) x.*y - y;            % várias variáveis

%% ---------- Vetores e matrizes ----------
v = [1 2 3];                    % linha
c = [1; 2; 3];                  % coluna
A = [4 3 2; 1 3 1; 2 1 3];      % ; separa linhas
A'                              % transposta
A(2,3), A(2,:), A(:,1)          % elemento, linha inteira, coluna inteira
A(k:n, k)                       % pedaço da coluna k
M = horzcat(A, b')              % matriz aumentada [A | b]  (= [A b'])
[valor, pos] = max(abs(v))      % maior valor e posição
zeros(1,n), ones(n), eye(n)     % vetor de zeros, matriz de uns, identidade
length(v), size(A)              % tamanhos

%% ---------- Intervalos ----------
x = linspace(-1, 1, 50);        % 50 pontos de -1 a 1
x = 0:0.1:4;                    % de 0 a 4 com passo 0.1
for i = n-1 : -1 : 1            % contagem regressiva

%% ---------- Controle (sintaxe Octave) ----------
for i = 1:10
  if abs(erro) < tol
    break
  elseif cond
    ...
  else
    ...
  endif
endfor

%% ---------- Saída formatada ----------
printf('%3i %11.6f\n', i, x);   % inteiro com 3 casas, real com 6 decimais
fprintf('x = %.6f\n', x);       % igual (funciona também no MATLAB)
disp(A)                         % imprime sem o nome da variável

%% ---------- Gráficos ----------
figure(1);
plot(x, F(x), 'r', 'linewidth', 2);
hold on; plot(x, G(x), 'b'); hold off;
grid on; xlabel('x'); ylabel('F(x)'); title('...');
legend('F', 'G');

%% ---------- Função em arquivo ----------
% arquivo MinhaFuncao.m (nome do arquivo = nome da função)
function [s1, s2] = MinhaFuncao(e1, e2)
  s1 = e1 * 2;
  s2 = e2 + 1;
endfunction
% uso: addpath('my_functions'); [a, b] = MinhaFuncao(1, 2)

%% ---------- Prontos do Octave (para conferir resultados) ----------
fzero(F, [a b])        roots([1 0 -2])       A \ b'        rref(M)
integral(F, a, b)      trapz(x, y)           diff(y)./diff(x)
ode45(f, [x0 xf], y0)  factorial(n)          rms(y)  mean(y)
```

---

## 🔄 Octave × MATLAB

Se precisar rodar no **MATLAB**, troque:

| Octave (usado aqui) | MATLAB |
|---|---|
| `endfor`, `endif`, `endfunction` | `end` |
| `printf(...)` | `fprintf(...)` |
| `i++`, `x += 1`, `x -= 1` | `i = i + 1`, `x = x + 1` |
| `## comentário` | `% comentário` |
| `e` (constante de Euler) | `exp(1)` |
| `"texto"` e `'texto'` | prefira `'texto'` |
| `pkg load signal` | apague (Signal Processing Toolbox já vem, se tiver licença) |

---

## ⚠️ Bugs conhecidos e pegadinhas

Coisas encontradas lendo o código — **para lembrar antes de confiar no resultado**:

| Arquivo | Problema | Correção |
|---|---|---|
| `PosFalsa.m`, `Raizes Reais/PosFalsa.m` | `Fa = F(a)` e `Fb = F(b)` são calculados **uma vez, fora do laço**; a fórmula do zero da reta usa valores velhos. Além disso o erro `(b−a)/2` é o da bisseção. | Recalcular `Fa = F(a); Fb = F(b);` no início de cada iteração; usar `abs(ZeroReta − ZeroRetaAnterior)` como erro. |
| `EscalonamentoGAUSS.m` | `x(n) = b(n)/A(n/n)` → `A(n/n)` é `A(1)`. | `x(n) = b(n)/A(n,n)` — ou use `EscalonamentoGAUSScomSubst.m`, que já está certo. |
| `EscalonamentoGAUSSJORDAM.m` | Nome diz Gauss-Jordan mas só faz Gauss (e imprime `i` no meio do laço). | Use `EscalonamentoGAUSSJORDAM2.m`. |
| `EscalonamentoGAUSSJORDAM2.m` | A troca de linhas quando o pivô é zero está quebrada (`a(z,:) = t`, sempre usa a linha 1). Funciona no exemplo porque nenhum pivô é zero. | Implementar o pivoteamento como em `EscalonamentoGAUSScomSubst.m`. |
| `Newton.m`, `Secante.m` (raiz) | Critério `abs(xi-x0)/x0` sem `abs` no denominador → falha para raiz negativa. | Use as versões de `Raizes Reais/`, que têm `abs((xi-x0)/x0)`. |
| Newton / Secante / Ponto Fixo | Se não convergir em `iMax` iterações, `xz` não existe e o `fprintf` final dá erro. | Inicializar `xz = NaN;` antes do laço. |
| `Coseno.m` | `if (i/2 == 0)` só é verdade para i = 0, então os sinais saem errados. | `if mod(i,2) == 0` — ou direto `Res += (-1)^i * X^(2*i)/factorial(2*i)`. |
| `Diferenciacao/Progressiva3Pontos.m` | `/2*h` multiplica por h; usa `u` em vez de `uVector`; precisa de `y` e `u` pré-definidos. | `dx(i) = (-3*u(i) + 4*u(i+1) - u(i+2)) / (2*h);` |
| `SItoING.m` (raiz) | Usa vírgula decimal (`2,205`) → o Octave entende como dois valores. | Use [`my_functions/SItoING.m`](my_functions/SItoING.m) (`2.205`, `39.37`). |
| `Octave/rmssaw.m` | `sawtooth(0.5)` gera um número só, não uma onda; `X` e `F` não são usados. | `Y = A*sawtooth(2*pi*F*X, 0.5);` |
| `Octave/Gerador de função.m` | Senoide de 1000 Hz amostrada a 1000 Hz → **aliasing** (sai tudo ≈ 0). | Amostrar a pelo menos 2× a frequência do sinal (Nyquist), de preferência 10×–20×. |
| `Octave/func_generator.m` | `Frequencia_Amostragem = e^9` (≈ 8103) é na verdade o **nº de pontos**; no fim usa `set(L3, ...)` em vez de `L4`; legenda da triangular força média "0". | Renomear para `N_pontos`; trocar para `L4`; ver comentário no código. |
| `pessoal/ConcentPlasm_ComAbs.m` | `horasmedia` não definida → erro na penúltima linha. | Trocar por `AmostragemMediaMovel_h`. |
| `EDO/edoEULERimpNewton.m` | Incompleto: hoje é Euler **explícito**, e `X(1)` aparece dentro de `F`. | Implementar yᵢ₊₁ = yᵢ + h·f(xᵢ₊₁, yᵢ₊₁) resolvendo com Newton a cada passo. |
| `testeRK.m` | Não roda sozinho: `fun`, `alpha`, `a`, `b`, `N` precisam existir. Só tem RK2, apesar do comentário falar em grau 2/3/4. | Ver exemplo de uso na seção [EDOs](#-edos--equações-diferenciais). |

**Pegadinhas gerais**

- **Acentos estranhos (`fun��o`)?** Os `.m` estão salvos em **ISO-8859-1 (Latin-1)**, não UTF-8. No Octave/VS Code, reabra com a codificação *Latin-1 / Windows-1252*. Para converter de vez:
  ```bash
  iconv -f latin1 -t utf-8 arquivo.m -o arquivo_utf8.m
  ```
- **Arquivo aparece numa linha só?** `Eq_BornMayer.m`, `EDO/edoEULERimpNewton.m` e `Integração/Simpson13composto.m` usam quebra de linha `CR` (Mac antigo). O Octave roda normal, mas alguns editores mostram tudo junto.
- **Esqueceu o ponto?** `x^2` em vez de `x.^2` quebra quando `x` é vetor (no `plot`, por exemplo).
- **`clear all` no topo** apaga tudo do workspace — por isso `testeRK.m` e `Progressiva3Pontos.m` (que dependem de variáveis prévias) **não** têm essa linha.
- **Tolerância relativa vs absoluta:** `abs(f(x)) < tol` mede o quão perto de zero está a função; `abs((x_novo − x)/x) < tol` mede a variação relativa de x. Os scripts usam um **ou** outro (o que vier primeiro).

---

<sub>Feito com GNU Octave 6.4 · Tyrion1606</sub>
