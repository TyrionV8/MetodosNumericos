# Pré-Implementação: Wiki Interativa de Métodos Numéricos

> **IMPORTANTE — DOCUMENTO DE PRÉ-PLANEJAMENTO (NÃO É UM PLANO DE EXECUÇÃO)**  
> Este documento registra as decisões de design, estética aprovada pelo usuário e a arquitetura dinâmica da interface.  
> O agente que assumir esta tarefa **deve obrigatoriamente revisar a branch de refatoração do código** antes de desenhar o plano de execução formal.

---

## 1. Contexto e Objetivos

* **Repositório da Aplicação:** [`/home/davi/repos/MetodosNumericos`](file:///home/davi/repos/MetodosNumericos).
* **Repositório Hospedeiro (Vitrine):** [`/home/davi/repos/HardCoded`](file:///home/davi/repos/HardCoded).
* **Propósito:** Criar uma Wiki viva e interativa que documenta e visualiza os algoritmos de cálculo numérico deste repositório, combinando códigos GNU Octave/MATLAB reais com animações e visualizadores interativos em **Remotion**.
* **Modelo Arquitetural de Integração:** Idêntico ao do **Panteão** ([`docs/PANTEAO-IGNATIUS.md`](file:///home/davi/repos/HardCoded/docs/PANTEAO-IGNATIUS.md)) e do **Mapa Cru** no HardCoded:
  * A aplicação web, o motor Remotion, os visualizadores e os scripts residem **100% neste repositório**.
  * O HardCoded apenas incorpora a versão publicada via `<iframe>` na rota `/apps/metodos-numericos/`, apontando para o Cloudflare Worker independente via Service Binding.
* **Capas Oficiais Padronizadas (`1024 × 1536 px`, proporção 2:3 idêntica ao Panteão):**
  * **Capa Principal Oficial (Variação B2 — Códice Da Vinci + Mandala de Silício):** [`docs/assets/capa.png`](file:///home/davi/repos/MetodosNumericos/docs/assets/capa.png) e [`docs/wiki-pre-implementacao/capa.png`](file:///home/davi/repos/MetodosNumericos/docs/wiki-pre-implementacao/capa.png) (também preservada em [`docs/assets/capa-v2-b-codice-cpu.png`](file:///home/davi/repos/MetodosNumericos/docs/assets/capa-v2-b-codice-cpu.png)).
  * **Capa Alternativa Salva (Variação AB — Síntese Arcana Digital):** [`docs/assets/capa-v2-ab-sintese-arcana.png`](file:///home/davi/repos/MetodosNumericos/docs/assets/capa-v2-ab-sintese-arcana.png) e [`docs/wiki-pre-implementacao/capa-ab-sintese-arcana.png`](file:///home/davi/repos/MetodosNumericos/docs/wiki-pre-implementacao/capa-ab-sintese-arcana.png).

---

## 2. Estética Aprovada pelo Usuário (Diretrizes Inegociáveis)

O usuário testou e **aprovou com entusiasmo** a identidade visual terminal/TUI desenvolvida no arquivo de referência:
👉 [`referencia-terminal.html`](file:///home/davi/repos/MetodosNumericos/docs/wiki-pre-implementacao/referencia-terminal.html) *(abra no navegador para inspecionar visual e comportamento)*.

### Pilares de Design Aprovados:
1. **Atmosfera Terminal TUI / Workstation Curses:**
   * Fundo breu de console (`#030508`), painéis em `#060b12` com bordas táticas em `#1a324d`.
   * Scanlines CRT sutis (`opacity: 0.35`, passo de 6px) que conferem textura analógica sem comprometer a leitura.
   * Paleta de alto contraste com efeito fósforo:
     * **Verde Elétrico (`#00ff66`):** Status nominais, nós ativos, alvos atingidos e opcodes.
     * **Ciano Fóssil (`#00f0ff`):** Curvas analíticas, retículos e variáveis selecionadas.
     * **Âmbar Técnico (`#ffb700`):** Pontos em processamento intermediário e títulos de seções.
     * **Carmim de Alerta (`#ff3366`):** Vetores tangentes e avisos de singularidade.
2. **Escala de Tipografia Ampla e Legível (Fontes Grandes):**
   * **Problema superado:** Versões com fontes pequenas (10–12px) causavam fadiga de leitura.
   * **Escala definitiva obrigatória:**
     * Fonte base monoespelhada (`JetBrains Mono`): `17px`.
     * Título principal do método: `24px bold` com chip identificador em `13px bold`.
     * Fórmula matemática em destaque: `20px bold` com subtexto em `14px`.
     * Telemetria do HUD ($x_k, \Delta x$, resíduos): `18px` a `20px bold` com brilho neon.
     * Linha de código Octave no HUD: `15px bold`.
     * Árvore do menu lateral: nós em `16px`, categorias em `15px bold`.
     * Script Octave e Tabela de Streaming: `15px` a `16px`.
     * Textos internos do gráfico vetorial: `15px bold`.
3. **Zero Prosa Vazia (Alta Densidade de Informação):**
   * Nada de redação longa, dissertativa ou prolixa.
   * O usuário quer foco absoluto em **dados, equações matemáticas compactas, código Octave real e visualização animada**.
4. **Osciloscópio Vetorial Remotion + HUD Sincronizado:**
   * O player simula um display CRT vetorial calibrado com retículo tático.
   * HUD lateral dinâmico exibindo variáveis em tempo real e a linha do script Octave sendo avaliada naquele passo.
5. **Divisão Inferior em Duas Abas (Split Screen):**
   * **Lado Esquerdo:** O script Octave real exibido com syntax highlighting terminal.
   * **Lado Direito:** Tabela de streaming numérico de convergência com a linha ativa iluminada em verde conforme o vídeo/passo avança.
6. **Navegação Rápida por Teclado (Hotkeys de Terminal):**
   * <kbd>SPACE</kbd>: Pausar / Retomar animação.
   * <kbd>[</kbd> e <kbd>]</kbd> (ou setas): Passo anterior / próximo passo.
   * <kbd>R</kbd>: Resetar para o chute inicial $x_0$.
   * <kbd>1</kbd> a <kbd>9</kbd>: Seleção direta de rotinas.

---

## 3. Arquitetura Dinâmica e Desacoplada (Decoupled Engine)

> [!IMPORTANT]
> **Atenção ao estado do código:**  
> Outro agente está trabalhando em uma branch paralela refatorando e reorganizando toda a estrutura dos arquivos de código Octave.  
> **Não fixe caminhos estáticos** (como `Raizes Reais/Newton.m` ou `EDO/edoEULER.m`). A arquitetura da wiki deve ser dinâmica e tolerante a reestruturações de pastas.

### Estratégia de Desacoplamento Proposta:

```
┌────────────────────────────────────────────────────────────────────────┐
│                   CÓDIGOS OCTAVE REORGANIZADOS                         │
│   (Localizados na pasta definida pela refatoração, ex: src/, alg/...)   │
└───────────────────────────────────┬────────────────────────────────────┘
                                    │
                         [Script Scanner / Manifest]
                                    │
                                    ▼
┌────────────────────────────────────────────────────────────────────────┐
│                    MANIFESTO DINÂMICO DE MÉTODOS                       │
│                       methods.manifest.json                            │
│  - slug: "newton-raphson"                                              │
│  - visualizerType: "scalar-root"                                       │
│  - sourceFile: resolvedPath (dinâmico via glob/scanner)                 │
│  - mathDef: { formula: "...", order: "O(h²)", stability: "..." }       │
│  - defaultInputs: { x0: 3.0, tol: 1e-6 }                               │
└───────────────────────────────────┬────────────────────────────────────┘
                                    │
                                    ▼
┌────────────────────────────────────────────────────────────────────────┐
│                 FRONTEND WIKI (TUI TERMINAL + REMOTION)                │
│  - Generic Visualizer Registry:                                        │
│    • ScalarRootVisualizer (Newton, Bisseção, Secante, Falsa Posição)    │
│    • MatrixEliminationVisualizer (Gauss, Gauss-Jordan, Pivoteamento)   │
│    • QuadratureVisualizer (Trapézios, Simpson 1/3)                     │
│    • ODEStepperVisualizer (Euler, RK2 / Heun)                          │
│  - Dynamic Code Loader: Lê o arquivo fonte via slug                    │
│  - CRT VectorScope Canvas / Remotion Player Component                  │
└────────────────────────────────────────────────────────────────────────┘
```

1. **Manifesto / Registro de Algoritmos (`methods.manifest.json`):**
   * Em vez do frontend saber o caminho físico de cada arquivo no disco, cada método possui um identificador estável (`slug`: `newton-raphson`, `bissecao`, `gauss-pivoteamento`, etc.).
   * Um script de build ou função de resolução associa o `slug` ao arquivo `.m` real por meio de busca dinâmica ou metadados de cabeçalho dos scripts.
2. **Visualizadores Modulares por Categoria:**
   * Em vez de criar um componente Remotion exclusivo para cada arquivo, criam-se **motores visuais genéricos** parametrizáveis:
     * **Visualizador de Raízes Escalares (`ScalarRootVisualizer`):** Recebe $f(x)$, intervalo ou chute inicial, e o histórico de iterações. Plota a curva, os feixes verticais, as tangentes (Newton), as secantes (Secante/Falsa Posição) ou o corte do intervalo (Bisseção).
     * **Visualizador Matricial (`MatrixEliminationVisualizer`):** Recebe $[A|b]$ e renderiza a matriz como tabuleiro de células com pivoteamento e eliminação em cascata.
     * **Visualizador de Quadratura (`QuadratureVisualizer`):** Fatiamento sob a curva (trapézios vs parábolas de Simpson).
     * **Visualizador de EDOs (`ODEStepperVisualizer`):** Comparação entre curva real e passo vetorial tangente/corretor.

---

## 4. Integração com o HardCoded (Protocolo Panteão)

1. **Build do Repositório MetodosNumericos:**
   * Utilizar Vite/React/Remotion gerando saída estática em `dist/`.
   * Possui seu próprio `wrangler.jsonc` como worker Cloudflare `metodos-numericos`.
2. **Arquivos no HardCoded (quando chegar a fase de integração):**
   * **Catálogo:** `HardCoded/src/content/wiki/metodos-numericos.md`
     ```markdown
     ---
     titulo: Métodos Numéricos — Terminal Interativo
     resumo: Análise numérica rigorosa, rotinas GNU Octave e visualizações dinâmicas em Remotion.
     categoria: Computação Científica
     apresentacao: metodos-numericos
     repositorio: https://github.com/TyrionV8/MetodosNumericos
     tags: [octave, matlab, calculo-numerico, remotion, algebra-linear]
     atualizado: 2026-11-03
     rascunho: false
     ---
     ```
   * **Componente de Apresentação:** `HardCoded/src/components/wiki/metodos-numericos/MetodosNumericos.astro`
     ```astro
     ---
     import Base from '../../../layouts/Base.astro';
     import type { Artigo } from '../../../lib/conteudo';

     interface Props { artigo: Artigo; }
     const { artigo } = Astro.props;
     ---

     <Base titulo={artigo.data.titulo} descricao={artigo.data.resumo} secao="wiki">
       <main class="metodos-wiki" aria-label={artigo.data.titulo}>
         <iframe
           class="metodos-frame"
           src="/apps/metodos-numericos/"
           title="Métodos Numéricos — Visualizador Interativo"
           allow="fullscreen"
         ></iframe>
       </main>
     </Base>

     <style>
       .metodos-wiki { background: #030508; }
       .metodos-frame { display: block; width: 100%; height: calc(100dvh - 61px); min-height: 720px; border: 0; background: #030508; }
     </style>
     ```
   * **Roteamento no Worker do HardCoded:** Encaminhamento do prefixo `/apps/metodos-numericos/*` para o worker `metodos-numericos` via binding.

---

## 5. Instruções para o Próximo Agente (Fase de Planejamento)

Antes de escrever qualquer plano de execução ou alterar código:

1. **Inspecionar a Branch de Refatoração:**
   * Verifique quais branches existem no repositório (`git branch -a`) e identifique onde a reorganização dos códigos Octave está ocorrendo.
   * Mapeie a nova taxonomia de diretórios (onde ficaram as raízes, os sistemas lineares, as integrações e as EDOs).
   * Verifique se as funções Octave têm cabeçalhos/comentários padronizados ou se precisam de um pequeno metadata docstring.
2. **Consultar o Arquivo de Referência:**
   * Abra e inspecione [`referencia-terminal.html`](referencia-terminal.html) para clonar com precisão o CSS, o esquema de cores, as proporções, a tipografia em 17px/24px e o comportamento reativo do player.
3. **Elaborar o Plano de Execução Formal:**
   * Estruture as tarefas divididas em:
     * **Marco 1:** Setup da aplicação frontend (Vite + React + Remotion) no repositório `MetodosNumericos`.
     * **Marco 2:** Criação do script de escaneamento e manifesto dinâmico dos códigos Octave refatorados.
     * **Marco 3:** Implementação dos visualizadores modulares Remotion / Canvas baseados no protótipo.
     * **Marco 4:** Integração dos controles de teclado, telemetria HUD e split view de código/log.
     * **Marco 5:** Configuração do deploy Cloudflare Worker e verificação da ligação com o HardCoded.
