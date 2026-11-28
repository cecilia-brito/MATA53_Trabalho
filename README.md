# Análise Estrutural, Complexidade Algorítmica e Robustez da Rede de Citações `cit-HepTh`

Repositório referente ao Trabalho Prático da disciplina **Teoria dos Grafos (MATA53)** — Instituto de Computação, Universidade Federal da Bahia (UFBA).

* **Autora:** Cecilia Brito Santos (Matrícula: `224116495`)
* **Conjunto de Dados Atribuído:** `cit-HepTh`
* **Fonte Original:** [Stanford Network Analysis Project (SNAP)](https://snap.stanford.edu/data/cit-HepTh.html)

---

## 📋 Sobre o Projeto

Este projeto realiza uma investigação estrutural e algorítmica sobre a rede real de citações científicas **cit-HepTh** (composta originalmente por $27.770$ vértices e $352.807$ arestas)[cite: 27], abrangendo:
1. **Pré-processamento e Tratamento de Dados:** Remoção de *auto-loops* e arestas múltiplas, conversão para grafo não-direcionado simples e extração da Maior Componente Conexa (LCC)[cite: 26, 27].
2. **Análise Estrutural:** Cálculo de métricas topológicas globais, distribuição de graus e visualização de subgrafo reduzido via Busca em Largura (BFS)[cite: 26, 35].
3. **Benchmark de Algoritmos:** Avaliação empírica ($n = 30$ execuções, com média, desvio padrão e Intervalo de Confiança de $95\%$) frente à complexidade assintótica teórica para BFS, DFS, Verificação de Eulerianidade, Dijkstra, Bellman-Ford, Floyd-Warshall, Tarjan (SCC) e Kruskal (AGM)[cite: 26, 35].
4. **Análise Estrutural Avançada:** Verificação manual da propriedade *small-world* (comparada ao modelo equivalente de Erdős-Rényi), análise de distribuição em Lei de Potência (*Power Law*) em escala log-log e avaliação de robustez sob remoção aleatória ($5\%$) versus ataque direcionado por centralidade de grau[cite: 26, 35].

---

## 📂 Estrutura do Repositório

```text
MATA53_Trabalho/
├── .gitignore					   # Filtrar arquivos do versionamento
├── analise_cit_hepth.py           # Script principal para análise e geração de gráficos
├── cit-HepTh_tratado_LCC.edgelist # Grafo tratado (Maior Componente Conexa simples e não-direcionada)
├── Cit-HepTh.txt 				   # Grafo original
└── README.md                      # Documentação e instruções de execução

```


## ⚙️ Pré-requisitos e Instalação

1. Clone este repositório:
```bash
git clone [https://github.com/cecilia-brito/MATA53_Trabalho.git](https://github.com/cecilia-brito/MATA53_Trabalho.git)
cd MATA53_Trabalho

```


2. Instale as bibliotecas necessárias:
```bash
pip install igraph numpy matplotlib scipy

```

---

## 🚀 Como Executar

1. Execute o script principal via terminal:
```bash
python analise_cit_hepth.py

```


3. O programa imprimirá no terminal:
* As estatísticas de pré-processamento e extração da LCC;
* Todas as métricas estruturais;
* A tabela de tempos de execução ($n = 30$) e intervalos de confiança de $95\%$;
* Os resultados manuais de *small-world* e a comparação de robustez.




4. Ao final da execução, os arquivos `cit-HepTh_tratado_LCC.edgelist`, `fig1_distribuicao_graus.png`, `fig2_visualizacao_grafo.png` e `fig3_robustez_comparacao.png` serão salvos automaticamente no diretório atual. A execução completa leva cerca de 5 minutos.


---

## 📊 Resumo dos Principais Resultados

| Métrica / Propriedade | Valor Obtido (LCC) |
| --- | --- |
| **Vértices ($\lvert V \rvert$) / Arestas ($\lvert E \rvert$)** | $27.400$ / $352.021$ |
| **Grau Mínimo / Médio / Máximo** | $1$ / $25{,}6950$ / $2.468$ |
| **Densidade** | $0{,}000938$ |
| **Diâmetro / Raio** | $15$ / $8$ |
| **Comprimento Médio dos Caminhos ($L$)** | $4{,}2787$ |
| **Coeficiente de Clusterização ($C$)** | $0{,}3139$ ($343\times$ maior que o modelo aleatório) |
| **Triângulos** | $1.478.698$ |

```

```