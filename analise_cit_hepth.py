import time
import random
import numpy as np
import igraph as ig
import matplotlib.pyplot as plt

# Semente para reprodutibilidade
random.seed(42)
np.random.seed(42)

print("--- 1. Carregando e Tratando Dados ---")
edges = []
with open("Cit-HepTh.txt", "r") as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        u, v = map(int, line.strip().split())
        edges.append((u, v))

#Mapeando IDs originais
unique_nodes = sorted(list(set([u for e in edges for u in e])))
node_map = {node_id: idx for idx, node_id in enumerate(unique_nodes)}
mapped_edges = [(node_map[u], node_map[v]) for u, v in edges]

#Grafo dirigido original usado apara Tarjan
G_dir = ig.Graph(n=len(unique_nodes), edges=mapped_edges, directed=True)
G_dir.simplify(multiple=True, loops=True)

# Grafo não direcionado simples
G_undir = G_dir.as_undirected(mode="collapse")
G_undir.simplify(multiple=True, loops=True)

print(f"Original (simples não-dirigido): V={G_undir.vcount()}, E={G_undir.ecount()}")

#Extraindo componentes conexas
components = G_undir.connected_components()
comp_sizes = sorted([len(c) for c in components], reverse=True)
print(f"Número de componentes conexas: {len(components)}")
print(f"Top 5 tamanhos de componentes: {comp_sizes[:5]}")

#LCC
G = components.giant()
print(f"Maior Componente Conexa (LCC): V={G.vcount()}, E={G.ecount()}")

# grafo tratado
G.write_edgelist("cit-HepTh_tratado_LCC.edgelist")

print("\n--- 2. Parte I: Métricas Estruturais ---")
V = G.vcount()
E = G.ecount()
degrees = G.degree()

print(f"Vértices (LCC): {V}")
print(f"Arestas (LCC): {E}")
print(f"Grau Mínimo: {min(degrees)} | Máximo: {max(degrees)} | Médio: {np.mean(degrees):.4f}")
print(f"Densidade: {G.density():.6f}")
print(f"Coeficiente de Clusterização Médio (Transitivity): {G.transitivity_avglocal_undirected(mode='zero')}")
print(f"Número de Triângulos: {len(G.list_triangles())}")

#diâmetro, raio e caminho médio
print("Calculando distâncias exatas (pode levar ~10s)...")
t0 = time.perf_counter()
avg_path_len = G.average_path_length(directed=False)
diameter = G.diameter(directed=False)
radius = G.radius()
print(f"Comprimento Médio dos Caminhos: {avg_path_len:.4f}")
print(f"Diâmetro: {diameter} | Raio: {radius} (calculado em {time.perf_counter()-t0:.2f}s)")

# gráfico 1: distribuição de graus
deg_values, deg_counts = np.unique(degrees, return_counts=True)
deg_probs = deg_counts / V

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
ax1.bar(deg_values, deg_counts, width=1.0, color="steelblue")
ax1.set_title("Distribuição de Graus (Escala Linear)")
ax1.set_xlabel("Grau (k)")
ax1.set_ylabel("Frequência")
ax1.set_xlim(0, 150)

ax2.loglog(deg_values, deg_probs, 'o', markersize=4, color="darkred", alpha=0.7)
ax2.set_title("Distribuição de Graus (Escala Log-Log - Power Law)")
ax2.set_xlabel("log(k)")
ax2.set_ylabel("log(P(k))")
ax2.grid(True, which="both", ls="--", alpha=0.5)
plt.tight_layout()
plt.savefig("fig1_distribuicao_graus.png", dpi=300)
plt.close()

# gráfico 2: visualização de subgrafo reduzido
sub_nodes = G.bfs(0)[0][:150]
G_vis = G.subgraph(sub_nodes)
layout = G_vis.layout_fruchterman_reingold()
fig, ax = plt.subplots(figsize=(8, 8))
ig.plot(G_vis, target=ax, layout=layout, vertex_size=5, vertex_color="royalblue", edge_width=0.4, edge_color="gray")
ax.set_title("Visualização de Subgrafo Reduzido (150 vértices via BFS)")
plt.savefig("fig2_visualizacao_grafo.png", dpi=300)
plt.close()


print("\n--- 3. Parte II: Benchmark de Algoritmos (n=30) ---")

def benchmark(func, n_runs=30):
    times = []
    for _ in range(n_runs):
        start = time.perf_counter()
        func()
        times.append(time.perf_counter() - start)
    mean = np.mean(times)
    std = np.std(times, ddof=1)
    # Normal Z (n >= 30) -> Z_0.025 = 1.96
    ci = 1.96 * (std / np.sqrt(n_runs))
    return mean, std, ci

# adicionando pesos unitários e subgrafo N=500 para Floyd-Warshall
G.es["weight"] = 1.0
G_sub500 = G.subgraph(G.bfs(0)[0][:500])
G_sub500.es["weight"] = 1.0

# implementação Bellman-Ford 
edges_bf = [(e.source, e.target, 1.0) for e in G.es] + [(e.target, e.source, 1.0) for e in G.es]
def run_bellman_ford(n_v, edge_list, src=0):
    dist = [float('inf')] * n_v
    dist[src] = 0.0
    for _ in range(n_v - 1):
        updated = False
        for u, v, w in edge_list:
            if dist[u] + w < dist[v]:
                dist[v] = dist[u] + w
                updated = True
        if not updated:
            break

# matriz inicial para Floyd-Warshall (N=500)
dist_fw_base = np.full((500, 500), np.inf)
np.fill_diagonal(dist_fw_base, 0.0)
for e in G_sub500.es:
    dist_fw_base[e.source, e.target] = 1.0
    dist_fw_base[e.target, e.source] = 1.0

def run_floyd_warshall(mat_base):
    d = mat_base.copy()
    n = d.shape[0]
    for k in range(n):
        d = np.minimum(d, d[:, k:k+1] + d[k:k+1, :])
    return d

algoritmos = {
    "BFS (1 origem)": lambda: G.bfs(0),
    "DFS (1 origem)": lambda: G.dfs(0),
    "Verificação Eulerianidade": lambda: all(d % 2 == 0 for d in G.degree()),
    "Dijkstra (1 origem)": lambda: G.distances(source=0, weights="weight"),
    "Bellman-Ford (1 origem)": lambda: run_bellman_ford(V, edges_bf, 0),
  	"Floyd-Warshall (Subgrafo N=500)": lambda: run_floyd_warshall(dist_fw_base),
    "Tarjan (SCC no Grafo Dirigido)": lambda: G_dir.connected_components(mode="strong"),
    "Árvore Geradora Mínima (Kruskal)": lambda: G.spanning_tree(weights="weight")
}

print(f"{'Algoritmo':<35} | {'Média (s)':<12} | {'Desvio (s)':<12} | {'IC 95% (±s)':<12}")
print("-" * 78)
for nome, fn in algoritmos.items():
    m, s, ic = benchmark(fn, n_runs=30)
    print(f"{nome:<35} | {m:<12.6f} | {s:<12.6f} | {ic:<12.6f}")

print("\n--- 4. Parte III: Small-World (Implementação Manual) ---")
def manual_small_world_metrics(graph):
    n = graph.vcount()
    adj = [set(graph.neighbors(i)) for i in range(n)]
    c_locals = []
    for i in range(n):
        viz = list(adj[i])
        k = len(viz)
        if k < 2:
            c_locals.append(0.0)
            continue
        links = 0
        for idx1 in range(k):
            u = viz[idx1]
            for idx2 in range(idx1 + 1, k):
                if viz[idx2] in adj[u]:
                    links += 1
        c_locals.append((2.0 * links) / (k * (k - 1)))
    C_obs = sum(c_locals) / n

    sample_sources = random.sample(range(n), min(500, n))
    total_dist = 0
    total_pairs = 0
    for src in sample_sources:
        visited = [-1] * n
        visited[src] = 0
        queue = [src]
        head = 0
        while head < len(queue):
            curr = queue[head]
            head += 1
            d_next = visited[curr] + 1
            for nxt in adj[curr]:
                if visited[nxt] == -1:
                    visited[nxt] = d_next
                    queue.append(nxt)
                    total_dist += d_next
                    total_pairs += 1
    L_obs = total_dist / total_pairs
    return L_obs, C_obs

L_G, C_G = manual_small_world_metrics(G)
# Grafo aleatório equivalente
G_rand = ig.Graph.Erdos_Renyi(n=V, m=E)
L_rand, C_rand = manual_small_world_metrics(G_rand)

print(f"Grafo Real (cit-HepTh): L(G) = {L_G:.4f} | C(G) = {C_G:.4f}")
print(f"Grafo Aleatório Equiv.: L_rand = {L_rand:.4f} | C_rand = {C_rand:.4f}")
print(f"Razão C(G)/C_rand = {C_G/C_rand:.2f}x | Razão L(G)/L_rand = {L_G/L_rand:.2f}x")

print("\n--- 5. Parte III: Robustez (Remoção Aleatória vs Ataque Direcionado) ---")
r = int(0.05 * V)  # 5% dos vértices
print(f"Removendo r = {r} vértices (5% de {V})...")

def avaliar_remocao(graph_base, nodes_to_remove):
    g_temp = graph_base.copy()
    g_temp.delete_vertices(nodes_to_remove)
    n_rem = g_temp.vcount()
    comps = g_temp.connected_components()
    sizes = [len(c) for c in comps]
    
    # Métrica A: Fração da maior componente em relação ao grafo restante
    met_A = max(sizes) / n_rem
    # Métrica B: Número de componentes
    met_B = len(comps)
    # Métrica C: Distância média na maior componente
    giant = comps.giant()
    sample_v = random.sample(range(giant.vcount()), min(100, giant.vcount()))
    dists = giant.distances(source=sample_v)
    met_C = np.mean([d for row in dists for d in row if d > 0])
    # Métrica D: Fração de nós isolados (grau 0)
    met_D = sum(1 for d in g_temp.degree() if d == 0) / n_rem
    return met_A, met_B, met_C, met_D

# remoção aleatória 30 repetições
rand_results = {"A": [], "B": [], "C": [], "D": []}
all_vertices = list(range(V))
for i in range(30):
    rem_nodes = random.sample(all_vertices, r)
    a, b, c, d = avaliar_remocao(G, rem_nodes)
    rand_results["A"].append(a)
    rand_results["B"].append(b)
    rand_results["C"].append(c)
    rand_results["D"].append(d)

# ataque direcionado 
top_r_degree = np.argsort(degrees)[::-1][:r].tolist()
att_A, att_B, att_C, att_D = avaliar_remocao(G, top_r_degree)

print("\nComparação de Robustez (Média ± Desvio Aleatório vs Ataque Direcionado):")
print(f"Métrica A (Maior Comp.) : Aleatório = {np.mean(rand_results['A']):.4f} ± {np.std(rand_results['A']):.4f} | Ataque = {att_A:.4f}")
print(f"Métrica B (Nº Comps.)   : Aleatório = {np.mean(rand_results['B']):.1f} ± {np.std(rand_results['B']):.1f} | Ataque = {att_B}")
print(f"Métrica C (Dist. Média) : Aleatório = {np.mean(rand_results['C']):.4f} ± {np.std(rand_results['C']):.4f} | Ataque = {att_C:.4f}")
print(f"Métrica D (Nós Isolados): Aleatório = {np.mean(rand_results['D']):.4f} ± {np.std(rand_results['D']):.4f} | Ataque = {att_D:.4f}")

# gráfico 3 robustez
fig, axes = plt.subplots(2, 2, figsize=(11, 8))
metricas_info = [
    ("A", "Métrica A: Fração Maior Componente", att_A, axes[0, 0]),
    ("B", "Métrica B: Número de Componentes", att_B, axes[0, 1]),
    ("C", "Métrica C: Distância Média (LCC)", att_C, axes[1, 0]),
    ("D", "Métrica D: Fração de Nós Isolados", att_D, axes[1, 1])
]

for key, titulo, val_ataque, ax in metricas_info:
	ax.boxplot(rand_results[key], tick_labels=["Remoção Aleatória (n=30)"])
	ax.axhline(val_ataque, color="red", linestyle="--", linewidth=2, label=f"Ataque Direcionado ({val_ataque:.3f})")
	ax.set_title(titulo)
	ax.legend()
	ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("fig3_robustez_comparacao.png", dpi=300)
plt.close()
print("\nConcluído! Gráficos salvos como fig1, fig2 e fig3.")