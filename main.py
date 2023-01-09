import networkx as nx
import matplotlib.pyplot as plt

G = nx.read_pajek("data/pajek/example.txt")

print(G.number_of_nodes())
print(G.number_of_edges())

print(G.edges())

pos = nx.spring_layout(G)
nx.draw_networkx_nodes(G, pos, cmap=plt.get_cmap('jet'), node_size = 500)
nx.draw_networkx_labels(G, pos)
nx.draw_networkx_edges(G, pos, edgelist=G.edges(), arrows=True)
plt.savefig("./img/example.png", format="PNG", dpi=300)
plt.show()