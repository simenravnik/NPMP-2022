from collections import defaultdict, OrderedDict

filename = "example.txt"

PATH_IN = "./data/original/"
PATH_OUT = "./data/pajek/"

with open(PATH_IN + filename) as f:
    lines = f.readlines()

nodes = set()
connections = defaultdict(list)
for line in lines:
    s = line.split()

    from_node = int(s[2])
    to_node = int(s[0])

    # Add both nodes to set of nodes
    nodes.add(from_node)
    nodes.add(to_node)

    # Append to connections
    connections[from_node].append(to_node)

# Order connections by node id inc
connections = dict(OrderedDict(sorted(connections.items())))

# Write to file
with open(PATH_OUT + filename, 'w') as f:

    f.write("*Vertices " + str(len(nodes)) + "\n")

    nodes_ordered = sorted(list(nodes))

    for node in nodes_ordered:
        f.write(str(node) + " " + str(node))
        f.write('\n')

    f.write("*arcs " + "\n")

    for from_node in connections:
        for to_node in connections[from_node]:
            f.write(str(from_node) + " " + str(to_node))
            f.write('\n')