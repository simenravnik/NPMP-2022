from collections import defaultdict, OrderedDict
import argparse

def main(args):
    filename = args.filename

    PATH_IN = args.path_in
    PATH_OUT = args.path_out

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
    with open(PATH_OUT + filename, 'w+') as f:

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


if __name__=='__main__':

    parser=argparse.ArgumentParser()
    parser.add_argument('--path_in', type=str, default='./data/original/')
    parser.add_argument('--path_out', type=str, default='./data/pajek/')
    parser.add_argument('--filename', type=str, default='example.txt')
    args=parser.parse_args()
    main(args)