import networkx as nx
import matplotlib.pyplot as plt


def prev_nodes(graph: nx.DiGraph, node: any, direction = True):
    if direction:
        return graph.predecessors(node)
    else:
        return graph.successors(node)
def next_nodes( graph: nx.DiGraph, node, direction = True):
    if direction:
        return graph.successors(node)
    else:
        return graph.predecessors(node)


def generate_test_graph():
    G = nx.DiGraph()
    G.add_nodes_from([1,2,3,4,5,6,7,8,9,10,11,12,13], type=False)
    G.add_nodes_from(["R1", "R2", "R3","R4"], type=True)
    G.add_edges_from([(1,"R1"),(2,"R1"),(3,"R1"),("R1",4),("R1",5),("R1",6)])
    G.add_edges_from([(4,"R2"),(5,"R2"),("R2",7)])
    G.add_edges_from([(7,"R2"),("R2",4),("R2",5)])
    G.add_edges_from([(6, "R3"),(8,"R3"),("R3",10),("R3",9)])
    G.add_edges_from([(12,"R4"),(13,"R4"),("R4",8),("R4",11)])

    return G



def bf_traverse(graph: nx.DiGraph, startNode):
    visited = [startNode]
    queue = [startNode]

    while queue:
        current = queue.pop(0)
        
        for suc in next_nodes(graph, current):
            if suc not in visited:
                visited.append(suc)
                queue.append(suc)
    return graph.subgraph(visited)

def main():
    testGraph = generate_test_graph()
    subgraph = bf_traverse(testGraph, 1)
    nx.draw(subgraph, with_labels=True)
    plt.show()


if __name__ == "__main__":
    main()