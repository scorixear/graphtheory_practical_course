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
    G.add_edges_from([(1,"R1"),(2,"R1"),(3,"R1"),("R1",4),("R1",5),("R1",6)], direction=1)
    G.add_edges_from([(4,"R2"),(5,"R2"),("R2",7)], direction=1)
    G.add_edges_from([(7,"R2"),("R2",4),("R2",5)], direction=2)
    G.add_edges_from([(6, "R3"),(8,"R3"),("R3",10),("R3",9)], direction=1)
    G.add_edges_from([(12,"R4"),(13,"R4"),("R4",8),("R4",11)], direction=1)

    return G

def add_reactions(queue: list[any], visited: set[any], reactions: set[any], startingNode: any, graph: nx.DiGraph, direction=True):
    for reaction in reactions:
        direction_label =graph.edges[startingNode, reaction]['direction']
        conditions = [c for c in prev_nodes(graph, reaction, direction) if graph.edges[c, reaction]['direction'] == direction_label]

        if set(conditions).issubset(visited) and reaction not in visited:
            queue.append(reaction)

def bf_traverse(graph: nx.DiGraph, startNode, essentialCompounds, direction=True):
    # mark all essential nodes and starting node as visited
    visited = set(essentialCompounds)
    visited.add(startNode)
    print(visited)
    queue = []

    # initialize queue
    for startingNode in visited:
        possible_reactions = next_nodes(graph, startingNode, direction)
        add_reactions(queue, visited, possible_reactions, startingNode, graph, direction)

    while queue:
        # dequeue current hyper edge (reaction)
        current = queue.pop(0)
        # add current hyper edge to visited
        visited.add(current)
        #enqueue
        for successor in next_nodes(graph, current, direction):
            visited.add(successor)
            # find all hyperedges of the products
            successor_reactions = next_nodes(graph, successor, direction)
            # add hyperedge to possible reactions
            add_reactions(queue, visited, set(successor_reactions), successor, graph, direction)
    return graph.subgraph(visited)

def main():
    testGraph = generate_test_graph()
    essential_compounds = [2,3]
    subgraph = bf_traverse(testGraph, 1, essential_compounds)
    
    nx.draw(subgraph, with_labels=True)
    plt.show()


if __name__ == "__main__":
    main()