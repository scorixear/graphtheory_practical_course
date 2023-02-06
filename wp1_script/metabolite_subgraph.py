import networkx as nx
import matplotlib.pyplot as plt
import pickle



def generate_test_graph():
    """Generates a Test Graph, described in 'Example_Graph.png'."""
    G = nx.DiGraph()
    # type = False means compound node
    # type = True means reaction node
    # this type annotation is not used later on
    G.add_nodes_from([1,2,3,4,5,6,7,8,9,10,11,12,13,19], type=False)
    G.add_nodes_from(["R1", "R2", "R3","R4", "R5"], type=True)
    # direction = 1 is default direction
    # direction = 2 is reversed direction
    G.add_edges_from([(1,"R1"),(2,"R1"),(3,"R1"),("R1",4),("R1",5),("R1",6)], direction=1)
    G.add_edges_from([(4,"R2"),(5,"R2"),("R2",7)], direction=1)
    G.add_edges_from([(7,"R2"),("R2",4),("R2",5)], direction=2)
    G.add_edges_from([(6, "R3"),(8,"R3"),("R3",10),("R3",9)], direction=1)
    G.add_edges_from([(12,"R4"),(13,"R4"),("R4",8),("R4",11)], direction=1)
    G.add_edges_from([(11,"R5"),(19,"R5"),("R5",6)], direction=1)
    return G

def bf_search(graph: nx.DiGraph, start_node: any, essential_compounds: list[any], direction=True):
    """Performs Breath First Search modified for Hyper Graphs.
    graph = given Graph to perform BFS on
    start_node = starting node
    essential_compounds = list of compounds assumed as 'given'"""

    # mark all essential nodes and starting node as visited
    # visited is a set, since no node can be visited twice
    # once seen, a node can always be reproduced by that reaction
    visited = set(essential_compounds)
    visited.add(start_node)
    # initialize queue as list (since sets do not keep location)
    queue = []

    # for each compound node in our starting position
    for node in visited:
        # find all possible reactions this compound is part of
        possible_reactions = next_nodes(graph, node, direction)
        # add reaction nodes to the queue
        # if the reaction has only educts, that are in the visited set.
        add_reactions(queue, visited, possible_reactions, node, graph, direction)

    while queue:
        # dequeue current reaction
        current = queue.pop(0)
        # Debug output
        # print(f"Queue: {len(queue)}")
        # print(f"Dequeu {len(visited)}")
        # add current reaction to visited (as in, this reaction is now performed)
        visited.add(current)
        # for each product of this reaction
        for successor in next_nodes(graph, current, direction):
            # add product compound to visited
            visited.add(successor)
            # find all reactions this product is part of as an educt
            successor_reactions = next_nodes(graph, successor, direction)
            # add reaction nodes to the queue
            # if reactions has only educts, that are in the visited set
            # and are not already visited
            # and are not already in queue
            add_reactions(queue, visited, set(successor_reactions), successor, graph, direction)
    # return the subgraph made by all visited nodes
    return graph.subgraph(visited)

def add_reactions(queue: list[any], visited: set[any], reactions: set[any], start_node: any, graph: nx.DiGraph, direction=True):
    """Adds given reactions to the queue if 
    the educt nodes are already visited
    the reaction is not in visited
    the reaction is not in queue
    
    queue = The queue the reactions to add on
    visited = The set of visited nodes
    start_node = The node we retrieved the reactions from
    graph = the graph these reactions are in
    direction = True if correct direciton, False if reversed"""
    # for each given possible reaction of one compound node
    for reaction in reactions:
        # if the current reaction has already been visited 
        #   (can happen for reversable reaction, example Graph Reaction R2)
        # if the current reaction has already been added to the queue 
        #   (can happen if multiple reactions form one product, resulting in same reactions being checked multiple times 
        #    for product as educt, example Graph Reaction R5 -> 6 -> R3)
        if reaction in visited or reaction in queue:
            continue
        # retrieve label "direction" of the edge between the given compound node and the reaction node
        direction_label =graph.edges[start_node, reaction]['direction']
        # initialize missmatch filter
        found_missmatch = False
        # for each educt of given reaction
        for educt in prev_nodes(graph, reaction, direction):
            # retrieve edge label 'direction'
            educt_direction_label = graph.edges[educt, reaction]['direction']
            # if current edge is part of a different direction, ignore it
            if educt_direction_label != direction_label:
                continue
            # if educt is not in visited, reaction cannot be formed (yet)
            if educt not in visited:
                found_missmatch = True
                break
        # of no missmatch was found, add reaction to queue
        if found_missmatch is False:
            queue.append(reaction)


def prev_nodes(graph: nx.DiGraph, node: any, direction = True):
    """Returns list of nodes, that are predecessors of the given node
    Direction = True means correct direction
    Direction = False means reversed direction"""
    if direction:
        return graph.predecessors(node)
    else:
        return graph.successors(node)
def next_nodes( graph: nx.DiGraph, node, direction = True):
    """Returns list of nodes, that are successors of the given node
    Direction = True means correct direction
    Direction = False means reversed direction"""
    if direction:
        return graph.successors(node)
    else:
        return graph.predecessors(node)

def read_file(file: str) -> list[str]:
    """reads in file, splits by new line, trims output."""
    input = open(file,"r")
    lines = input.readlines()
    return [line.strip() for line in lines]



def main():
    # read in / generate graph
    #graph = generate_test_graph()
    graph = pickle.load(open("acacae_adamCRN.pi", "rb"))    
    
    # read in essential compounds
    essential_compounds = read_file("essential_compounds.txt")
    
    # create subgraph via breath first search
    subgraph = bf_search(graph, "D-glucose", essential_compounds)
    #subgraph = bf_search(graph, 1, essential_compounds=[1,2,3,19,12,13])
    
    # generate position dictionary for graph nodes
    pos = nx.spring_layout(subgraph)
    
    # draw network on plot
    nx.draw(subgraph,pos=pos, with_labels=True)
    
    # show plot
    plt.show()


if __name__ == "__main__":
    main()