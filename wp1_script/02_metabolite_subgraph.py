import networkx as nx
import matplotlib.pyplot as plt
import pickle


def generate_test_graph():
    """Generates a Test Graph, described in 'Example_Graph.png'."""
    G = nx.DiGraph()
    # type = False means compound node
    # type = True means reaction node
    # this type annotation is not used later on
    G.add_nodes_from(["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "19"], nodeType=0)
    G.add_nodes_from(["R1", "R3", "R4", "R5"], nodeType=1, reversible=False)
    G.add_node("R2", nodeType=1, reversible=True)
    # direction = 1 is default direction
    # direction = 2 is reversed direction
    G.add_edges_from(
        [("1", "R1"), ("2", "R1"), ("3", "R1"), ("R1", "4"), ("R1", "5"), ("R1", "6")], direction=1, multiplicity=1
    )
    G.add_edges_from([("4", "R2"), ("5", "R2")], direction=1, multiplicity=1)
    G.add_edge("R2", "7", multiplicity=2, direction=1)
    G.add_edges_from([("7", "R2"), ("R2", "4"), ("R2", "5")], direction=2, multiplicity=1)
    G.add_edges_from([("6", "R3"), ("8", "R3"), ("R3", "10"), ("R3", "9")], direction=1, multiplicity=1)
    G.add_edges_from([("12", "R4"), ("13", "R4"), ("R4", "11")], direction=1, multiplicity=1)
    G.add_edge("R4", "8", multiplicity=2, direction=1)
    G.add_edges_from([("11", "R5"), ("19", "R5"), ("R5", "6")], direction=1, multiplicity=1)
    return G


def bf_search(
    graph: nx.DiGraph, start_node: any, essential_compounds: list[any], direction=True
):
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
            add_reactions(
                queue, visited, set(successor_reactions), successor, graph, direction
            )
    # return the subgraph made by all visited nodes
    return graph.subgraph(visited)


def reverse_bf_search(graph: nx.DiGraph, start_node: any):
    """Returns a subgraph composed of all predecessors of a given start graph by performing breath first search"""
    # initialize visited with start node
    visited = set([start_node])
    # and empty queue (will only be holding reaction nodes)
    queue = []

    # retrieve reaction nodes that result in start_node as product
    possible_reactions = prev_nodes(graph, start_node)
    # add reaction nodes to queue
    queue.extend(possible_reactions)
    # possibly extend this graph with additional products of the given reaction

    while queue:
        # select oldest reaction
        current = queue.pop(0)
        # add to visited
        visited.add(current)
        # look at all required educts of reaction
        for predecessor in prev_nodes(graph, current):
            # add educt to visited
            visited.add(predecessor)
            # look at reactions forming educt
            predecessor_reactions = prev_nodes(graph, predecessor)
            # add reaction to queue if reaction was not visited earlier and is not in queue already
            for reaction in predecessor_reactions:
                if reaction not in visited and reaction not in queue:
                    queue.append(reaction)
                # possibly extend graph with additional product from reaction that produces given educt
    # return formed sub graph of all visited educts.
    return graph.subgraph(visited)


def add_reactions(
    queue: list[any],
    visited: set[any],
    reactions: set[any],
    start_node: any,
    graph: nx.DiGraph,
    direction=True,
):
    """Adds given reactions to the queue if
    the educt nodes are already visited
    the reaction is not in visited
    the reaction is not in queue

    queue = The queue the reactions to add on
    visited = The set of visited nodes
    start_node = The node we retrieved the reactions from
    graph = the graph these reactions are in
    direction = True if correct direction, False if reversed"""
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
        direction_label = graph.edges[start_node, reaction]["direction"]
        # initialize missmatch filter
        found_missmatch = False
        # for each educt of given reaction
        for educt in prev_nodes(graph, reaction, direction):
            # retrieve edge label 'direction'
            educt_direction_label = graph.edges[educt, reaction]["direction"]
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


def prev_nodes(graph: nx.DiGraph, node: any, direction=True):
    """Returns list of nodes, that are predecessors of the given node
    Direction = True means correct direction
    Direction = False means reversed direction"""
    if direction:
        return graph.predecessors(node)
    else:
        return graph.successors(node)


def next_nodes(graph: nx.DiGraph, node, direction=True):
    """Returns list of nodes, that are successors of the given node
    Direction = True means correct direction
    Direction = False means reversed direction"""
    if direction:
        return graph.successors(node)
    else:
        return graph.predecessors(node)


def read_file(file: str) -> list[str]:
    """reads in file, splits by new line, trims output."""
    with open(file, "r", encoding="UTF-8") as reader:
        lines = reader.readlines()
    return [line.strip() for line in lines]

def show_graph(graph: nx.Graph):
    """Draws the graph on a plot of MatplotLib"""
    nx.draw(graph, pos=nx.spectral_layout(graph), with_labels=True)
    plt.show()

def build_example_aminoacid(graph: nx.DiGraph) -> nx.DiGraph:
    essential_compounds = ["2","3","19","12","13"]
    amino_acids = ["7","9","10"]
    glucose_graph: nx.DiGraph = bf_search(graph, 1, essential_compounds)
    reverse = nx.DiGraph()
    reachable_amino_acids = []
    for aa in amino_acids:
        if aa in glucose_graph.nodes:
            reachable_amino_acids.append(aa)
    for amino_acid in reachable_amino_acids:
        try:
            current_backward = reverse_bf_search(graph, amino_acid)
            reverse = nx.compose(reverse, current_backward)
        except KeyError:
            print(amino_acid + " was not found in the Digraph")
        except nx.NetworkXError:
            print(amino_acid + " was not found in the Digraph")

    # intersection
    final = glucose_graph.copy()
    final.remove_nodes_from(n for n in glucose_graph if n not in reverse)
    final.remove_edges_from(e for e in glucose_graph.edges if e not in reverse.edges)
    return final


def build_aminoacid_graph(graph: nx.DiGraph) -> nx.DiGraph:
    # read in essential compounds
    essential_compounds = read_file("wp1_script/essential_compounds.txt")
    # read in target amino acid
    amino_acids = read_file("wp1_script/amino_acids.txt")

    #filter essential compounds that are not in the graph
    ec_clean = [ec for ec in essential_compounds if graph.has_node(ec)]
    essential_compounds = ec_clean

    # create subgraph via breadth first search
    glucose_graph: nx.DiGraph = bf_search(graph, "D-glucose", essential_compounds)
    #minus_subgraph: nx.DiGraph = bf_search(graph, essential_compounds[0], essential_compounds)

    #glucose_graph = subgraph.copy()
    #glucose_graph.remove_nodes_from(n[0] for n in subgraph.nodes(data=True) if n[0] in minus_subgraph and n[1]['nodeType']==1)

    # in the second step a reverse search is performed starting from the amino acids
    reverse = nx.DiGraph()

    # some amino acids cannot be reached from glucose.
    reachable_amino_acids = []
    for aa in amino_acids:
        if aa in glucose_graph.nodes:
            reachable_amino_acids.append(aa)

    for amino_acid in reachable_amino_acids:
        #print("\t " + amino_acid)
        try:
            current_backward = reverse_bf_search(graph, amino_acid)
            #print(f"{amino_acid}: {len(current_backward.nodes())}")
            #print(nx.shortest_path(current_backward, "L-proline", "L-histidine"))
            reverse = nx.compose(reverse, current_backward)
            #print(f"{amino_acid}: {len(reverse.nodes())}")
        except KeyError:
            print(amino_acid + " was not found in the Digraph")
            pass
        except nx.NetworkXError:
            print(amino_acid + " was not found in the Digraph")
            pass

    # intersection
    final = glucose_graph.copy()
    final.remove_nodes_from(n for n in glucose_graph if n not in reverse)
    final.remove_edges_from(e for e in glucose_graph.edges if e not in reverse.edges)
    return final
