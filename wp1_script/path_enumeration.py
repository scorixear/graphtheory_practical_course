# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 13:51:54 2023

@author: franz
"""
import networkx as nx
import pickle
import matplotlib.pyplot as plt

def read_file(file: str) -> list[str]:
    """reads in file, splits by new line, trims output."""
    input = open(file, "r")
    lines = input.readlines()
    return [line.strip() for line in lines]

def reverse_bf_search(graph: nx.DiGraph, start_node=any):
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

def build_reaction_graph(G: nx.DiGraph, acids_present: list[str], common_compounds: list[str]) -> nx.DiGraph:
    """create a reaction graph consiting only of reactions, glucose and amino acids.
        We are taking all products (except the common ones) and create edges from them.    
    """
    reactionGraph = nx.DiGraph()
    reactionData = nx.get_node_attributes(G, "nodeType")
    for reactionNode in G:
        #only using reaction nodes
        if reactionData[reactionNode] == 1:
            for product in G.successors(reactionNode):

                #we only want to follow our carbon therefore exclude links via H, ATP ect.
                if product not in common_compounds:
                    continue
                for followReaction in G.successors(product):
                    reactionGraph.add_edge(reactionNode, followReaction)

    # add source and target( glucose and aa)
    for reaction in G.successors(source):
        reactionGraph.add_edge(source, reaction)

    for a in acids_present:
        for reaction in G.successors(a):
            reactionGraph.add_edge(reaction, a)
    
    return reactionGraph

source = "D-glucose"
G = pickle.load(open(
    "C:/Users/franz/graphen_praktikum/graphtheory_practical_course/data/amino_reaction_cycle/ecoli_cimIV_aa_cycle.pi", "rb"))
amino_acids = read_file("amino_acids.txt")
common_compounds = read_file("essential_compounds.txt")
acids_present = [aa for aa in amino_acids if G.has_node(aa)]

reactionGraph = build_reaction_graph(G, acids_present, common_compounds)
print(reactionGraph.size())
#TODO filter subgraphs for each aa


#TODO enumerate simple paths of the subgraphs
#TODO enumerate with branch and bound on the subgraphs
#TODO resore reaction paths from the reaction nodes


# pathsIter = nx.all_simple_paths(reactionGraph, source, "L-proline", cutoff=50)
# pathThreshold = 20
# paths = []

# for path in pathsIter:
#     if len(paths) <= pathThreshold:
#         paths.append(path)
#     else:
#        break


