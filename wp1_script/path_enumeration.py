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

source = "D-glucose"
G = pickle.load(open(
    "C:/Users/franz/graphen_praktikum/graphtheory_practical_course/data/amino_reaction_cycle/blongum_cimIV_aa_cycle.pi", "rb"))
amino_acids = read_file("wp1_script/amino_acids.txt")
acids_present = [aa for aa in amino_acids if G.has_node(aa)]

reactionData = nx.get_node_attributes(G, "nodeType")



# build a graph only consiting of reaction nodes
reactionGraph = nx.DiGraph()

for reactionNode in G:
    #only using reaction nodes
    if reactionData[reactionNode] == 1:
        for product in G.successors(reactionNode):
            for followReaction in G.successors(product):
                reactionGraph.add_edge(reactionNode, followReaction)


# add source and target( glucose and aa)
for reaction in G.successors(source):
    reactionGraph.add_edge(source, reaction)

for a in acids_present:
    for reaction in G.successors(a):
        reactionGraph.add_edge(reaction, a)

#TODO filter subgraphs for each aa
#TODO enumerate simple paths of the subgraphs
#TODO enumerate with branch and bound on the subgraphs
#TODO resore reaction paths from the reaction nodes


# pathsIter = nx.all_simple_paths(reactionGraph, source, "L-proline", cutoff=50)
pathThreshold = 20
paths = []

for path in pathsIter:
    if len(paths) <= pathThreshold:
        paths.append(path)
    else:
        break


