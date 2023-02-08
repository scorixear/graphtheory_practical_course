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


G = pickle.load(open(
    "C:/Users/franz/graphen_praktikum/graphtheory_practical_course/data/amino_reaction_cycle/blongum_cimIV_aa_cycle.pi", "rb"))
amino_acids = read_file("amino_acids.txt")

source = "D-glucose"

# test if amino acid is in graph
acids_present = [aa for aa in amino_acids if G.has_node(aa)]
print("Anzahl vorhandener AA: " + str(len(acids_present)))

# build a graph only consiting of reaction nodes
reactionGraph = nx.DiGraph()

for reactionNode in G:
    if reactionNode.startswith("R_"):
        for succ1 in G.successors(reactionNode):
            for succ2 in G.successors(succ1):
                reactionGraph.add_edge(reactionNode, succ2)

# add source and target( glucose and aa)
for reaction in G.successors(source):
    reactionGraph.add_edge(source, reaction)

for a in acids_present:
    reactionGraph.add_node(a)
    for reaction in G.successors(a):
        reactionGraph.add_edge(reaction, a)

pathsIter = nx.all_simple_paths(reactionGraph, source, "L-proline", cutoff=50)
pathThreshold = 20
paths = []

for path in pathsIter:
    if len(paths) <= pathThreshold:
        paths.append(path)
    else:
        break

pathsArg = []
pathsIter = nx.all_simple_paths(reactionGraph, source, "L-arginine", cutoff=50)

for path in pathsIter:
    if len(pathsArg) <= pathThreshold:
        pathsArg.append(path)
    else:
        break
