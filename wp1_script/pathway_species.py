import os
import pickle
import networkx as nx

def main():
    datadir = "data/amino_reaction_cycle/"
    adam_species = []
    cimIV_species = []
    for entry in os.scandir(datadir):
        if entry.is_file() and entry.name.endswith(".pi"):
            if "adam" in entry.name:
                adam_species.append(pickle.load(open(entry.path, "rb")))
            else:
                cimIV_species.append(pickle.load(open(entry.path, "rb")))
    for specie in adam_species:
        for other_specie in adam_species:
            if specie != other_specie:
                compare(specie, other_specie)

def difference(graph: nx.DiGraph, other: nx.DiGraph):
    final: nx.DiGraph = graph.copy()
    final.remove_nodes_from(n for n in graph if n in other)
    return final

def compare(graph, other):
    

if __name__ == "__main__":
    main()