import os
import sys
import pickle
import networkx as nx

sys.path.append("./library")
file_handler = __import__("file_handler")

def run(
    datadir: str = "data/amino_reaction_cycle/",
    resultdir: str = "data/pathway_species/",
):
    # data folder holding the subgraph pickle files
    adam_species = []
    cimIV_species = []
    # for each subgraph
    for entry in os.scandir(datadir):
        if entry.is_file() and entry.name.endswith(".pi"):
            # read in graph and sort them in both mediums
            with open(entry.path, "rb") as reader:
                graph_object = {"graph": pickle.load(reader), "name": entry.name.split(".")[0]}
                if "adam" in entry.name:
                    adam_species.append(graph_object)
                else:
                    cimIV_species.append(graph_object)
    results = []
    # for each medium
    # compare species to every other species
    for index in range(len(adam_species)):
        for other_index in range(index + 1, len(adam_species)):
            results.append(compare(
                    adam_species[index], 
                    adam_species[other_index], 
                    "adam", 
                    resultdir
                ))
    for index in range(len(cimIV_species)):
        for other_index in range(index + 1, len(cimIV_species)):
            results.append(compare(
                    cimIV_species[index], 
                    cimIV_species[other_index], 
                    "cimIV", 
                    resultdir
                ))
    file_handler.write_json(results, resultdir+"results.json")


def difference(graph: nx.DiGraph, other: nx.DiGraph):
    final: nx.DiGraph = graph.copy()
    final.remove_nodes_from(n for n in graph if n in other)
    return final


def compare(graph, other, medium, resultdir):
    # retrieve one-sided difference between the given graph an another
    graph_diff = difference(graph["graph"], other["graph"]).nodes(data=True)
    # and the other side
    other_diff = difference(other["graph"], graph["graph"]).nodes(data=True)
    graph_reaction_nodes = 0
    other_reaction_nodes = 0
    for node in graph["graph"].nodes(data=True):
        if "nodeType" in node[1] and node[1]["nodeType"] == 1:
            graph_reaction_nodes += 1
    for node in other["graph"].nodes(data=True):
        if "nodeType" in node[1] and node[1]["nodeType"] == 1:
            other_reaction_nodes += 1
    final_graph_nodes = []
    final_other_nodes = []
    # go through all nodes, if nodes are Reaction nodes, save meta attribute
    for node in graph_diff:
        if "nodeType" in node[1] and node[1]["nodeType"] == 1:
            final_graph_nodes.append(node[0] + ": " + node[1]["meta"])
    for node in other_diff:
        if "nodeType" in node[1] and node[1]["nodeType"] == 1:
            final_other_nodes.append(node[0] + ": " + node[1]["meta"])
    # define output file path

    outfile = resultdir + medium + "/" + graph["name"] + "_vs_" + other["name"]+ ".json"
    node_results = {
        graph["name"]: final_graph_nodes,
        other["name"]: final_other_nodes
    }

    file_handler.write_json(node_results, outfile)
    # count unique reaction nodes
    results = {
        graph["name"]: {
            "reaction_nodes": graph_reaction_nodes,
            "final_reaction_nodes": final_graph_nodes
        },
        other["name"]: {
            "reaction_nodes": other_reaction_nodes,
            "final_reaction_nodes": final_other_nodes
        }
    }
    return results


if __name__ == "__main__":
    run()
