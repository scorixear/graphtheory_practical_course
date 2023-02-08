import os
import pickle
import networkx as nx

def main():
    # data folder holding the subgraph pickle files
    datadir = "data/amino_reaction_cycle/"
    adam_species = []
    cimIV_species = []
    # for each subgraph
    for entry in os.scandir(datadir):
        if entry.is_file() and entry.name.endswith(".pi"):
            # read in graph and sort them in both mediums
            with open(entry.path, "rb") as reader:
                graph_object = {'graph': pickle.load(reader), 'name': entry.name}
                if "adam" in entry.name:
                    adam_species.append(graph_object)
                else:
                    cimIV_species.append(graph_object)
    result_str = ""
    # for each medium
    # compare species to every other species
    for index in range(len(adam_species)):
        for other_index in range(index+1, len(adam_species)):
            result_str = result_str + "\n" + compare(adam_species[index], adam_species[other_index], "adam")
    for index in range(len(cimIV_species)):
        for other_index in range(index+1, len(cimIV_species)):
            result_str = result_str + "\n" + compare(cimIV_species[index], cimIV_species[other_index], "cimIV")
    with open("data/pathway_species/results.txt", "w", encoding="UTF-8") as writer:
        writer.write(result_str)

def difference(graph: nx.DiGraph, other: nx.DiGraph):
    final: nx.DiGraph = graph.copy()
    final.remove_nodes_from(n for n in graph if n in other)
    return final

def compare(graph, other, medium):
    # retrieve one-sided difference between the given graph an another
    graph_diff = difference(graph['graph'], other['graph']).nodes(data=True)
    # and the other side
    other_diff = difference(other['graph'], graph['graph']).nodes(data=True)
    final_graph_nodes = set()
    final_other_nodes = set()
    # go through all nodes, if nodes are Reaction nodes, save meta attribute
    for node in graph_diff:
        if 'nodeType' in node[1] and node[1]['nodeType'] == 1:
            final_graph_nodes.add(node[1]['meta'])
    for node in other_diff:
        if 'nodeType' in node[1] and node[1]['nodeType'] == 1:
            final_other_nodes.add(node[1]['meta'])
    # define output file path
    outfile = "data/pathway_species/"+medium+"/"+graph['name']+"_vs_"+other['name']
    # write metanames to file
    outstring = graph['name']+":\n"+"\n".join(final_graph_nodes)+"\n\n"+other['name']+":\n"+"\n".join(final_other_nodes)
    with open(outfile, "w", encoding="UTF-8") as writer:
        writer.write(outstring)
    # count unique reaction nodes
    return_str = f"{graph['name']} vs {other['name']}: {len(final_graph_nodes)} - {len(final_other_nodes)}"
    print(return_str)
    return return_str

if __name__ == "__main__":
    main()
    