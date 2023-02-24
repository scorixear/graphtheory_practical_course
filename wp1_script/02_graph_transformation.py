import os
import sys
import pickle
import networkx as nx
import metabolite_subgraph

sys.path.append("./library")
file_handler = __import__("file_handler")

def export_example_graph():
    resultdir = "data/example/"
    graph = metabolite_subgraph.generate_test_graph()
    amino_acid_graph = metabolite_subgraph.build_example_aminoacid(graph)
    outfile = resultdir + "example_CRN.pi"
    with open(outfile, "wb") as writer:
        pickle.dump(amino_acid_graph, writer)
    metabolite_subgraph.show_graph(graph)
    metabolite_subgraph.show_graph(amino_acid_graph)


def run(
    datadir: str = "data/crn/", 
    resultdir: str = "data/amino_reaction_cycle/",
    inputdir: str = "input/"):

    for entry in os.scandir(datadir):
        if entry.is_file() and entry.name.endswith(".pi"):
            graph: nx.DiGraph = pickle.load(open(entry.path, "rb"))
            amino_acid_graph = metabolite_subgraph.build_aminoacid_graph(graph, inputdir)
            species, medium = entry.name.split("/")[-1].split("_")[0:2]
            with open(f"{resultdir}{species}_{medium}_aa_cycle.pi", "wb") as writer:
                pickle.dump(amino_acid_graph, writer)


if __name__ == "__main__":
    run()
