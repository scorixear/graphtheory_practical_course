import os
import sys
import pickle
import networkx as nx

sys.path.append("./input")
file_handler = __import__("file_handler")

metabolite_subgraph = __import__("02_metabolite_subgraph")


def export_example_graph():
    resultdir = "data/example/"
    graph = metabolite_subgraph.generate_test_graph()
    amino_acid_graph = metabolite_subgraph.build_example_aminoacid(graph)
    outfile = resultdir + "example_CRN.pi"
    with open(outfile, "wb") as writer:
        pickle.dump(amino_acid_graph, writer)
    metabolite_subgraph.show_graph(graph)
    metabolite_subgraph.show_graph(amino_acid_graph)


def run(datadir: str = "data/crn/", resultdir: str = "data/amino_reaction_cycle/"):

    for entry in os.scandir(datadir):
        if entry.is_file() and entry.name.endswith(".pi"):
            print(entry.path)
            graph: nx.DiGraph = pickle.load(open(entry.path, "rb"))
            amino_acid_graph = metabolite_subgraph.build_aminoacid_graph(graph)
            outfile = resultdir + entry.name.split("/")[-1].replace("CRN", "aa_cycle")
            pickle.dump(amino_acid_graph, open(outfile, "wb"))

            amino_acids = file_handler.read_json("input/amino_acids.json")
            acids_present = [aa for aa in amino_acids if amino_acid_graph.has_node(aa)]
            acids_not_present = [
                aa for aa in amino_acids if amino_acid_graph.has_node(aa) is False
            ]
            print(f"Acids present: {len(acids_present)}")
            delimiter = " "
            print(f"Acids not present: {delimiter.join(acids_not_present)}")


if __name__ == "__main__":
    run()
