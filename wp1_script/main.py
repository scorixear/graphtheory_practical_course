# main script
import os
import metabolite_subgraph
import pickle
import networkx as nx


def main():
    datadir = "data/crn"
    for entry in os.scandir(datadir):
        if entry.is_file() and entry.name.endswith(".pi"):
            print(entry.path)
            graph: nx.DiGraph = pickle.load(open(entry.path, "rb"))
            amino_acid_graph = metabolite_subgraph.build_aminoacid_graph(graph)
            #metabolite_subgraph.show_graph(amino_acid_graph)
            amino_acids = metabolite_subgraph.read_file("wp1_script/amino_acids.txt")
            acids_present = [aa for aa in amino_acids if amino_acid_graph.has_node(aa)]
            print(f"Acids present: {len(acids_present)}")

if __name__ == "__main__":
    main()
