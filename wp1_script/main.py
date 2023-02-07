# main script
import os
import metabolite_subgraph
import pickle


def main():
    datadir = "../../graph_pra"
    for entry in os.scandir(datadir):
        if entry.is_file() and entry.name.endswith(".pi"):
            print(entry.path)
            graph = pickle.load(open(entry.path, "rb"))
            amino_acid_graph = metabolite_subgraph.build_aminoacid_graph(graph)


if __name__ == "__main__":
    main()
