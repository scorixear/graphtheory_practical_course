# main script
import os
import metabolite_subgraph
import pickle


def main():
    datadir = "../data/crn/"
    resultdir = "../data/amino_reaction_cycle/"
    for entry in os.scandir(datadir):
        if entry.is_file() and entry.name.endswith(".pi"):
            print(entry.path)
            graph = pickle.load(open(entry.path, "rb"))
            amino_acid_graph = metabolite_subgraph.build_aminoacid_graph(graph)
            
            outfile = resultdir + entry.name.split("/")[-1].replace("CRN", "aa_cycle")
            pickle.dump(amino_acid_graph, open(outfile, "wb"))
            
if __name__ == "__main__":
    main()
