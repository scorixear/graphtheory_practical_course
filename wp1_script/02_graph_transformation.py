# main script
import os
import metabolite_subgraph
import pickle
import networkx as nx

def export_example_graph():
    resultdir = "data/example/"
    graph = metabolite_subgraph.generate_test_graph()
    amino_acid_graph = metabolite_subgraph.build_example_aminoacid(graph)
    outfile = resultdir+ "example_CRN.pi"
    with open(outfile, "wb") as writer:
        pickle.dump(amino_acid_graph, writer)
    metabolite_subgraph.show_graph(graph)
    metabolite_subgraph.show_graph(amino_acid_graph)


def main():
    datadir = "data/crn/"
    resultdir = "data/amino_reaction_cycle/"
    for entry in os.scandir(datadir):
        if entry.is_file() and entry.name.endswith(".pi"):
            print(entry.path)
            graph: nx.DiGraph = pickle.load(open(entry.path, "rb"))
            amino_acid_graph = metabolite_subgraph.build_aminoacid_graph(graph)
            outfile = resultdir + entry.name.split("/")[-1].replace("CRN", "aa_cycle")
            pickle.dump(amino_acid_graph, open(outfile, "wb"))
            
            amino_acids = metabolite_subgraph.read_file("wp1_script/amino_acids.txt")
            acids_present = [aa for aa in amino_acids if amino_acid_graph.has_node(aa)]
            acids_not_present = [aa for aa in amino_acids if amino_acid_graph.has_node(aa) == False]
            print(f"Acids present: {len(acids_present)}")
            delimiter = " "
            print(f"Acids not present: {delimiter.join(acids_not_present)}")
            
if __name__ == "__main__":
    main()
