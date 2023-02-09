import networkx as nx
import metabolite_subgraph
import os
import pickle

def main():
    datadir = "data/crn/"
    resultdir = "data/enumeration/"
    for entry in os.scandir(datadir):
        if entry.is_file() and entry.name.endswith(".pi"):
            print(entry.path)
            graph: nx.DiGraph = pickle.load(open(entry.path, "rb"))
            essential_compounds = metabolite_subgraph.read_file("wp1_script/essential_compounds.txt")
            glucose_graph: nx.DiGraph = metabolite_subgraph.bf_search(graph, 'D-glucose', essential_compounds)

            #minus_graph = metabolite_subgraph.bf_search(graph, essential_compounds[0], essential_compounds)
            enumeration_graph = glucose_graph.copy()
            enumeration_graph.remove_nodes_from(n for n in essential_compounds)

            amino_acids = metabolite_subgraph.read_file("wp1_script/amino_acids.txt")
            for acid in amino_acids:
                try:
                    reversed = metabolite_subgraph.reverse_bf_search(enumeration_graph, acid)
                    output_graph = nx.DiGraph()
                    output_graph = nx.compose(reversed, output_graph)
                    print(f"{acid}: {len(reversed.nodes())}")
                    #try:
                    #    print(nx.shortest_path(output_graph, "D-glucose", acid))
                    #except:
                    #    print(f"No path found between D-glucose to {acid}")
                    file_path = f"{resultdir}/{entry.name.split('.')[0]}_{acid}.pi"
                    with open(file_path,"wb") as writer:
                        pickle.dump(output_graph, writer)
                except nx.NetworkXError:
                    print(f"Acid missing: {acid}")
                
            


if __name__ == "__main__":
    main()