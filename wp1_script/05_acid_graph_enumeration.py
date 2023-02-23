import networkx as nx
import os
import pickle

metabolite_subgraph = __import__("02_metabolite_subgraph")


def run(datadir: str = "data/crn/", resultdir: str = "data/enumeration/"):
    for entry in os.scandir(datadir):
        if entry.is_file() and entry.name.endswith(".pi"):
            print(entry.path)
            graph: nx.DiGraph = pickle.load(open(entry.path, "rb"))
            essential_compounds = metabolite_subgraph.read_file(
                "wp1_script/essential_compounds.txt"
            )
            ec_clean = [ec for ec in essential_compounds if graph.has_node(ec)]
            essential_compounds = ec_clean
            glucose_graph: nx.DiGraph = metabolite_subgraph.bf_search(
                graph, "D-glucose", essential_compounds
            )

            enumeration_graph = glucose_graph.copy()
            enumeration_graph.remove_nodes_from(n for n in essential_compounds)

            amino_acids = metabolite_subgraph.read_file("wp1_script/amino_acids.txt")
            for acid in amino_acids:
                try:
                    reversed_graph = metabolite_subgraph.reverse_bf_search(
                        enumeration_graph, acid
                    )
                    output_graph = nx.DiGraph()
                    output_graph = nx.compose(reversed_graph, output_graph)
                    print(f"{acid}: {len(reversed_graph.nodes())}")
                    # try:
                    #    print(nx.shortest_path(output_graph, "D-glucose", acid))
                    # except:
                    #    print(f"No path found between D-glucose to {acid}")
                    file_path = f"{resultdir}/{entry.name.split('.')[0]}_{acid}.pi"
                    with open(file_path, "wb") as writer:
                        pickle.dump(output_graph, writer)
                except nx.NetworkXError:
                    print(f"Acid missing: {acid}")


if __name__ == "__main__":
    run()
