import os
import sys
import pickle
import networkx as nx

sys.path.append("./library")
file_handler = __import__("file_handler")
metabolite_subgraph = __import__("02_metabolite_subgraph")


def run(datadir: str = "data/crn/", resultdir: str = "data/enumeration/"):
    for entry in os.scandir(datadir):
        if entry.is_file() and entry.name.endswith(".pi"):
            #print(entry.path)
            graph: nx.DiGraph = pickle.load(open(entry.path, "rb"))
            essential_compounds = file_handler.read_json(
                "input/essential_compounds.json"
            )
            ec_clean = [ec for ec in essential_compounds if graph.has_node(ec)]
            essential_compounds = ec_clean
            glucose_graph: nx.DiGraph = metabolite_subgraph.bf_search(
                graph, "D-glucose", essential_compounds
            )

            enumeration_graph = glucose_graph.copy()
            enumeration_graph.remove_nodes_from(n for n in essential_compounds)

            amino_acids = file_handler.read_json("input/amino_acids.json")
            for acid in amino_acids:
                try:
                    reversed_graph = metabolite_subgraph.reverse_bf_search(
                        enumeration_graph, acid
                    )
                    output_graph = nx.DiGraph()
                    output_graph = nx.compose(reversed_graph, output_graph)
                    # print(f"{acid}: {len(reversed_graph.nodes())}")
                    file_path = f"{resultdir}/{entry.name.split('.')[0]}_{acid}.pi"
                    with open(file_path, "wb") as writer:
                        pickle.dump(output_graph, writer)
                except nx.NetworkXError:
                    pass
                    # print(f"Acid missing: {acid}")


if __name__ == "__main__":
    run()
