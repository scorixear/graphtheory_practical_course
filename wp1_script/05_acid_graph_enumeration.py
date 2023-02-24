import os
import sys
import pickle
import networkx as nx
import metabolite_subgraph

sys.path.append("./library")
file_handler = __import__("file_handler")

def run(
    datadir: str = "data/crn/", 
    resultdir: str = "data/enumeration/",
    inputdir: str = "input/"):
    for entry in os.scandir(datadir):
        if entry.is_file() and entry.name.endswith(".pi"):
            #print(entry.path)
            graph: nx.DiGraph = pickle.load(open(entry.path, "rb"))
            essential_compounds = file_handler.read_json(
                inputdir+"essential_compounds.json"
            )
            ec_clean = [ec for ec in essential_compounds if graph.has_node(ec)]
            essential_compounds = ec_clean
            glucose_graph: nx.DiGraph = metabolite_subgraph.bf_search(
                graph, "D-glucose", essential_compounds
            )

            enumeration_graph = glucose_graph.copy()
            enumeration_graph.remove_nodes_from(n for n in essential_compounds)

            amino_acids = file_handler.read_json(inputdir+"amino_acids.json")
            for acid in amino_acids:
                try:
                    reversed_graph = metabolite_subgraph.reverse_bf_search(
                        enumeration_graph, acid
                    )
                    output_graph = nx.DiGraph()
                    output_graph = nx.compose(reversed_graph, output_graph)
                    # print(f"{acid}: {len(reversed_graph.nodes())}")
                    species, medium = entry.name.split("/")[-1].split("_")[0:2]
                    file_path = f"{resultdir}/{species}_{medium}_{acid}.pi"
                    with open(file_path, "wb") as writer:
                        pickle.dump(output_graph, writer)
                except nx.NetworkXError:
                    pass
                    # print(f"Acid missing: {acid}")


if __name__ == "__main__":
    run()
