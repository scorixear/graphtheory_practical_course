import os
import sys
import pickle
import networkx as nx

sys.path.append("./input")
file_handler = __import__("file_handler")
metabolite_subgraph = __import__("02_metabolite_subgraph")

def get_graph_information(graph: nx.DiGraph) -> None:
    nodeTypes = list(dict(graph.nodes(data="nodeType")).values())
    return {
        "node_number": graph.number_of_nodes(),
        "reaction_nodes": nodeTypes.count(1),
        "compound_nodes": nodeTypes.count(0)
    }

def run(
    orig_graph_dir: str = "data/crn/",
    pathway_graph_dir: str = "data/amino_reaction_cycle/",
    output_dir: str = "data/graph_sizes/"
):
    
    summary = {}
    # iterate over species
    essential_compounds = file_handler.read_json("input/essential_compounds.json")
    for entry in os.scandir(orig_graph_dir):
        if entry.is_file() and entry.name.endswith(".pi"):
            species, medium = entry.name.split("_")[0:2]
            with open(entry.path,"rb",) as reader:
                original_graph: nx.DiGraph = pickle.load(reader)
                glucose_graph: nx.DiGraph = metabolite_subgraph.bf_search(original_graph, "D-glucose", essential_compounds)
                data = get_graph_information(original_graph)
                glucose_data = get_graph_information(glucose_graph)
                if species in summary:
                    summary[species][medium]['original'] = data
                    summary[species][medium]['glucose'] = glucose_data
                else:
                    summary[species] = {
                        medium: {
                            'original': data,
                            'glucose' : glucose_data
                        }
                    }
                
    for entry in os.scandir(pathway_graph_dir):
        if entry.is_file() and entry.name.endswith(".pi"):
            species, medium = entry.name.split("_")[0:2]
            with open(entry.path,"rb",) as reader:
                original_graph: nx.DiGraph = pickle.load(reader)
                data = get_graph_information(original_graph)
                if species in summary:
                    summary[species][medium]['amino_acid'] = data
                else:
                    summary[species] = {
                        medium: {
                            'amino_acid': data
                        }
                    }
    file_handler.write_json(output_dir+"graph_information.json")

if __name__ == "__main__":
    run()
