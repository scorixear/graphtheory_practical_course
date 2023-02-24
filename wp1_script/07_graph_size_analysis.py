import pickle
import networkx as nx
import os


def get_graph_information(graph: nx.DiGraph, file, prefix) -> None:
    file.write("\t \t ")
    file.write(f"no_nodes {graph.number_of_nodes()}, ")
    nodeTypes = list(dict(graph.nodes(data="nodeType")).values())
    file.write(f"no_reaction_nodes  {nodeTypes.count(1)}, ")
    file.write(f"no_compound_nodes {nodeTypes.count(0)} \n")

def run(
    orig_graph_dir: str = "data_clean/01_crn_clean/",
    pathway_graph_dir: str = "data_clean/02_amino_reaction_cycle_clean/",
    output_dir: str = "data_clean/06_graph_sizes"
):
    
    summary = {}
    # iterate over species
    
    for entry in os.scandir(orig_graph_dir):
        if entry.is_file() and entry.name.endswith(".pi"):
            species, medium = entry.name.split("_")[0:2]
            with open(entry.path,"rb",) as reader:
                original_graph: nx.DiGraph = pickle.load(reader)
                data = get_graph_information(original_graph)
                if species in summary:
                    summary[species][medium]['original'] = data
                else:
                    summary[species] = {
                        medium: {
                            'original': data
                        }
                    }
                
    for entry in os.scandir(orig_graph_dir):
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


if __name__ == "__main__":
    summarize_graph_sizes(
        orig_graph_dir="data_clean/01_crn_clean/",
        pathway_graph_dir="data_clean/02_amino_reaction_cycle_clean/",
        summary_file_path="data_clean/summary_cleaned.txt",
        suffix="_cleaned",
    )
    summarize_graph_sizes(
        orig_graph_dir="data/crn/",
        pathway_graph_dir="data/amino_reaction_cycle/",
        summary_file_path="data/summary.txt",
        suffix="",
    )
