import pickle
import networkx as nx
import os


def write_graph_information(graph: nx.DiGraph, file, prefix) -> None:
    file.write("\t \t ")
    file.write(f"no_nodes {graph.number_of_nodes()}, ")
    nodeTypes = list(dict(graph.nodes(data="nodeType")).values())
    file.write(f"no_reaction_nodes  {nodeTypes.count(1)}, ")
    file.write(f"no_compound_nodes {nodeTypes.count(0)} \n")


def get_species_media(path: str) -> tuple[list[str], list[str]]:
    species: list[str] = []
    media: list[str] = []
    for entry in os.scandir(path):
        if entry.is_file() and entry.name.endswith(".pi"):
            species.append(entry.name.split("_")[0])
            media.append(entry.name.split("_")[1])
    return list(set(species)), list(set(media))


def summarize_graph_sizes(
    orig_graph_dir: str = "data_clean/01_crn_clean/",
    pathway_graph_dir: str = "data_clean/02_amino_reaction_cycle_clean/",
    summary_file_path: str = "data_clean/summary.txt",
    suffix: str = "cleaned",
):
    with open(summary_file_path, "w+") as summary_file:
        # iterate over species
        all_species, all_media = get_species_media(orig_graph_dir)
        for species in all_species:
            summary_file.write(f"{species}------------------\n")
            for medium in all_media:
                # get graph size of the original graph:
                summary_file.write("\t " + medium + "\n")
                with open(
                    orig_graph_dir + species + "_" + medium + suffix + "_CRN.pi",
                    "rb",
                ) as gfile:
                    original_graph: nx.DiGraph = pickle.load(gfile)
                    write_graph_information(
                        original_graph, summary_file, "Original Graph: "
                    )
                with open(
                    pathway_graph_dir
                    + species
                    + "_"
                    + medium
                    + suffix
                    + "_aa_cycle.pi",
                    "rb",
                ) as pathfile:
                    pathway_graph: nx.DiGraph = pickle.load(pathfile)
                    write_graph_information(
                        pathway_graph, summary_file, "Pathway Graph: "
                    )


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
