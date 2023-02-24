import os
import sys
import pickle
import networkx as nx
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

sys.path.append("./library")
file_handler = __import__("file_handler")
metabolite_subgraph = __import__("02_metabolite_subgraph.py")


def get_aa_dict(graph_dir: str, medium: str, amino_acids: list[str]) -> dict:
    results: dict = {}

    for entry in os.scandir(graph_dir):
        if entry.is_file() and entry.name.endswith(".pi") and (medium in entry.name):
            fname: str = entry.name
            species, _ = fname.split("_")[0:2]
            with open(entry.path, "rb") as pickled_graph:
                current_aa_graph: nx.DiGraph = pickle.load(pickled_graph)
            aa_can_synth = np.array(
                [current_aa_graph.has_node(aa) for aa in amino_acids]
            )
            results[species] = aa_can_synth
    return results

def get_glucose_aa_dict(graph_dir: str, medium: str, amino_acids: list[str], essential_compounds: list[str]) -> dict:
    results: dict = {}
    for entry in os.scandir(graph_dir):
        if entry.is_file() and entry.name.endswith(".pi") and (medium in entry.name):
            fname: str = entry.name
            species, _ = fname.split("_")[0:2]
            with open(entry.path, "rb") as reader:
                graph: nx.DiGraph = pickle.load(reader)
            glucose_graph: nx.DiGraph = metabolite_subgraph.bf_search(graph, "D-glucose", essential_compounds)
            aa_can_synth = np.array(
                [glucose_graph.has_node(aa) for aa in amino_acids]
            )
            results[species] = aa_can_synth
    return results

def visualize_dict(
    aa_annotation: dict, title: str, output_file: str
):
    amino_acids = file_handler.read_json("input/amino_acids.json")
    matrix = np.zeros((len(aa_annotation.keys()), len(amino_acids)))
    species: list[str] = []
    for row_idx, key in enumerate(aa_annotation.keys()):
        matrix[row_idx, :] = aa_annotation[key]
        species.append(key[0])
    axes = sns.heatmap(
        matrix,
        xticklabels=amino_acids,
        yticklabels=species,
        linewidths=0.5,
        square=True,
        cbar=False,
        vmin=0,
        vmax=1.0,
    )
    axes.set_title(title)
    plt.savefig(output_file)
    plt.close()


def run(crn_input: str = "data/crn/", amino_input: str = "data/amino_reaction_cycle/", output_folder: str = "data/aa_synthesis/"):
    amino_acids = file_handler.read_json("input/amino_acids.json")
    dict_adam_orig = get_aa_dict(crn_input, "adam", amino_acids)
    visualize_dict(dict_adam_orig, "Adam - Original Data", output_folder+ "adam_original.png")

    dict_cim_orig = get_aa_dict(crn_input, "cimIV", amino_acids)
    visualize_dict(dict_cim_orig, "cimIV - Original Data", output_folder+"cimIV_original.png")

    essential_compounds = file_handler.read_json("input/essential_compounds.json")
    dict_adam_glucose = get_glucose_aa_dict(crn_input, "adam", amino_acids, essential_compounds)
    visualize_dict(dict_adam_glucose, "Adam - Glucose Graph", output_folder+"adam_glucose.png")

    dict_cim_glucose = get_glucose_aa_dict(crn_input, "cimIV", amino_acids, essential_compounds)
    visualize_dict(dict_cim_glucose, "cimIV - Glucose Graph", output_folder+"cimIV_glucose.png")

    dict_adam_orig = get_aa_dict(amino_input, "adam", amino_acids)
    visualize_dict(dict_adam_orig, "Adam - Amino Acid Synthesis Pathway", output_folder+"adam_amino_acid.png")

    dict_cim = get_aa_dict(amino_input, "cimIV", amino_acids)
    visualize_dict(dict_cim, "cimIV - Amino Acid Synthesis Pathway", output_folder+"cimIV_amino_acid.png")

if __name__ == "__main__":
    run()
