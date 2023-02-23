# -*- coding: utf-8 -*-

import networkx as nx
import os
import pickle
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def get_aa_list() -> list[str]:
    amino_acids: list[str] = []
    with open("wp1_script/amino_acids.txt", "r") as file:
        for line in file.readlines():
            amino_acids.append(line.replace("\n", ""))
    return amino_acids


def get_aa_dict(
    amino_reaction_cycle_dir: str = "data/amino_reaction_cycle/", medium: str = "adam"
) -> dict:

    amino_acids = get_aa_list()
    results: dict = dict()

    for entry in os.scandir(amino_reaction_cycle_dir):
        if entry.is_file() and entry.name.endswith(".pi") and (medium in entry.name):
            fname: str = entry.name
            species, file_medium = fname.split("_")[0:2]
            with open(entry.path, "rb") as pickled_graph:
                current_aa_graph: nx.DiGraph = pickle.load(pickled_graph)
            aa_can_synth = np.array(
                [current_aa_graph.has_node(aa) for aa in amino_acids]
            )
            results[(species, medium)] = aa_can_synth
    return results


def visualize_dict(aa_annotation: dict, ax: plt.Axes):
    amino_acids = get_aa_list()
    matrix = np.zeros((len(aa_annotation.keys()), len(amino_acids)))
    species: list[str] = []
    for row_idx, key in enumerate(aa_annotation.keys()):
        matrix[row_idx, :] = aa_annotation[key]
        species.append(key[0])
    sns.heatmap(
        matrix,
        xticklabels=amino_acids,
        yticklabels=species,
        linewidths=0.5,
        ax=ax,
        cbar=False,
        square=True,
    )


def run_original_data():
    dict_adam = get_aa_dict("data/amino_reaction_cycle/", "adam")
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(16, 9))

    mat = visualize_dict(dict_adam, axes[0])
    axes[0].set_title("Adam")
    dict_cim = get_aa_dict("data/amino_reaction_cycle/", "cimIV")
    axes[1].set_title("cimIV")

    mat = visualize_dict(dict_cim, axes[1])
    plt.tight_layout()
    plt.show()


def run_clean_data():
    dict_adam = get_aa_dict("data_clean/02_amino_reaction_cycle_clean", "adam")
    fig, axes = plt.subplots(figsize=(8, 9))

    mat = visualize_dict(dict_adam, axes)
    axes.set_title("Adam")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    run_original_data()
    run_clean_data()
