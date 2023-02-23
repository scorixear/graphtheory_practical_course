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


def visualize_dict(
    aa_annotation: dict, ax: plt.Axes, title: str = "Amino Acid Heatmap"
):
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
        square=True,
        cbar=False,
        vmin=0,
        vmax=1.0,
    )
    ax.set_title(title)


def run_original_data():
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(16, 9))

    dict_adam_orig = get_aa_dict("data/crn/", "adam")
    visualize_dict(dict_adam_orig, axes.flat[0], "Adam - Original Data")

    dict_cim_orig = get_aa_dict("data/crn/", "cimIV")
    visualize_dict(dict_cim_orig, axes.flat[1], "cimIV - Original Data")

    dict_adam_orig = get_aa_dict("data/amino_reaction_cycle/", "adam")
    visualize_dict(dict_adam_orig, axes.flat[2], "Adam - Amino Acid Synthesis Pathway")

    dict_cim = get_aa_dict("data/amino_reaction_cycle/", "cimIV")
    visualize_dict(dict_cim, axes.flat[3], "cimIV - Amino Acid Synthesis Pathway")

    plt.tight_layout()
    plt.show()


def run_clean_data():
    fig, axes = plt.subplots(figsize=(8, 9), nrows=2, ncols=1)

    dict_adam_orig = get_aa_dict("data_clean/01_crn_clean/", "adam")
    visualize_dict(dict_adam_orig, axes.flat[0], "Adam - Original Data (clean)")

    dict_adam = get_aa_dict("data_clean/02_amino_reaction_cycle_clean/", "adam")
    visualize_dict(dict_adam, axes.flat[1], "Adam - Amino Acid Synthesis Pathway")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    run_original_data()
    run_clean_data()
