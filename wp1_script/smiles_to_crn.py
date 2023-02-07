# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 10:06:01 2023

@author: franz

Converts the smiles list into the König representation of a hypergraph
Reaction nodes have nodeType 1 and compounds nodeType 0
edges have a multiplicity attribute that includes the stochiometric coefficient
output is a pickle file
"""

import networkx as nx
from collections import Counter
import matplotlib.pyplot as plt
import numpy
import scipy
import pickle
import os
import sys

# define variables

wd = "C:/Users/franz/graphen_praktikum/graphtheory_practical_course/"
dataDir = "data/smiles_list/"
outdir = "data/crn/"
# infile = "C:/Users/franz/graphen_praktikum/graph_theory/sihumix/acacae_adam/acacae_adam.smiles_list"

os.chdir(wd)
fileList = os.listdir(dataDir)


def parseLines(lines):
    reactionData = dict()

    linesTmp = [line.replace("\n", "") for line in lines]
    lines = linesTmp

    if lines[2] == "None":
        return None

    # parse ID and reverability
    reactionData["bigg"] = lines[0].split("Bigg ID:")[1].split(" ")[1]
    reactionData["meta"] = lines[0].split("MetaNetXId:")[1].split(" ")[1]

    if lines[0].split("Reversible: ")[1].replace("\n", "") == "True":
        reactionData["reversible"] = True
    elif lines[0].split("Reversible: ")[1].replace("\n", "") == "False":
        reactionData["reversible"] = False
    else:
        reactionData["reversible"] = None
    # parse reaction
    reactionData["in"] = [
        line.replace(" ", "") for line in lines[2].split("=")[0].split(" + ")
    ]
    reactionData["out"] = [
        line.replace(" ", "") for line in lines[2].split("=")[1].split(" + ")
    ]

    # collect the multiplicities (stochiometric coeffcents)

    reactionData["in_multiplicity"] = Counter(reactionData["in"])
    reactionData["out_multiplicity"] = Counter(reactionData["out"])

    # parse the smiles
    smilesList = lines[3].replace(">>", ".").split(".")
    compoundList = [
        compound.replace(" ", "")
        for compound in lines[2].replace("=", " + ").split(" + ")
    ]

    smilesDir = dict()
    for compound, smiles in zip(compoundList, smilesList):
        smilesDir[compound] = smiles

    reactionData["smiles"] = smilesDir

    return reactionData


for file in fileList:
    # initialize grpah
    G = nx.DiGraph()

    infile = wd + dataDir + file
    fileName = infile.split("/")[-1].split(".")[0]
    with open(infile, "r") as inData:
        for line in inData:
            if line.startswith("Bigg ID"):
                lines = []
                lines.append(line)

                lines.append(next(inData))
                lines.append(next(inData))
                lines.append(next(inData))

                reactionData = parseLines(lines)

                # add to graph
                if reactionData == None:
                    continue

                # add reaction node
                reactionID = reactionData["bigg"]
                G.add_node(
                    reactionID,
                    nodeType=1,
                    meta=reactionData["meta"],
                    reversible=reactionData["reversible"],
                )

                # add edges to educts and the reverse reaction if possible
                for pre in reactionData["in_multiplicity"].keys():
                    if not G.has_node(pre):
                        G.add_node(pre, nodeType=0, smiles=reactionData["smiles"][pre])
                    G.add_edge(
                        pre,
                        reactionID,
                        multiplicity=reactionData["in_multiplicity"][pre],
                        direction=1,
                    )
                    if reactionData["reversible"]:
                        G.add_edge(
                            reactionID,
                            pre,
                            multiplicity=reactionData["in_multiplicity"][pre],
                            direction=2,
                        )

                # add the products as succsessor nodes to the reaction
                for succ in reactionData["out_multiplicity"].keys():
                    if not G.has_node(pre):
                        G.add_node(
                            succ, nodeType=0, smiles=reactionData["smiles"][succ]
                        )
                    G.add_edge(
                        reactionID,
                        succ,
                        multiplicity=reactionData["in_multiplicity"][succ],
                        direction=1,
                    )
                    if reactionData["reversible"]:
                        G.add_edge(
                            succ,
                            reactionID,
                            multiplicity=reactionData["in_multiplicity"][succ],
                            direction=2,
                        )

    # testing for isolated reactions
    # node = "R_2AGPGAT180"
    # nodes = [node]
    # nodes = nodes + [n for n in G.successors(node)]
    # nodes = nodes + [n for n in G.predecessors(node)
    # nx.draw(G.subgraph(nodes), with_labels=True)
    # plt.show()

    # save output with pickle
    outfile = open(outdir + fileName + "_CRN.pi", "wb")
    pickle.dump(G, outfile)
    outfile.close()

print("Finished creating CRN")
