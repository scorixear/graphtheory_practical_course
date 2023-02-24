import networkx as nx
from collections import Counter
import pickle
import os


def parseLines(lines):
    reactionData = dict()

    linesTmp = [line.replace("\n", "") for line in lines]
    lines = linesTmp

    if lines[2] == "None":
        return None

    # parse ID and reverability
    reactionData["bigg"] = lines[0].split("Bigg ID:")[1].split(" ")[1]
    reactionData["meta"] = lines[0].split("MetaNetXId:")[1].split(" ")[1]

    if lines[0].split("Reversible: ")[1].split(" ")[0].replace("\n", "") == "True":
        reactionData["reversible"] = True
    elif lines[0].split("Reversible: ")[1].split(" ")[0].replace("\n", "") == "False":
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

    for v in reactionData["in_multiplicity"].values():
        if v == 0:
            print("error")
    for v in reactionData["out_multiplicity"].values():
        if v == 0:
            print("error")
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


def build_graph_from_file(fname: str) -> nx.DiGraph:
    G = nx.DiGraph()
    with open(fname, "r", encoding="UTF-8") as inData:
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
                    if not G.has_node(succ):
                        G.add_node(
                            succ, nodeType=0, smiles=reactionData["smiles"][succ]
                        )
                    G.add_edge(
                        reactionID,
                        succ,
                        multiplicity=reactionData["out_multiplicity"][succ],
                        direction=1,
                    )
                    if reactionData["reversible"]:
                        G.add_edge(
                            succ,
                            reactionID,
                            multiplicity=reactionData["out_multiplicity"][succ],
                            direction=2,
                        )
    return G


def run(
    smiles_list_directory: str = "input/01_smiles_list/",
    crn_save_directory: str = "data/02_crn_clean/",
):

    for entry in os.scandir(smiles_list_directory):
        if entry.is_file() and entry.name.endswith(".smiles_list"):
            parsed_graph: nx.DiGraph = build_graph_from_file(
                smiles_list_directory + entry.name
            )
            filename = entry.name.split(".")[0]
            with open(crn_save_directory + filename + "_CRN.pi", "wb") as save_file:
                pickle.dump(parsed_graph, save_file)


if __name__ == "__main__":
    run()
