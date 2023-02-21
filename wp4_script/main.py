import sys
import networkx as nx
import os
import pickle

sys.path.append("./wp1_script")
smiles_to_crn = __import__("01_smiles_to_crn")


# step 1: parse all files, write graphs as pickle files
smiles_list_directory: str = "./data/smiles_list_clean/"
crn_save_directory: str = "./data/crn_clean/"

for entry in os.scandir(smiles_list_directory):
    if entry.is_file() and entry.name.endswith(".smiles_list"):
        parsed_graph: nx.DiGraph = smiles_to_crn.build_graph_from_file(
            smiles_list_directory + entry.name
        )
        filename = entry.name.split(".")[0]
        with open(crn_save_directory + filename + "_CRN.pi", "wb") as save_file:
            pickle.dump(parsed_graph, save_file)
