import sys
import networkx as nx
import os
import pickle

sys.path.append("./wp1_script")
smiles_to_crn = __import__("01_smiles_to_crn")
metabolite_subgraph = __import__("02_metabolite_subgraph")

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

# step 2: analysis scripts

datadir: str = "data/crn_clean/"
resultdir: str = "data/amino_reaction_cycle_clean/"

for entry in os.scandir(datadir):
    if entry.is_file() and entry.name.endswith(".pi"):
        print(entry.path)
        with open(entry.path, "rb") as graph_file:
            graph: nx.DiGraph = pickle.load(graph_file)
        amino_acid_graph = metabolite_subgraph.build_aminoacid_graph(graph)
        outfile = resultdir + entry.name.split("/")[-1].replace("CRN", "aa_cycle")
        with open(outfile, "wb") as save_file:
            pickle.dump(amino_acid_graph, save_file)

        amino_acids = metabolite_subgraph.read_file("./wp1_script/amino_acids.txt")
        acids_present = [aa for aa in amino_acids if amino_acid_graph.has_node(aa)]
        acids_not_present = [
            aa for aa in amino_acids if amino_acid_graph.has_node(aa) == False
        ]
        print(f"Acids present: {len(acids_present)}")
        delimiter = " "
        print(f"Acids not present: {delimiter.join(acids_not_present)}")
