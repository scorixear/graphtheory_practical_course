import sys
import networkx as nx
import os
import pickle

sys.path.append("./wp1_script")
sys.path.append("./wp2_script")
sys.path.append("./wp3_script")


smiles_to_crn = __import__("01_smiles_to_crn")

# step 1: parse all files, write graphs as pickle files
smiles_list_directory: str = "data/smiles_list_clean/"
crn_save_directory: str = "data/crn_clean/"

smiles_to_crn.run(smiles_list_directory, crn_save_directory)

# ------------------------------------------------------------------------------
# step 2: analysis scripts

graph_transformation = __import__("02_graph_transformation")

datadir: str = "data/crn_clean/"
resultdir: str = "data/amino_reaction_cycle_clean/"
graph_transformation.run(datadir, resultdir)

# ------------------------------------------------------------------------------

analysis_media_comparison = __import__("03_analysis_media_comparison")
analysis_media_comparison.run()


# ------------------------------------------------------------------------------
pathway_species = __import__("04_pathway_species")
datadir: str = "data/amino_reaction_cycle/"
resultdir: str = "data/pathway_species/"
pathway_species.run(datadir, resultdir)

# ------------------------------------------------------------------------------
acid_graph_enumeration = __import__("05_acid_graph_enumeration")
datadir: str = "data/crn/"
resultdir: str = "data/enumeration/"
acid_graph_enumeration.run(datadir, resultdir)

# ------------------------------------------------------------------------------
path_enumeration = __import__("05_path_enumeration")
infile = "data/enumeration/acacae_adam_CRN_glycine.pi"
path_enumeration.run(infile)
# ------------------------------------------------------------------------------
amino_acid_ratios = __import__("01_amino_acid_ratios")
"data/proteins/"
amino_acid_ratios.run()

# ------------------------------------------------------------------------------
pulp_solve = __import__("02_pulp_solve")
datadir: str = "data/amino_reaction_cycle/"
resultsdir: str = "data/flux_results/"
pulp_solve.run(datadir, resultdir)

# ------------------------------------------------------------------------------
flux_visualization = __import__("03_flux_visualization")
fluxDir: str = "data/flux_results/"
# flux_visualization.run(fluxDir)
# ------------------------------------------------------------------------------

atm_analysis = __import__("atm_analysis")
species_parent_dir: str = "data/atn_graphs/"
atm_analysis.run(species_parent_dir)
