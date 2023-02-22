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

rdir_transformation: str = "data/amino_reaction_cycle_clean/"
graph_transformation.run(crn_save_directory, rdir_transformation)

# ------------------------------------------------------------------------------

# analysis_media_comparison = __import__("03_analysis_media_comparison")
# analysis_media_comparison.run(rdir_transformation)

# ------------------------------------------------------------------------------
pathway_species = __import__("04_pathway_species")
resultdir: str = "data/pathway_species_clean/"
pathway_species.run(rdir_transformation, resultdir)

# ------------------------------------------------------------------------------
acid_graph_enumeration = __import__("05_acid_graph_enumeration")
resultdir: str = "data/enumeration_clean/"
acid_graph_enumeration.run(crn_save_directory, resultdir)

# ------------------------------------------------------------------------------
path_enumeration = __import__("05_path_enumeration")
infile = resultdir + "blongum_adam_cleaned_CRN_glycine.pi"
path_enumeration.run(infile)

# ------------------------------------------------------------------------------
pulp_solve = __import__("02_pulp_solve")
datadir: str = "data/amino_reaction_cycle_clean/"
resultsdir_psolve: str = "data/flux_results_clean/"
pulp_solve.run(datadir, resultsdir_psolve)

# ------------------------------------------------------------------------------
flux_visualization = __import__("03_flux_visualization")
graphics_save_directory: str = "graphics/"
aaCycleDir: str = "data/amino_reaction_cycle_clean/"
flux_visualization.run(resultsdir_psolve, aaCycleDir, graphics_save_directory)
# ------------------------------------------------------------------------------

atm_analysis = __import__("atm_analysis")
graph_dir: str = "data/atn_graphs/"
analysis_dir: str = "data/atn_analysis/"
atm_analysis.run(graph_dir, analysis_dir)
