import sys
import os

sys.path.append("./wp1_script")
sys.path.append("./wp2_script")
sys.path.append("./wp3_script")

def mkdir(path: str):
    try:
        os.makedirs(path)
    except OSError:
        pass


def main():

    data_directory = "data_clean/"
    mkdir(data_directory)

    # step 1: parse all files, write graphs as pickle files
    smiles_to_crn = __import__("01_smiles_to_crn")
    smiles_list_directory: str = data_directory+"00_smiles_list_clean/"
    mkdir(smiles_list_directory)
    crn_save_directory: str = data_directory+"01_crn_clean/"
    mkdir(crn_save_directory)
    smiles_to_crn.run(smiles_list_directory, crn_save_directory)

    # ------------------------------------------------------------------------------
    # step 2: analysis scripts
    graph_transformation = __import__("02_graph_transformation")
    rdir_transformation: str = data_directory+"02_amino_reaction_cycle_clean/"
    mkdir(rdir_transformation)
    graph_transformation.run(crn_save_directory, rdir_transformation)

    # ------------------------------------------------------------------------------
    # not done since cleaned data do only have one media
    # analysis_media_comparison = __import__("03_analysis_media_comparison")
    # analysis_media_comparison.run(rdir_transformation)

    # ------------------------------------------------------------------------------
    # compare species in same medium to each other
    pathway_species = __import__("04_pathway_species")
    resultdir: str = data_directory+"04_pathway_species_clean/"
    mkdir(resultdir)
    mkdir(resultdir+"adam")
    mkdir(resultdir+"cimIV")
    pathway_species.run(rdir_transformation, resultdir)

    # ------------------------------------------------------------------------------
    # create graphs from glucose forward bfs and single amino acids backward bfs
    acid_graph_enumeration = __import__("05_acid_graph_enumeration")
    resultdir: str = data_directory+"05_enumeration_clean/"
    mkdir(resultdir)
    acid_graph_enumeration.run(crn_save_directory, resultdir)

    # ------------------------------------------------------------------------------
    # enumerate paths (only console output, not save)
    path_enumeration = __import__("05_path_enumeration")
    infile = resultdir + "blongum_adam_cleaned_CRN_glycine.pi"
    path_enumeration.run(infile)

    # ------------------------------------------------------------------------------
    # solve flux model
    pulp_solve = __import__("01_pulp_solve")
    resultsdir_psolve: str = data_directory+"06_flux_results_clean/"
    proteindir: str = data_directory+"06_proteins/"
    mkdir(resultsdir_psolve)
    mkdir(proteindir)
    pulp_solve.run(rdir_transformation, resultsdir_psolve, proteindir)

    # ------------------------------------------------------------------------------
    # visualize flux model
    flux_visualization = __import__("02_flux_visualization")
    graphics_save_directory: str = data_directory+"07_flux_visulation/"
    mkdir(graphics_save_directory)
    flux_visualization.run(resultsdir_psolve, rdir_transformation, graphics_save_directory)
    # ------------------------------------------------------------------------------
    # generate statistics about ATN
    atm_analysis = __import__("atm_analysis")
    graph_dir: str = data_directory+"08_atn_graphs/"
    analysis_dir: str = data_directory+"09_atn_analysis/"
    mkdir(graph_dir)
    mkdir(analysis_dir)
    mkdir(analysis_dir+"no_essentials")
    mkdir(analysis_dir+"plots")
    atm_analysis.run(graph_dir, analysis_dir)


if __name__ == "__main__":
    main()
