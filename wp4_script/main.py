import sys
import os

sys.path.append("./wp1_script")
sys.path.append("./wp2_script")
sys.path.append("./wp3_script")
sys.path.append("./input")

def mkdir(path: str):
    try:
        os.makedirs(path)
    except OSError:
        pass


def main():

    data_directory = "data_clean/"
    input_directory = "input/"
    mkdir(data_directory)
    mkdir(input_directory)
    
    

    # parse all files, write graphs as pickle files
    
    smiles_to_crn = __import__("01_smiles_to_crn")
    smiles_list_directory: str = input_directory+"01_smiles_list/"
    mkdir(smiles_list_directory)
    rdir_crn: str = data_directory+"01_crn/"
    mkdir(rdir_crn)
    print("\n==== Smile to CRN ====")
    print("! Please provide all .smiles file in "+smiles_list_directory)
    print("! Please provide amino_acis.json under "+input_directory)
    print("! Please provide essential_compounds.json under "+input_directory)
    print("> Output will be under "+rdir_crn)
    input("< Press ENTER to continue...")
    smiles_to_crn.run(smiles_list_directory, rdir_crn, input_directory)
    print("> Parsing successfull")

    # ------------------------------------------------------------------------------
    # analysis scripts
    graph_transformation = __import__("02_graph_transformation")
    rdir_amino_reaction_cycle: str = data_directory+"02_amino_reaction_cycle/"
    mkdir(rdir_amino_reaction_cycle)
    print("\n==== CRN to Amino Acid Intersection ====")
    print("> Output will be under "+rdir_amino_reaction_cycle)
    input("< Press ENTER to continue...")
    graph_transformation.run(rdir_crn, rdir_amino_reaction_cycle)
    print("> Graph transformation successfull")

    # ------------------------------------------------------------------------------
    analysis_medium_comparison = __import__("03_analysis_media_comparison")
    rdir_medium_comparsion = data_directory+"03_medium_comparison/"
    mkdir(rdir_medium_comparsion)
    print("\n==== Medium Comparison ====")
    print("> Output will be under "+rdir_medium_comparsion)
    input("< Press ENTER to continue...")
    analysis_medium_comparison.run(rdir_amino_reaction_cycle, rdir_medium_comparsion)
    print("> Medium comparison successfull")

    # ------------------------------------------------------------------------------
    # compare species in same medium to each other
    
    pathway_species = __import__("04_pathway_species")
    rdir_species_comparison: str = data_directory+"04_pathway_species/"
    mkdir(rdir_species_comparison)
    mkdir(rdir_species_comparison+"adam")
    mkdir(rdir_species_comparison+"cimIV")
    print("\n====  Species Comparison ====")
    print("> Output will be under "+rdir_species_comparison)
    input("< Press ENTER to continue...")
    pathway_species.run(rdir_amino_reaction_cycle, rdir_species_comparison)
    print("> Species comparison successfull")

    # ------------------------------------------------------------------------------
    # create graphs from glucose forward bfs and single amino acids backward bfs
    acid_graph_enumeration = __import__("05_acid_graph_enumeration")
    rdir_acid_graphs: str = data_directory+"05_acid_graphs/"
    mkdir(rdir_acid_graphs)
    print("\n==== CRN to Single Acid Intersection")
    print("> Output will be under "+rdir_acid_graphs)
    input("< Press ENTER to continue...")
    acid_graph_enumeration.run(rdir_crn, rdir_acid_graphs)
    print("> Graph transformation successfull")

    # ------------------------------------------------------------------------------
    # enumerate paths (only console output, not save)
    path_enumeration = __import__("06_path_enumeration")
    infile = rdir_acid_graphs + "blongum_adam_glycine.pi"
    rdir_path_enumeration = data_directory+"06_path_enumeration/"
    mkdir(rdir_path_enumeration)
    print("\n==== Path Enumeration ====")
    print("> Enumeration paths of graph "+infile)
    print("! Please check if file exists")
    input("< Press ENTER to continue...")
    path_enumeration.run(infile, rdir_path_enumeration, input_directory)
    print("> Path enumeration successfull")
    
    # ------------------------------------------------------------------------------
    # generate amino acid synthesis visualitzations
    aa_synthesis_vis = __import__("07_aa_synthesis_vis")
    rdir_aa_synthesis_vis = data_directory+"07_aa_synthesis/"
    mkdir(rdir_aa_synthesis_vis)
    print("\n==== Amino acid synthesis visualization ====")
    print("> Output will be under "+rdir_aa_synthesis_vis)
    input("< Press ENTER to continue...")
    aa_synthesis_vis.run(rdir_crn, rdir_amino_reaction_cycle)
    print("> Amino acid synthesis visualization successfull")
    
    # ------------------------------------------------------------------------------
    # generate graph size visualization
    graph_size_analysis = __import__("08_graph_size_analysis")
    rdir_graph_size_analysis = data_directory+"08_graph_sizes"
    mkdir(rdir_graph_size_analysis)
    print("\n==== Graph size analysis ====")
    print("> Output will be under "+rdir_graph_size_analysis)
    input("< Press ENTER to continue...")
    graph_size_analysis.run(rdir_crn, rdir_amino_reaction_cycle, rdir_graph_size_analysis, input_directory)
    print("> Graph size analysis successfull")

    # ------------------------------------------------------------------------------
    # solve flux model
    pulp_solve = __import__("01_pulp_solve")
    rdir_pulp_solve: str = data_directory+"09_flux_results/"
    proteomdir: str = input_directory+"02_proteoms/"
    mkdir(rdir_pulp_solve)
    mkdir(proteomdir)
    print("\n==== Flux analysis ====")
    print("> Output will be under "+rdir_pulp_solve)
    print("! Please provide proteom files (.faa) under "+proteomdir)
    input("< Press ENTER to continue...")
    pulp_solve.run(rdir_amino_reaction_cycle, rdir_pulp_solve, proteomdir, input_directory)
    print("> Flux analysis successfull")

    # ------------------------------------------------------------------------------
    # Flux variation
    flux_variation = __import__("02_flux_variation")
    rdir_flux_variation = data_directory+"10_flux_variation"
    mkdir(rdir_flux_variation)
    print("\n==== Flux input constraint variation ====")
    print("> Output will be under "+rdir_flux_variation)
    input("< Press ENTER to continue...")
    flux_variation.run(rdir_amino_reaction_cycle, proteomdir, rdir_flux_variation, input_directory)
    print("> Flux variation successfull")
    
    # ------------------------------------------------------------------------------
    # Flux glucose variation
    flux_glucose_variation = __import__("03_flux_glucose_variation")
    rdir_flux_glucose_variation = data_directory+"11_flux_glucose_variation"
    mkdir(rdir_flux_glucose_variation)
    print("\n==== Flux glucose constraint variation ====")
    print("> Output will be under "+rdir_flux_glucose_variation)
    input("< Press ENTER to continue...")
    flux_glucose_variation.run(rdir_amino_reaction_cycle, proteomdir, rdir_flux_glucose_variation, input_directory)
    print("> Flux glucose variation successfull")

    # ------------------------------------------------------------------------------
    # Kegg pathways extraction
    extract_kegg_pathways = __import__("04_extract_kegg_pathways")
    rdir_extract_kegg_pathways = data_directory+"12_kegg_pathways"
    mkdir(rdir_extract_kegg_pathways)
    print("\n==== Kegg pathways extraction ====")
    print("> Output will be under "+rdir_extract_kegg_pathways)
    input("< Press ENTER to continue...")
    extract_kegg_pathways.run(rdir_extract_kegg_pathways)
    print("> Kegg pathways extraction successfull")

    # ------------------------------------------------------------------------------
    # Enrichment test
    enrichment_test = __import__("05_enrichment_test")
    rdir_enrichment_test= data_directory+"13_enrichment_test"
    mkdir(rdir_enrichment_test)
    print("\n==== Enrichment test ====")
    print("> Output will be under "+rdir_enrichment_test)
    print("! Please provide bigg_models_reactions.txt under "+input_directory)
    input("< Press ENTER to continue...")
    enrichment_test.run(rdir_pulp_solve, rdir_extract_kegg_pathways, rdir_enrichment_test, input_directory)
    print("> Enrichment test successfull")

    # ------------------------------------------------------------------------------
    # generate statistics about ATN
    atm_analysis = __import__("01_atn_analysis")
    graph_dir: str = input_directory+"03_atn_graphs/"
    rdir_atm_analysis: str = data_directory+"14_atn_analysis/"
    molecule_start_compound = "D-glucose"
    molecule_start_atom = 0
    molecule_end_compound = "L-leucine"
    endpoint_start_compound = "D-glucose"
    endpoint_start_element = "C"
    mkdir(graph_dir)
    mkdir(rdir_atm_analysis)
    mkdir(rdir_atm_analysis+"no_essentials")
    mkdir(rdir_atm_analysis+"plots")
    print("\n==== ATN Analysis ====")
    print("> Molecule Settings:")
    print("  > Molecule Start Compound: "+molecule_start_compound)
    print("  > Molecule Start Atom Number: "+molecule_start_atom)
    print("  > Molecule End Compound: "+molecule_end_compound)
    print("> Endpoint Settings:")
    print("  > Endpoint Start Compound: "+endpoint_start_compound)
    print("  > Endpoint Start Element: "+endpoint_start_element)
    print("> Output will be under "+rdir_atm_analysis)
    print("! Please provide ATN graphs (.gml) under "+graph_dir)
    input("< Press ENTER to continue...")
    atm_analysis.run(
        graph_dir, 
        rdir_atm_analysis,
        molecule_start_compound,
        molecule_start_atom,
        molecule_end_compound,
        endpoint_start_compound,
        endpoint_start_element,
        input_directory)
    print("> ATN analysis successfull")


if __name__ == "__main__":
    main()
