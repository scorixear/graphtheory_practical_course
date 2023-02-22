import sys
import networkx as nx
import os
import pickle

sys.path.append("./wp1_script")
sys.path.append("./wp2_script")
sys.path.append("./wp3_script")


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

# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------

analysis_media_comparison = __import__("03_analysis_media_comparison")

data_parent_directory = "data/amino_reaction_cycle/"
organisms = [
    "acacae",
    "blongum",
    "bproducta",
    "btheta",
    "cbuty",
    "ecoli",
    "eramosum",
    "lplantarum",
]

for organism in organisms:
    print(f"-------------- {organism} --------------")
    a_fpath = data_parent_directory + organism + "_adam_aa_cycle.pi"
    c_fpath = data_parent_directory + organism + "_cimIV_aa_cycle.pi"
    with open(a_fpath, "rb") as a_reader:
        agraph = pickle.load(a_reader)
    with open(c_fpath, "rb") as c_reader:
        cgraph = pickle.load(c_reader)
    analysis_media_comparison.compare_graphs(agraph, cgraph)


# ------------------------------------------------------------------------------


pathway_species = __import__("04_pathway_species")

# data folder holding the subgraph pickle files
datadir = "data/amino_reaction_cycle/"
adam_species = []
cimIV_species = []
# for each subgraph
for entry in os.scandir(datadir):
    if entry.is_file() and entry.name.endswith(".pi"):
        # read in graph and sort them in both mediums
        with open(entry.path, "rb") as reader:
            graph_object = {"graph": pickle.load(reader), "name": entry.name}
            if "adam" in entry.name:
                adam_species.append(graph_object)
            else:
                cimIV_species.append(graph_object)
result_str = ""
# for each medium
# compare species to every other species
for index in range(len(adam_species)):
    for other_index in range(index + 1, len(adam_species)):
        result_str = (
            result_str
            + "\n"
            + pathway_species.compare(
                adam_species[index], adam_species[other_index], "adam"
            )
        )
for index in range(len(cimIV_species)):
    for other_index in range(index + 1, len(cimIV_species)):
        result_str = (
            result_str
            + "\n"
            + pathway_species.compare(
                cimIV_species[index], cimIV_species[other_index], "cimIV"
            )
        )
with open("data/pathway_species/results.txt", "w", encoding="UTF-8") as writer:
    writer.write(result_str)

# ------------------------------------------------------------------------------
datadir = "data/crn/"
resultdir = "data/enumeration/"
for entry in os.scandir(datadir):
    if entry.is_file() and entry.name.endswith(".pi"):
        print(entry.path)
        with open(entry.path, "rb") as graph_file:
            graph: nx.DiGraph = pickle.load(graph_file)
        essential_compounds = metabolite_subgraph.read_file(
            "wp1_script/essential_compounds.txt"
        )
        ec_clean = [ec for ec in essential_compounds if graph.has_node(ec)]
        essential_compounds = ec_clean
        glucose_graph: nx.DiGraph = metabolite_subgraph.bf_search(
            graph, "D-glucose", essential_compounds
        )

        enumeration_graph = glucose_graph.copy()
        enumeration_graph.remove_nodes_from(n for n in essential_compounds)

        amino_acids = metabolite_subgraph.read_file("wp1_script/amino_acids.txt")
        for acid in amino_acids:
            try:
                reversed_graph = metabolite_subgraph.reverse_bf_search(
                    enumeration_graph, acid
                )
                output_graph = nx.DiGraph()
                output_graph = nx.compose(reversed_graph, output_graph)
                print(f"{acid}: {len(reversed_graph.nodes())}")
                # try:
                #    print(nx.shortest_path(output_graph, "D-glucose", acid))
                # except:
                #    print(f"No path found between D-glucose to {acid}")
                file_path = f"{resultdir}/{entry.name.split('.')[0]}_{acid}.pi"
                with open(file_path, "wb") as writer:
                    pickle.dump(output_graph, writer)
            except nx.NetworkXError:
                print(f"Acid missing: {acid}")

# ------------------------------------------------------------------------------
path_enumeration = __import__("05_path_enumeration")
infile = "data/enumeration/acacae_adam_CRN_glycine.pi"
source = "D-glucose"
G = pickle.load(open(infile, "rb"))
amino_acids = path_enumeration.read_file("wp1_script/amino_acids.txt")
common_compounds = path_enumeration.read_file("wp1_script/essential_compounds.txt")
acids_present = [aa for aa in amino_acids if G.has_node(aa)]

reactionGraph = path_enumeration.build_reaction_graph(
    G, source, acids_present, common_compounds
)
# TODO filter subgraphs for each aa

aaSimplePaths = dict()
pathThreshold = 1000
for aa in acids_present:
    print("enumerating ", aa)
    paths = path_enumeration.enumerate_simple_paths(
        reactionGraph, source, aa, pathThreshold
    )
    aaSimplePaths[aa] = paths
    # try:
    #     print("shortest graph: ", nx.shortest_path(reactionGraph, source, aa))
    # except nx.NetworkXNoPath:
    #     print("no shortest path found for: ", aa)
    #     pass
for aa in aaSimplePaths:
    print(f"path number for {aa}: {len(aaSimplePaths[aa])}")

    pathLenghts = [len(path) for path in aaSimplePaths[aa]]
    # plt.hist(pathLenghts)
    # plt.show()

    # extract all paths of equal size
    pathSizes = set(pathLenghts)
    paths_equal_length = dict()
    for pathSize in pathSizes:
        paths_equal_length[pathSize] = []
    for path in aaSimplePaths[aa]:
        paths_equal_length[len(path)].append(path)

    for pathLength in paths_equal_length:
        pathSet = [set(path) for path in paths_equal_length[pathLength]]
        sharedPath = set.intersection(*pathSet)
        print(
            f"für AS {aa} haben die Pfade der Länge {pathLength} {len(sharedPath)} gemeinsame Reaktionen"
        )
        print(f"gemeinsam sind die folgenden Reaktionen: {sharedPath}")
        print()

# ------------------------------------------------------------------------------
amino_acid_ratios = __import__("01_amino_acid_ratios")
dataDir = "data/proteins/"
filepath = dataDir + "proteom_ecoli_uniprot.fasta"
proteomString = amino_acid_ratios.read_fasta(filepath)
ratioDir = amino_acid_ratios.calculate_ratios(proteomString, ["L-cysteine"])
print(ratioDir)

# ------------------------------------------------------------------------------
pulp_solve = __import__("02_pulp_solve")
pulp_solve.main()

# ------------------------------------------------------------------------------
flux_visualization = __import__("03_flux_visualization")
# flux_visualization.main()
# ------------------------------------------------------------------------------
atm_analysis = __import__("atm_analysis")
atm_analysis.main()
