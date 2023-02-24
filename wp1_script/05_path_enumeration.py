import sys
import pickle
import networkx as nx

sys.path.append("./input")
file_handler = __import__("file_handler")


def read_file(file: str) -> list[str]:
    """reads in file, splits by new line, trims output."""
    input_str = open(file, "r", encoding="UTF-8")
    lines = input_str.readlines()
    return [line.strip() for line in lines]


def build_reaction_graph(
    G: nx.DiGraph, source: str, acids_present: list[str], common_compounds: list[str]
) -> nx.DiGraph:
    """create a reaction graph consiting only of reactions, glucose and amino acids.
    We are taking all products (except the common ones) and create edges from them.
    """
    reactionGraph = nx.DiGraph()
    reactionData = nx.get_node_attributes(G, "nodeType")
    for reactionNode in G:
        # only using reaction nodes
        if reactionData[reactionNode] == 1:
            for product in G.successors(reactionNode):
                # we only want to follow our carbon therefore exclude links via H, ATP ect.
                if product in common_compounds or reactionData[product] == 1:
                    continue
                for followReaction in G.successors(product):
                    if reactionData[followReaction] == 1:
                        reactionGraph.add_edge(reactionNode, followReaction)
        else:
            continue
    # add source and target( glucose and aa)
    for reaction in G.successors(source):
        reactionGraph.add_edge(source, reaction)

    for a in acids_present:
        for reaction in G.predecessors(a):
            reactionGraph.add_edge(reaction, a)

    return reactionGraph


def enumerate_simple_paths(
    reactionGraph: nx.DiGraph(), source: str, target: str, maxPaths: int
) -> list[str]:
    reverseG = nx.reverse(reactionGraph)
    # pathsIter = nx.all_simple_paths(reactionGraph, source, target, cutoff=50)
    pathsIter = nx.all_simple_paths(reverseG, target, source, cutoff=100)
    paths = []

    for path in pathsIter:
        if len(paths) <= maxPaths:
            paths.append(list(reversed(path)))
        else:
            break
    return paths


def enumerate_simple_paths_deletion(
    reactionGraph: nx.DiGraph(), source: str, target: str, maxPaths: int
) -> list[str]:
    reverseG = nx.reverse(reactionGraph)
    # pathsIter = nx.all_simple_paths(reactionGraph, source, target, cutoff=50)
    pathsIter = nx.all_simple_paths(reverseG, target, source, cutoff=100)
    paths = []

    for path in pathsIter:
        if len(paths) <= maxPaths:
            paths.append(list(reversed(path)))
        else:
            break
    return paths


def run(
    infile: str = "data/enumeration/acacae_adam_CRN_glycine.pi",
    output_folder: str = "data/paths/"):
    source = "D-glucose"
    with open(infile, "rb") as reader:
        graph = pickle.load(reader)
    
    amino_acids = file_handler.read_json("input/amino_acids.json")
    essential_compounds = file_handler.read_json("input/essential_compounds.json")
    
    acids_present = [aa for aa in amino_acids if graph.has_node(aa)]

    reactionGraph = build_reaction_graph(graph, source, acids_present, essential_compounds)

    aaSimplePaths = {}
    pathThreshold = 1000
    results = {}
    for aa in acids_present:
        paths = enumerate_simple_paths(reactionGraph, source, aa, pathThreshold)
        aaSimplePaths[aa] = paths
        
    for [acid_key, paths] in aaSimplePaths.items():

        results[acid_key] = {}
        pathLenghts = [len(path) for path in paths]

        # extract all paths of equal size
        pathSizes = set(pathLenghts)
        paths_equal_length = {}
        for pathSize in pathSizes:
            paths_equal_length[pathSize] = []
        for path in paths:
            paths_equal_length[len(path)].append(path)

        for [length, equal_paths] in paths_equal_length.items():
            pathSet = [set(path) for path in equal_paths]
            sharedPath = set.intersection(*pathSet)
            results[acid_key][length] = sharedPath
    file_handler.write_json(results, output_folder+"acid_paths.json")

if __name__ == "__main__":
    run()
