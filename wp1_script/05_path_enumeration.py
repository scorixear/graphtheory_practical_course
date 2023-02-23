import networkx as nx
import pickle


def read_file(file: str) -> list[str]:
    """reads in file, splits by new line, trims output."""
    input = open(file, "r")
    lines = input.readlines()
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


def run(infile: str = "data/enumeration/acacae_adam_CRN_glycine.pi"):
    source = "D-glucose"
    G = pickle.load(open(infile, "rb"))
    amino_acids = read_file("wp1_script/amino_acids.txt")
    common_compounds = read_file("wp1_script/essential_compounds.txt")
    acids_present = [aa for aa in amino_acids if G.has_node(aa)]

    reactionGraph = build_reaction_graph(G, source, acids_present, common_compounds)
    # TODO filter subgraphs for each aa

    aaSimplePaths = dict()
    pathThreshold = 1000
    for aa in acids_present:
        print("enumerating ", aa)
        paths = enumerate_simple_paths(reactionGraph, source, aa, pathThreshold)
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


if __name__ == "__main__":
    run()
