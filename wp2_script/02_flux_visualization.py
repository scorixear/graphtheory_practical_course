import os
import sys
import pickle
import networkx as nx
import matplotlib.pyplot as plt

sys.path.append("./library")
file_handler = __import__("file_handler")


def build_reaction_graph(G: nx.DiGraph) -> nx.DiGraph:
    """create a reaction graph consiting only of reactions"""
    reactionGraph = nx.DiGraph()
    reactionData = nx.get_node_attributes(G, "nodeType")
    for reactionNode in G:
        # only using reaction nodes
        if reactionData[reactionNode] == 1:
            for product in G.successors(reactionNode):

                if reactionData[product] == 1:
                    continue
                for followReaction in G.successors(product):
                    if reactionData[followReaction] == 1:
                        reactionGraph.add_edge(reactionNode, followReaction)
        else:
            continue

    return reactionGraph


def read_file(
    file: str,
    fluxDir: str,
    aaCycleDir: str,
) -> tuple[dict, nx.DiGraph]:

    flux = file_handler.read_json(fluxDir + file)

    name = file.replace("flux", "aa_cycle").replace(".json",".pi")
    with open(aaCycleDir+name, "rb") as reader:
        graph = pickle.load(reader)

    return flux, graph


def run(
    fluxDir: str = "data/flux_results/",
    aaCycleDir: str = "data/amino_reaction_cycle/",
    outdir: str = "data/flux_visualization",
):
    for file in os.listdir(fluxDir):
        flux, graph = read_file(file, fluxDir, aaCycleDir)

        # select a subgraph from all reactions that are not 0
        subgraphNodes = set()
        for [node, activation_value] in flux.items():
            # filter out irrelevant reactions and export import nodes
            if activation_value == 0 or not graph.has_node(node):
                continue

            subgraphNodes.add(node)
            # find the connected compounds
            for pre in graph.predecessors(node):
                subgraphNodes.add(pre)
            for succ in graph.successors(node):
                subgraphNodes.add(succ)

        # visualization of the whole graph (including reactions and compounds)
        aaSynthesisSubgraph = nx.subgraph(graph, subgraphNodes)
        colours = []
        for node in aaSynthesisSubgraph:
            if node in flux:
                colours.append(flux[node])
            else:
                colours.append(0)

        # visualize only the reactions
        reactionSubgraph = build_reaction_graph(aaSynthesisSubgraph)
        nx.draw(
            reactionSubgraph,
            pos=nx.shell_layout(reactionSubgraph),
            with_labels=True,
            node_size=20,
            font_size=4,
            width=0.3,
            arrowsize=5,
            node_color=[flux[node] for node in reactionSubgraph],
            cmap=plt.cm.viridis,
        )

        del reactionSubgraph
        del graph
        del aaSynthesisSubgraph
        outfile = outdir + file.replace("_flux.json", "_aa_synthesis.png")
        plt.savefig(outfile, dpi=200)
        plt.close()


if __name__ == "__main__":
    run()
