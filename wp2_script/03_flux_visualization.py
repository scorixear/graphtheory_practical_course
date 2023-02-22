import networkx as nx
import pickle
import os
import matplotlib.pyplot as plt


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
    fluxDir: str = "data/flux_results/",
    aaCycleDir: str = "data/amino_reaction_cycle/",
) -> tuple([dict(), nx.DiGraph()]):

    flux = pickle.load(open(fluxDir + file, "rb"))

    name = file.replace("flux", "aa_cycle")
    graph = pickle.load(open(aaCycleDir + name, "rb"))

    return flux, graph


def run(fluxDir: str = "data/flux_results/", outdir: str = "./"):
    for file in os.listdir(fluxDir):
        flux, graph = read_file(file)

        # select a subgraph from all reactions that are not 0
        subgraphNodes = set()
        for node in flux:
            # filter out irrelevant reactions and export import nodes
            if flux[node] == 0 or not graph.has_node(node):
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

        reactionNodes = [
            node
            for node in aaSynthesisSubgraph
            if aaSynthesisSubgraph.nodes(data=True)[node]["nodeType"] == 1
        ]
        # nx.draw(aaSynthesisSubgraph, pos=nx.bipartite_layout(aaSynthesisSubgraph, reactionNodes),
        #        with_labels=True, node_size=20, font_size=4, width=0.3, arrowsize=5, node_color=colours, cmap=plt.cm.viridis)
        # plt.show()

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
        outfile = outdir + "graphics/" + file.replace("_flux.pi", "_aa_synthesis.png")
        plt.savefig(outfile, dpi=200)
        plt.close()


if __name__ == "__main__":
    run()
