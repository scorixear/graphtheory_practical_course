import networkx as nx
import pickle


def graph_intersection_data(G: nx.DiGraph, H: nx.DiGraph) -> nx.DiGraph:
    R = G.copy()
    R.remove_nodes_from(n for n in G if n not in H)
    R.remove_edges_from(e for e in G.edges if e not in H.edges)
    return R


def graph_information(graph: nx.DiGraph, identifier: str):
    nt_ann = list(nx.get_node_attributes(graph, "nodeType").values())
    print(
        f"The {identifier} graph has {len(nt_ann)} nodes [{nt_ann.count(0)} compounds / {nt_ann.count(1)} reactions]"
    )


def compare_graphs(a_graph: nx.DiGraph, c_graph: nx.DiGraph):
    graph_information(a_graph, "adam")
    graph_information(c_graph, "cimIV")
    graph_information(graph_intersection_data(a_graph, c_graph), "intersection")

    nt_ann_1 = nx.get_node_attributes(a_graph, "nodeType")
    nt_ann_2 = nx.get_node_attributes(c_graph, "nodeType")
    nt_merged = nt_ann_1 | nt_ann_2

    unique_reaction_nodes = []
    unique_compound_nodes = []

    for key in nt_merged.keys():
        if (not (a_graph.has_node(key))) or (not (c_graph.has_node(key))):
            print(key, nt_merged[key])
            if nt_merged[key] == 1:
                unique_reaction_nodes.append(key)
            elif nt_merged[key] == 0:
                unique_compound_nodes.append(key)
            else:
                raise ValueError


def main():
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
        compare_graphs(agraph, cgraph)


if __name__ == "__main__":
    main()
