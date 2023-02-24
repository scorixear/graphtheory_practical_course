import sys
import os.path as path
import pickle
import networkx as nx

sys.path.append("./library")
file_handler = __import__("file_handler")

def graph_intersection_data(G: nx.DiGraph, H: nx.DiGraph) -> nx.DiGraph:
    R = G.copy()
    R.remove_nodes_from(n for n in G if n not in H)
    R.remove_edges_from(e for e in G.edges if e not in H.edges)
    return R


def graph_information(graph: nx.DiGraph):
    nt_ann = list(nx.get_node_attributes(graph, "nodeType").values())
    return {
        "nodes": len(nt_ann),
        "compounds": nt_ann.count(0),
        "reactions": nt_ann.count(1)
    }

def compare_graphs(a_graph: nx.DiGraph, c_graph: nx.DiGraph):
    result_data = {}
    result_data["adam"] = graph_information(a_graph)
    result_data["cimIV"] = graph_information(c_graph)
    result_data["intersection"] = graph_information(graph_intersection_data(a_graph, c_graph))

    nt_ann_1 = nx.get_node_attributes(a_graph, "nodeType")
    nt_ann_2 = nx.get_node_attributes(c_graph, "nodeType")
    nt_merged = nt_ann_1 | nt_ann_2

    unique_adam_reaction_nodes = []
    unique_adam_compound_nodes = []
    unique_cimIV_reaction_nodes = []
    unique_cimIV_compound_nodes = []

    for [key, value] in nt_merged.items():
        if not a_graph.has_node(key):
            if value == 1:
                unique_adam_reaction_nodes.append(key)
            elif value == 0:
                unique_adam_compound_nodes.append(key)
            else:
                raise ValueError
        elif not c_graph.has_node(key):
            if value == 1:
                unique_cimIV_reaction_nodes.append(key)
            elif value == 0:
                unique_cimIV_compound_nodes.append(key)
            else:
                raise ValueError
    
    result_data["unique"] = {
        "reaction": {
            "adam": unique_adam_reaction_nodes,
            "cimIV": unique_cimIV_reaction_nodes
            },
        "compound": {
            "adam": unique_adam_compound_nodes,
            "cimIV": unique_cimIV_compound_nodes
        }}
    return result_data


def run(
    datadir: str = "data/amino_reaction_cycle/", 
    result_dir: str = "data/media_comparison/"):
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
        a_fpath = datadir + organism + "_adam_aa_cycle.pi"
        c_fpath = datadir + organism + "_cimIV_aa_cycle.pi"
        if path.exists(a_fpath) and path.exists(c_fpath):
            with open(a_fpath, "rb") as a_reader:
                agraph = pickle.load(a_reader)
            with open(c_fpath, "rb") as c_reader:
                cgraph = pickle.load(c_reader)
            comparison = compare_graphs(agraph, cgraph)
            output_file = result_dir+organism+".json"
            file_handler.write_json(comparison, output_file)


if __name__ == "__main__":
    run()
