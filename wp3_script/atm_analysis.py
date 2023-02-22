import os
import networkx as nx
from typing import Callable
import all_cycles
import matplotlib.pyplot as plt
import math


def write_output(
    input_dir: str, output_dir: str, graph_generation: Callable[[str], nx.Graph]
) -> None:
    for entry in os.scandir(input_dir):
        if entry.is_file() and entry.name.endswith(".gml"):
            graph = graph_generation(input_dir + "/" + entry.name)
            output_file = output_dir + "/" + entry.name.split("/")[-1].split(".")[0]
            with open(output_file + ".txt", "w", encoding="UTF-8") as writer:
                writer.write(entry.name + "\n")
                writer.write(f"Component Number: {get_number_of_compounds(graph)}\n")
                writer.write(f"Density {density(graph)}\n\n")
                # component_sizes = '\n'.join(get_component_size(graph))
                # writer.write(f"Component Sizes: \n{component_sizes}\n\n")
                cycle_data = get_cycle_lengths(graph)
                for key_value in cycle_data.items():
                    writer.write(f"Compound: {key_value[1]['compound']}\n")
                    writer.write(f"Element: {key_value[1]['element']}\n")
                    writer.write(f"Compound Size: {key_value[1]['len']}\n")
                    writer.write(f"Cycle-Basis: {key_value[1]['basis']}\n")
                    if key_value[1]["cycles_found"] > -1:
                        writer.write(f"Possible Cycles: {key_value[1]['possible']}\n")
                        writer.write(f"Cycles found: {key_value[1]['cycles_found']}\n")
                        writer.write(f"Max Cycle Length: {key_value[1]['max_cycle']}\n")
                        writer.write(
                            f"Average Cycle Length: {key_value[1]['cycles']}\n"
                        )
                        writer.write(f"Min Cycle Length: {key_value[1]['min_cycle']}\n")
                    writer.write("\n")
            compound_sizes = [cycle[1]["len"] for cycle in cycle_data.items()]
            plt.yscale("log")
            plt.hist(
                compound_sizes,
                bins=(max(compound_sizes) - min(compound_sizes))
                // int(math.sqrt(len(compound_sizes))),
            )
            plt.savefig(output_file + ".png")
            plt.close()


def read_file(file: str) -> list[str]:
    """reads in file, splits by new line, trims output."""
    with open(file, "r", encoding="UTF-8") as reader:
        lines = reader.readlines()
    return [line.strip() for line in lines]


def filter_graph(filepath, filtered_items):

    graph: nx.Graph = nx.read_gml(filepath, label=None)
    to_be_removed = []
    for node in graph.nodes.items():
        if node[1]["compound_name"] in filtered_items:
            to_be_removed.append(node[0])
    for node in to_be_removed:
        graph.remove_node(node)
    return graph


def main():
    write_output(
        "data/atn_graphs",
        "data/atn_analysis",
        lambda filepath: nx.read_gml(filepath, label=None),
    )
    essential_compounds = read_file("wp1_script/essential_compounds.txt")
    write_output(
        "data/atn_graphs",
        "data/atn_analysis/no_essentials",
        lambda filepath: filter_graph(filepath, essential_compounds),
    )


# number of compounds
def get_number_of_compounds(t_graph: nx.Graph) -> int:
    new_graph = t_graph.copy()
    for u, v, edge_type in t_graph.edges(data="transition"):
        if edge_type == "TransitionType.NO_TRANSITION":
            new_graph.remove_edge(u, v)
    return nx.number_connected_components(new_graph)


def get_component_size(t_graph: nx.Graph) -> list[int]:
    new_graph = t_graph.copy()
    for u, v, edge_type in t_graph.edges(data="transition"):
        if edge_type == "TransitionType.NO_TRANSITION":
            new_graph.remove_edge(u, v)
    conn_comp = sorted(nx.connected_components(new_graph), key=len, reverse=True)
    # print(t_graph.nodes[list(conn_comp[0])[0]])
    return [
        f"{len(comp)}: {t_graph.nodes[list(comp)[0]]['element']} - {t_graph.nodes[list(comp)[0]]['compound_name']}"
        for comp in conn_comp
    ]


def density(t_graph: nx.DiGraph) -> float:
    return nx.number_of_edges(t_graph) / nx.number_of_nodes(t_graph)


def get_cycle_lengths(graph: nx.DiGraph) -> dict[int, dict[str, any]]:
    new_graph = graph.copy()
    cycles = dict()
    for u, v, edge_type in graph.edges(data="transition"):
        if edge_type == "TransitionType.NO_TRANSITION":
            new_graph.remove_edge(u, v)
    for i, conn_comp in enumerate(
        sorted(nx.connected_components(new_graph), key=len, reverse=True)
    ):
        basis = nx.cycle_basis(new_graph.subgraph(conn_comp))
        cycles[i] = {
            "len": len(conn_comp),
            "element": new_graph.nodes[list(conn_comp)[0]]["element"],
            "compound": new_graph.nodes[list(conn_comp)[0]]["compound_name"],
            "basis": len(basis),
            "possible": 2 ** (len(basis) - 1),
            "cycles_found": -1,
            "cycles": -1,
            "max_cycle": -1,
            "min_cycle": -1,
        }
        if 2 ** (len(basis) - 1) < 50000 and len(basis) > 0:
            cycles_found = [
                len(c)
                for c in all_cycles.find_all_cycles(new_graph.subgraph(conn_comp))
            ]  # nx.find_cycle(new_graph.subgraph(conn_comp))]
            cycles[i]["cycles_found"] = len(cycles_found)
            cycles[i]["cycles"] = sum(cycles_found) / len(cycles_found)
            cycles[i]["max_cycle"] = max(cycles_found)
            cycles[i]["min_cycle"] = min(cycles_found)
    return cycles


if __name__ == "__main__":
    run()
