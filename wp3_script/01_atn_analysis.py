import os
import sys
import math
from typing import Callable
import networkx as nx
import all_cycles
import endpoint_bfs
import matplotlib.pyplot as plt
import numpy as np

sys.path.append("./library")

file_handler = __import__("file_handler")

class PlotData:
    def __init__(self) -> None:
        pass

class HistPlotData(PlotData):
    def __init__(self, values: list, bins: int, label: str) -> None:
        PlotData.__init__(self)
        self.values = values
        self.bins = bins
        self.label = label

class BarPlotData(PlotData):
    def __init__(self, data: list, labels: list[str], title: str, outputfile: str) -> None:
        PlotData.__init__(self)
        self.data = data
        self.labels = labels
        self.title = title
        self.outputfile = outputfile

def write_output(
        input_dir: str,
        output_dir: str,
        molecule_start_compound: str,
        molecule_start_atom: int,
        molecule_end_compound: str,
        endpoint_start_compound: str,
        endpoint_start_element: str,
        graph_generation: Callable[[str], nx.Graph]) -> list[PlotData]:

    endpoint_analysis: list[tuple[str, int]] = []
    molecule_paths: list[tuple[str, tuple[int, str]]] = []
    density_analysis: list[tuple[str, float]] = []
    component_numbers: list[tuple[str, int]] = []
    plots: list[PlotData] = []
    print("\n---- Starting Graph Reading ----")
    for entry in os.scandir(input_dir):
        if entry.is_file() and entry.name.endswith(".gml"):
            graph = graph_generation(input_dir+entry.name)
            output_file = output_dir+entry.name.split("/")[-1].split(".")[0]
            species_name = entry.name.split(".")[0]
            print(f"\nSpecies: {species_name}")
            print("0% Analysis values")
            component_number = get_number_of_compounds(graph)
            component_numbers.append((species_name, component_number))
            dens = density(graph)
            density_analysis.append((species_name, dens))
            print("25% Molecule Paths")
            molecule_path = find_paths(graph, molecule_start_compound, molecule_start_atom, molecule_end_compound)
            molecule_paths.append((species_name, (len(molecule_path[0].nodes), graph.nodes[molecule_path[1]]['element'])))
            print("50% Endpoints")
            endpoints = endpoint_bfs.bfs_endpoint(graph, endpoint_start_compound, endpoint_start_element)
            endpoint_analysis.append((species_name, len(endpoints)))
            print("75% Cycles")
            cycle_data = get_cycle_lengths(graph)
            results = {
                "conn_comp_number": component_number,
                "density": dens,
                "molecule_path": {
                    "start_compound": molecule_start_compound,
                    "start_atom_number": molecule_start_atom,
                    "start_atom_element": graph.nodes[molecule_path[1]]['element'],
                    "end_compound": molecule_end_compound,
                    "graph_size": len(molecule_path[0].nodes)
                },
                "endpoints": {
                    "amount": len(endpoints),
                    "data": endpoints
                },
                "connected_components": []
            }
            
            for _, value in cycle_data.items():
                component = {
                    "compound": value["compound"],
                    "element": value["element"],
                    "component_size": value['len'],
                    "cycle_basis": value["basis"],
                }
                if value["cycles_found"] > -1:
                    component["possible_cycles"] = value["possible"]
                    component["cycles_found"] = value["cycles_found"]
                    component["max_cycle_length"] = value["max_cycle"]
                    component["average_cycle_length"] = value["cycles"]
                    component["min_cycle_length"] = value["min_cycle"]
                results["connected_components"].append(component)
                        
            file_handler.write_json(results, output_file+".json")
            
            
            conn_comp_sizes = [cycle[1]['len'] for cycle in cycle_data.items()]

            plots.append(HistPlotData(
                conn_comp_sizes,
                (max(conn_comp_sizes)-min(conn_comp_sizes))//int(math.sqrt(len(conn_comp_sizes))),
                species_name
            ))
    
    plots.append(BarPlotData(
        [data[1] for data in endpoint_analysis], 
        [data[0] for data in endpoint_analysis],
        f"Endpoint Analysis: {endpoint_start_compound}_{endpoint_start_element}",
        "01_endpoint_analysis.png"))
    plots.append(BarPlotData(
        [data[1][0] for data in molecule_paths],
        [data[0] for data in molecule_paths],
        f"Molecule Paths: {molecule_start_compound}_{molecule_start_atom}_{molecule_paths[0][1][1]} -> {molecule_end_compound}",
        "02_molecule_paths.png"
    ))
    plots.append(BarPlotData(
        [data[1] for data in density_analysis],
        [data[0] for data in density_analysis],
        "Density Analysis",
        "03_density_analysis.png"
    ))
    plots.append(BarPlotData(
        [data[1] for data in component_numbers],
        [data[0] for data in component_numbers],
        "Connected Component Numbers",
        "04_component_numbers.png"
    ))
    return plots

def find_paths(graph: nx.Graph, molecule_a: str, identifier_a: int, molecule_b: str) -> tuple[nx.Graph, any]:
    new_graph = graph.copy()
    
    for u, v, edge_type in graph.edges(data="transition"):
        if edge_type == "TransitionType.NO_TRANSITION":
            new_graph.remove_edge(u,v)
    
    return bfs_to_molecule(new_graph, molecule_a, identifier_a, molecule_b)

def bfs_to_molecule(graph: nx.Graph, molecule_a: str, identifier_a: int, molecule_b: str) -> tuple[nx.Graph, any]:
    a_start_atom = None
    for key, data in graph.nodes(data=True):
        if data['compound_name'] == molecule_a and int(data['label'].split("_")[1]) == identifier_a:
            a_start_atom = key
            break
    a_graph = bfs_atn(graph, a_start_atom)
    a_element = graph.nodes[a_start_atom]["element"]
    b_graph = nx.Graph()
    for key, data in graph.nodes(data=True):
        if data["compound_name"] == molecule_b and data["element"] == a_element:
            atom_graph = bfs_atn(graph, key)
            b_graph = nx.compose(b_graph, atom_graph)
    return (nx.intersection(a_graph, b_graph), a_start_atom)
    

def bfs_atn(graph: nx.Graph, start_node: any):
    visited = set()
    queue = [start_node]
    
    while queue:
        next_node = queue.pop(0)
        if next_node in visited:
            continue
        visited.add(next_node)
        for neighbor in graph.neighbors(next_node):
            if neighbor not in visited and neighbor not in queue:
                queue.append(neighbor)
    return graph.subgraph(visited)

def filter_graph(filepath, filtered_items): 
    graph: nx.Graph = nx.read_gml(filepath, label=None)
    to_be_removed = []
    for node in graph.nodes.items():
        if node[1]["compound_name"] in filtered_items:
            to_be_removed.append(node[0])
    for node in to_be_removed:
        graph.remove_node(node)
    return graph

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
    return nx.density(t_graph)


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

def run(
    graphdir: str = "data/atn_graphs/",
    output_dir: str = "data/atn_analysis/",
    molecule_start_compound: str = "D-glucose",
    molecule_start_atom: int = 0,
    molecule_end_compound: str = "L-leucine",
    endpoint_start_compound: str = "D-glucose",
    endpoint_start_element: str = "C",
    inputdir: str = "input/"):
    data_full = write_output(
        graphdir,
        output_dir,
        molecule_start_compound,
        molecule_start_atom,
        molecule_end_compound,
        endpoint_start_compound,
        endpoint_start_element,
        lambda filepath: nx.read_gml(filepath, label=None))
    essential_compounds = file_handler.read_json(inputdir+"essential_compounds.json")
    data_cleaned = write_output(
        graphdir,
        output_dir,
        molecule_start_compound,
        molecule_start_atom,
        molecule_end_compound,
        endpoint_start_compound,
        endpoint_start_element,
        lambda filepath: filter_graph(filepath, essential_compounds))

    default_width = 0.5
    for i, full in enumerate(data_full):
        if isinstance(full, BarPlotData):
            full_bar: BarPlotData = full
            cleaned: BarPlotData = data_cleaned[i]
            _, ax = plt.subplots(figsize=(10,6))
            _ = ax.bar(
                np.arange(len(full_bar.labels))*2,
                full_bar.data,
                default_width,
                label="Full")
            _ = ax.bar(
                np.arange(len(cleaned.labels))*2+default_width,
                cleaned.data,
                default_width,
                label="No essential compounds"
            )
            
            ax.set_xlabel("Species")
            ax.set_title(full_bar.title)
            ax.set_xticks(np.arange(len(full_bar.labels))*2 + default_width/2)
            ax.set_xticklabels(full_bar.labels)
            ax.legend()
            
            plt.savefig(output_dir+"plots/"+full_bar.outputfile)
            plt.close()
        elif isinstance(full, HistPlotData):
            full_hist: HistPlotData = full
            cleaned: HistPlotData = data_cleaned[i]
            plt.yscale("log")
            plt.hist([full_hist.values, cleaned.values], bins=full_hist.bins)
            plt.legend(["Full", "No essential compounds"])
            plt.title(f"Connected Component Sizes: {full_hist.label}")
            plt.savefig(output_dir+"plots/"+full_hist.label+".png")
            plt.close()

if __name__ == "__main__":
    run()
