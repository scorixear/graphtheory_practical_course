import os
import math
import networkx as nx
from typing import Callable
import all_cycles
import endpoint_bfs
import matplotlib.pyplot as plt
import numpy as np

def main():
    data_full = write_output(
        "data/atn_graphs",
        "data/atn_analysis",
        "D-glucose",
        0,
        "L-leucine",
        "D-glucose",
        "C",
        lambda filepath: nx.read_gml(filepath, label=None))
    essential_compounds = read_file("wp1_script/essential_compounds.txt")
    data_cleaned = write_output(
        "data/atn_graphs",
        "data/atn_analysis/no_essentials",
        "D-glucose",
        0,
        "L-leucine",
        "D-glucose",
        "C",
        lambda filepath: filter_graph(filepath, essential_compounds))
    default_width = 0.5
    for i, full in enumerate(data_full):
        cleaned = data_cleaned[i]
        _, ax = plt.subplots(figsize=(10,6))
        _ = ax.bar(
            np.arange(len(full[0]))*2,
            full[1],
            default_width,
            label="Full")
        _ = ax.bar(
            np.arange(len(cleaned[0]))*2+default_width,
            cleaned[1],
            default_width,
            label="No essential compounds"
        )
        
        ax.set_xlabel("Species")
        ax.set_title(full[2])
        ax.set_xticks(np.arange(len(full[0]))*2 + default_width/2)
        ax.set_xticklabels(full[0])
        ax.legend()
        
        plt.savefig("data/atn_analysis/"+full[3])
        plt.close()
        
def write_output(
        input_dir: str,
        output_dir: str,
        molecule_start_compound: str,
        molecule_start_atom: int,
        molecule_end_compound: str,
        endpoint_start_compound: str,
        endpoint_start_element: str,
        graph_generation: Callable[[str], nx.Graph]) -> list[tuple[list[str], list, str, str]]:
    endpoint_analysis: list[tuple[str, int]] = []
    molecule_paths: list[tuple[str, tuple[int, str]]] = []
    density_analysis: list[tuple[str, float]] = []
    component_numbers: list[tuple[str, int]] = []
    print("Starting Graph Reading")
    for entry in os.scandir(input_dir):
        if entry.is_file() and entry.name.endswith(".gml"):
            graph = graph_generation(input_dir+"/"+entry.name)
            output_file = output_dir+"/"+entry.name.split("/")[-1].split(".")[0]
            species_name = entry.name.split(".")[0]
            
            print(f"Generating analysis values for species {species_name}")
            component_number = get_number_of_compounds(graph)
            component_numbers.append((species_name, component_number))
            dens = density(graph)
            density_analysis.append((species_name, dens))
            print(f"Calculating Molecule Paths for species {species_name}")
            molecule_path = find_paths(graph, molecule_start_compound, molecule_start_atom, molecule_end_compound)
            molecule_paths.append((species_name, (len(molecule_path[0].nodes), graph.nodes[molecule_path[1]]['element'])))
            print(f"Calculating Endpoints for species {species_name}")
            endpoints = endpoint_bfs.bfs_endpoint(graph, endpoint_start_compound, endpoint_start_element)
            #endpoint_str = "\n".join([f"{key}: {len(data)}" for [key, data] in endpoints.items()])
            endpoint_analysis.append((species_name, len(endpoints)))
            print(f"Calculating Cycles for species {species_name}")
            cycle_data = get_cycle_lengths(graph)
            
            with open(output_file+".txt", "w", encoding="UTF-8") as writer:
                writer.write(species_name+"\n")    
                writer.write(f"Component Number: {component_number}\n")
                writer.write(f"Density {dens}\n")
                writer.write(f"Molecule Paths:\n{molecule_start_compound}_{molecule_start_atom}_{graph.nodes[molecule_path[1]]['element']} -> {molecule_end_compound} = {len(molecule_path[0].nodes)}\n\n")
                writer.write(f"Endpoint analysis: {len(endpoints)}\n\n")#\n{endpoint_str}\n\n")
                for key_value in cycle_data.items():
                    writer.write(f"Compound: {key_value[1]['compound']}\n")
                    writer.write(f"Element: {key_value[1]['element']}\n")
                    writer.write(f"Compound Size: {key_value[1]['len']}\n")
                    writer.write(f"Cycle-Basis: {key_value[1]['basis']}\n")
                    if key_value[1]['cycles_found'] > -1:
                        writer.write(f"Possible Cycles: {key_value[1]['possible']}\n")
                        writer.write(f"Cycles found: {key_value[1]['cycles_found']}\n")
                        writer.write(f"Max Cycle Length: {key_value[1]['max_cycle']}\n")
                        writer.write(f"Average Cycle Length: {key_value[1]['cycles']}\n")
                        writer.write(f"Min Cycle Length: {key_value[1]['min_cycle']}\n")
                    writer.write("\n")
            compound_sizes = [cycle[1]['len'] for cycle in cycle_data.items()]
            plt.yscale('log')
            plt.hist(
                compound_sizes, 
                bins=(max(compound_sizes)-min(compound_sizes))//int(math.sqrt(len(compound_sizes))),
                label="Cycle Analysis: "+species_name)
            plt.savefig(output_file+".png")
            plt.close()
    
    #plt.xticks(rotation=90)
    return [
        (
            [data[0] for data in endpoint_analysis], 
            [data[1] for data in endpoint_analysis], 
            f"Endpoint Analysis: {endpoint_start_compound}_{endpoint_start_element}",
            "01_endpoint_analysis.png"
        ),
        (
            [data[0] for data in molecule_paths],
            [data[1][0] for data in molecule_paths],
            f"Molecule Paths: {molecule_start_compound}_{molecule_start_atom}_{molecule_paths[0][1][1]} -> {molecule_end_compound}",
            "02_molecule_paths.png"
        ),
        (
            [data[0] for data in density_analysis],
            [data[1] for data in density_analysis],
            "Density Analysis",
            "03_density_analysis.png"
        ),
        (
            [data[0] for data in component_numbers],
            [data[1] for data in component_numbers],
            "Component Numbers",
            "04_component_numbers.png"
        )]
def read_file(file: str) -> list[str]:
    """reads in file, splits by new line, trims output."""
    with open(file, "r", encoding="UTF-8") as reader:
        lines = reader.readlines()
    return [line.strip() for line in lines]

def find_paths(graph: nx.Graph, molecule_a: str, identifier_a: int, molecule_b: str) -> tuple[nx.Graph, any]:
    new_graph = graph.copy()
    
    for u, v, edge_type in graph.edges(data="transition"):
        if edge_type == "TransitionType.NO_TRANSITION":
            new_graph.remove_edge(u,v)
    
    return bfs_to_molecule(new_graph, molecule_a, identifier_a, molecule_b)

def bfs_to_molecule(graph: nx.Graph, molecule_a: str, identifier_a: int, molecule_b: str) -> tuple[nx.Graph, any]:
    atom_a = None
    molecule_b_nodes = []
    final_molecule_b_nodes = []

    for [key, data] in graph.nodes(data=True):
        if data['compound_name'] == molecule_a and int(data['label'].split("_")[1]) == identifier_a:
            atom_a = key
        elif data['compound_name'] == molecule_b:
            molecule_b_nodes.append(key)
    for atom_b in molecule_b_nodes:
        if graph.nodes[atom_b]['element'] == graph.nodes[atom_a]['element']:
            final_molecule_b_nodes.append(atom_b)
    
    visited = set()
    queue = [atom_a]

    while queue:
        next_node = queue.pop(0)
        if next_node in final_molecule_b_nodes or next_node in visited:
            continue
        visited.add(next_node)
        for neighbor in graph.neighbors(next_node):
            if neighbor not in visited and neighbor not in queue:
                queue.append(neighbor)
    return (graph.subgraph(visited), atom_a)

def filter_graph(filepath, filtered_items): 
    graph: nx.Graph = nx.read_gml(filepath, label=None)
    to_be_removed = []
    for node in graph.nodes.items():
        if node[1]['compound_name'] in filtered_items:
            to_be_removed.append(node[0])
    for node in to_be_removed:
        graph.remove_node(node)
    return graph

# number of compounds
def get_number_of_compounds(t_graph: nx.Graph) -> int:
    new_graph = t_graph.copy()
    for u, v, edge_type in t_graph.edges(data="transition"):
        if edge_type == "TransitionType.NO_TRANSITION":
            new_graph.remove_edge(u,v)
    return nx.number_connected_components(new_graph)

def get_component_size(t_graph: nx.Graph) -> list[int]:
    new_graph = t_graph.copy()
    for u, v, edge_type in t_graph.edges(data="transition"):
        if edge_type == "TransitionType.NO_TRANSITION":
            new_graph.remove_edge(u,v)
    conn_comp = sorted(nx.connected_components(new_graph), key=len, reverse=True)
    #print(t_graph.nodes[list(conn_comp[0])[0]])
    return [f"{len(comp)}: {t_graph.nodes[list(comp)[0]]['element']} - {t_graph.nodes[list(comp)[0]]['compound_name']}" for comp in conn_comp]

def density(t_graph: nx.DiGraph) -> float:
    return nx.number_of_edges(t_graph)/nx.number_of_nodes(t_graph)

def get_cycle_lengths(graph: nx.DiGraph) -> dict[int, dict[str, any]]:
    new_graph = graph.copy()
    cycles = dict()
    for u, v, edge_type in graph.edges(data="transition"):
        if edge_type == "TransitionType.NO_TRANSITION":
            new_graph.remove_edge(u,v)
    for i, conn_comp in enumerate(sorted(nx.connected_components(new_graph), key=len, reverse=True)):
        basis = nx.cycle_basis(new_graph.subgraph(conn_comp))
        cycles[i] = {
                     'len': len(conn_comp),
                     'element': new_graph.nodes[list(conn_comp)[0]]['element'],
                     'compound': new_graph.nodes[list(conn_comp)[0]]['compound_name'],
                     'basis': len(basis),
                     'possible': 2**(len(basis)-1),
                     'cycles_found': -1,
                     'cycles': -1,
                     'max_cycle': -1,
                     'min_cycle': -1,
                    }
        if 2**(len(basis)-1) < 50000 and len(basis) > 0:
            cycles_found = [len(c) for c in all_cycles.find_all_cycles(new_graph.subgraph(conn_comp))]# nx.find_cycle(new_graph.subgraph(conn_comp))]
            cycles[i]['cycles_found'] = len(cycles_found)
            cycles[i]['cycles'] = sum(cycles_found)/len(cycles_found)
            cycles[i]['max_cycle'] = max(cycles_found)
            cycles[i]['min_cycle'] = min(cycles_found)
    return cycles
    
if __name__ == "__main__":
    main()