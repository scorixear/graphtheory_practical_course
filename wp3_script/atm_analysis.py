import networkx as nx
import matplotlib.pyplot as plt
import os 
import all_cycles

def main():
    species_parent_dir = "data/atn_graphs/"
    graph = nx.read_gml(species_parent_dir+"blongum_adam_cleaned.gml", label=None)
    for entry in os.scandir(species_parent_dir):
        if entry.is_file() and entry.name.endswith(".gml"):
            print(entry.name)
            graph = nx.read_gml(species_parent_dir+entry.name, label=None)
            print(f"Component Number: {get_number_of_compounds(graph)}")
            print(f"Component Sizes: {get_component_size(graph)}")
            print(f"Density {density(graph)}")
            cycle_basiseseses = get_cycle_lengths(graph)
            for cycle_basis in cycle_basiseseses:
                if len(cycle_basiseseses[cycle_basis]) != 0:
                    lcb = len(cycle_basiseseses[cycle_basis])
                    print(f"Cycle-Basis: {cycle_basis}: {lcb}")
                    print(f"Possible Cycles: {2**lcb}")
            #print(f"Cycle Length {len(cycles)}")
            print("\n")

# number of compounds
def get_number_of_compounds(t_graph: nx.Graph) -> int:
    new_graph = t_graph.copy()
    for u, v, edge_type in t_graph.edges(data="transition"):
        if edge_type == "TransitionType.NO_TRANSITION":
            new_graph.remove_edge(u,v)
    return nx.number_connected_components(new_graph)

def get_component_size(t_graph: nx.Graph) -> list[int]:
    new_graph = t_graph.copy()
    componentSizes = []
    for u, v, edge_type in t_graph.edges(data="transition"):
        if edge_type == "TransitionType.NO_TRANSITION":
            new_graph.remove_edge(u,v)
    for conn_comp in nx.connected_components(new_graph):
        componentSizes.append(len(conn_comp))
    return componentSizes

def density(t_graph: nx.DiGraph) -> float:
    return nx.number_of_edges(t_graph)/nx.number_of_nodes(t_graph)

def get_cycle_lengths(graph: nx.DiGraph) -> dict[int, list]:
    new_graph = graph.copy()
    cycles = dict()
    for u, v, edge_type in graph.edges(data="transition"):
        if edge_type == "TransitionType.NO_TRANSITION":
            new_graph.remove_edge(u,v)
    for i, conn_comp in enumerate(nx.connected_components(new_graph)):
        cycles[i] = nx.cycle_basis(new_graph.subgraph(conn_comp))#all_cycles.find_all_cycles(new_graph.subgraph(conn_comp))
        #print(f"{i}: {len(cycles[i])}")
    return cycles
    
if __name__ == "__main__":
    main()