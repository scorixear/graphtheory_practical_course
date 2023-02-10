import pulp
import os
import networkx as nx
import pickle

INFINITY=1000

def add_acid_export_reactions(graph: nx.DiGraph, acids: dict):
    graph.add_node("BIOMASS", nodeType=0)
    for acid in acids:
        if acid in graph.nodes():
            graph.add_node(f"R_reaction_{acid}", nodeType=1, reversible=False)
            graph.add_edge(acid, f"R_reaction_{acid}", multiplicity=1, direction=1)
            graph.add_edge(f"R_reaction_{acid}", "BIOMASS", multiplicity=acids[acid], direction=1)
        else:
            print(f"Acid missing {acid}. Stopping LP")
            raise ValueError()

def read_file(file: str) -> list[str]:
    """reads in file, splits by new line, trims output."""
    input = open(file, "r")
    lines = input.readlines()
    return [line.strip() for line in lines]

def add_exchange_reactions(graph: nx.DiGraph):
    nodes_to_be_added = []
    for node in graph.nodes(data=True):
        if node[1]['nodeType'] == 0:
            if len(list(graph.predecessors(node[0]))) == 0:
                nodes_to_be_added.append(f"R_input_{node[0]}")
            if len(list(graph.successors(node[0]))) == 0:
                nodes_to_be_added.append(f"R_output_{node[0]}")

    for node in nodes_to_be_added:
        graph.add_node(node, nodeType = 1, reversible=False)
        if "input" in node:
            compound = node.split("_input_")[1]
            graph.add_edge(node, compound, multiplicity=1, direction=1)
        else:
            compound = node.split("_output_")[1]
            graph.add_edge(compound, node, multiplicity=1, direction=1)

def main():
    datadir = "data/amino_reaction_cycle/"
    for entry in os.scandir(datadir):
        if entry.is_file() and entry.name.endswith("example_CRN.pi"):
            print(entry.path)
            graph: nx.DiGraph = pickle.load(open(entry.path, "rb"))
            model = pulp.LpProblem("Maximising_Problem", pulp.LpMaximize)
            
            amino_acids=read_file("wp1_script/amino_acids.txt")
            acids = dict()
            acids[7] = 0.33
            acids[9] = 0.33
            acids[10] = 0.33
            # extend graph
            add_acid_export_reactions(graph, acids)
            add_exchange_reactions(graph)

            variables = dict()
            for node in graph.nodes(data=True):
                if node[1]['nodeType'] == 1:
                    lowerBound = 0
                    if node[0].startswith ("R_input_"):
                        if graph.has_successor(node[0], 1):#if graph.has_successor(node[0], "D-glucose"):
                            lowerBound = -10
                        else:
                            lowerBound = -INFINITY
                    elif node[0].startswith("R_output_"):
                        lowerBound = 0
                    elif node[1]['reversible']:
                        lowerBound = -INFINITY
                    else:
                        lowerBound = 0
                    variables[node[0]]=pulp.LpVariable(node[0], lowBound=lowerBound, upBound=INFINITY, cat='Continuous')
            
            # objective function
            model += variables["R_output_BIOMASS"], "Profit"

            for node in graph.nodes(data=True):
                if node[1]['nodeType'] == 0:
                    constraints = []
                    for predeccessor in graph.predecessors(node[0]):
                        node_direction = graph.get_edge_data(predeccessor, node[0])['direction']
                        if node_direction == 2:
                            continue
                        coefficient = graph.get_edge_data(predeccessor, node[0])['multiplicity']
                        constraints.append(coefficient*variables[predeccessor])
                    for successor in graph.successors(node[0]):
                        node_direction = graph.get_edge_data(node[0], successor)['direction']
                        if node_direction == 2:
                            continue
                        coefficient = -1 * graph.get_edge_data(node[0], successor)['multiplicity']
                        constraints.append(coefficient*variables[successor])
                    model += pulp.lpSum(constraints) == 0
            
            model.solve()
            pulp.LpStatus[model.status]

            for v in variables:
                print(variables[v].name, ":", variables[v].varValue)
            print(pulp.value(model.objective))
            break

if __name__ == "__main__":
    main()