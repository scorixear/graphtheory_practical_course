import pulp
import os
import networkx as nx
import pickle
import amino_acid_ratios

INFINITY=1000

def add_acid_export_reactions(graph: nx.DiGraph, acid_ratios: dict):
    graph.add_node("R_output_BIOMASS", nodeType=1, reversible=False)
    for acid in acid_ratios:
        if acid in graph.nodes():
            # graph.add_node(f"R_reaction_{acid}", nodeType=1, reversible=False)
            graph.add_edge(acid, f"R_output_BIOMASS", multiplicity=acid_ratios[acid], direction=1)
            # graph.add_edge(f"R_reaction_{acid}", "BIOMASS", multiplicity=acids[acid], direction=1)
        else:
            print(f"Acid missing {acid}. Stopping LP")
            raise ValueError()

def read_file(file: str) -> list[str]:
    """reads in file, splits by new line, trims output."""
    input = open(file, "r")
    lines = input.readlines()
    return [line.strip() for line in lines]

def add_exchange_reactions(graph: nx.DiGraph, to_be_added: list[str]):
    nodes_to_be_added = []
    for node in graph.nodes(data=True):
        if node[1]['nodeType'] == 0:
            predecessors = list(graph.predecessors(node[0]))
            new_predecessors = [p for p in predecessors if graph.get_edge_data(p, node[0])['direction']==1]
            # if len(new_predecessors) == 0:
                # nodes_to_be_added.append(f"R_input_{node[0]}")
            successors = list(graph.successors(node[0]))
            new_successors = [s for s in successors if graph.get_edge_data(node[0], s)['direction']==1]
            if len(new_successors) == 0:
                nodes_to_be_added.append(f"R_output_{node[0]}")
    for node in to_be_added:
        if graph.has_node(node):
            nodes_to_be_added.append(f"R_input_{node}")

    for node in nodes_to_be_added:
        graph.add_node(node, nodeType = 1, reversible=False)
        if "input" in node:
            compound = node.split("_input_")[1]
            graph.add_edge(compound, node, multiplicity=1, direction=1)
        else:
            compound = node.split("_output_")[1]
            graph.add_edge(compound, node, multiplicity=1, direction=1)

def main():
    datadir = "data/amino_reaction_cycle/"
    proteom = amino_acid_ratios.read_fasta("data/proteins/proteom_ecoli_uniprot.fasta")
    amino_acids=read_file("wp1_script/amino_acids.txt")
    essential_compounds = read_file("wp1_script/essential_compounds.txt")
    for entry in os.scandir(datadir):
        if entry.is_file() and entry.name.endswith(".pi"):
            print(entry.path)
            graph: nx.DiGraph = pickle.load(open(entry.path, "rb"))
            model = pulp.LpProblem("Maximising_Problem", pulp.LpMaximize)
            
            missing_acids = [aa for aa in amino_acids if graph.has_node(aa) == False]
            print("Missing Acids:")
            print(missing_acids)
            acid_ratios = amino_acid_ratios.calculate_ratios(proteom, missing_acids)
            print(acid_ratios)
            # extend graph
            add_acid_export_reactions(graph, acid_ratios)
            add_exchange_reactions(graph, essential_compounds+ ['D-glucose'])

            variables = dict()
            for node in graph.nodes(data=True):
                if node[1]['nodeType'] == 1:
                    lowerBound = 0
                    if node[0].startswith ("R_input_"):
                        if graph.has_predecessor(node[0], "D-glucose"):
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


            #print(graph.get_edge_data("L-proline", "R_PROabc"))
            for node in graph.nodes(data=True):
                if node[1]['nodeType'] == 0:
                    constraints = []
                    for predeccessor in graph.predecessors(node[0]):
                        node_direction = graph.get_edge_data(predeccessor, node[0])['direction']
                        if node_direction == 2:
                            continue
                        coefficient = graph.get_edge_data(predeccessor, node[0])['multiplicity']
                        if coefficient == 0:
                            continue
                        constraints.append(coefficient*variables[predeccessor])
                    for successor in graph.successors(node[0]):
                        node_direction = graph.get_edge_data(node[0], successor)['direction']
                        if node_direction == 2:
                            continue
                        coefficient = -1 * graph.get_edge_data(node[0], successor)['multiplicity']#
                        if coefficient == 0:
                            continue
                        constraints.append(coefficient*variables[successor])
                    #if len(constraints) <=1:
                    #    continue
                    model += pulp.lpSum(constraints) == 0
            
            model.solve()
            pulp.LpStatus[model.status]
            input("ENTER for Variables...")
            relevant_counter = 0
            for v in variables:
                print(variables[v].name, ":", variables[v].varValue)
                if variables[v].varValue != 0:
                    relevant_counter+=1
            print(f"Relevant Reaction: {relevant_counter}")
            print(f"Biomass Output: {pulp.value(model.objective)}")
            input("ENTER for model...")
            print(model)
            break

if __name__ == "__main__":
    main()