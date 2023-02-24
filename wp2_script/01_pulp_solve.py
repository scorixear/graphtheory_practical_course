import os
import sys
import pickle
import pulp
import networkx as nx

sys.path.append("./library")

file_handler = __import__("file_handler")
amino_acid_ratios = __import__("01_amino_acid_ratios")

INFINITY = 1000
DEFAULT_CONSTRAINT = -1


def add_acid_export_reactions(graph: nx.DiGraph, acid_ratios: dict):
    """Adds Reaction Node for Biomass called "R_output_BIOMASS" with the given multiplicities

    Args:
        graph (nx.DiGraph): the given graph to add the reaction node
        acid_ratios (dict): the acids and their rations needed to generate the biomass

    Raises:
        ValueError: Raised when an acid is not in the graph
    """
    graph.add_node("R_output_BIOMASS", nodeType=1, reversible=False)
    for acid in acid_ratios:
        if acid in graph.nodes():
            graph.add_edge(
                acid, 
                "R_output_BIOMASS",
                multiplicity=acid_ratios[acid], 
                direction=1)
        else:
            print(f"Acid missing {acid}. Stopping LP")
            raise ValueError()


def add_exchange_reactions(graph: nx.DiGraph, to_be_added: list[str]):
    """Adds R_output_[node] reaction nodes to the graph for every compound node that does not successors.
    Edges are specified as [node] -> R_output_[node]
    Also adds R_input_[node] reactions for every given node in to_be_added.
    Edges are specified as [node] -> R_input_[node]

    Args:
        graph (nx.DiGraph): the given graph to add nodes to
        to_be_added (list[str]): a list of additional nodes that require R_input_[node]'s
    """
    nodes_to_be_added = []
    for node in graph.nodes(data=True):
        # for every node that is a compound
        if node[1]["nodeType"] == 0:
            # collect all successors
            successors = list(graph.successors(node[0]))
            # filter successors to only be connected via an edge with direction = 1 label
            new_successors = [s 
                              for s in successors
                                if graph.get_edge_data(node[0], s)['direction'] == 1
                            ]
            # if no successor is given, add an output reaction node
            if len(new_successors) == 0:
                nodes_to_be_added.append(f"R_output_{node[0]}")
    # for every given node, add an input reaction node
    for node in to_be_added:
        if graph.has_node(node):
            nodes_to_be_added.append(f"R_input_{node}")
    # for every node that needs to be added
    for node in nodes_to_be_added:
        # add the node to the graph
        graph.add_node(node, nodeType=1, reversible=False)
        # and the required edges
        if "input" in node:
            compound = node.split("_input_")[1]
            graph.add_edge(compound, node, multiplicity=1, direction=1)
        else:
            compound = node.split("_output_")[1]
            graph.add_edge(compound, node, multiplicity=1, direction=1)

def generate_variables(graph: nx.DiGraph, constrains: dict[str, int]) -> dict[str, pulp.LpVariable]:
    """generates pulp.LpVariables for every reaction node

    Args:
        graph (nx.DiGraph): the graph containing the reaction nodes
        constrains (dict[str, int]): additional constrains for the essential compounds (and glucose)

    Returns:
        dict[str, pulp.LpVariable]: the dictionary containing reaction node names and their corresponding pulp variables
    """

    variables: dict[str, pulp.LpVariable] = {}
    # for every node that is a reaction node
    for node in graph.nodes(data=True):
        if node[1]["nodeType"] == 1:
            # figure out lower_bound
            lower_bound = 0
            # if input variable, set to -infinity -> infinity or to the specified constraint
            # except glucose, set to -10 -> infinity
            if node[0].startswith("R_input_"):
                compound = node[0].split("_input_")[1]
                
                if compound in constrains:
                    lower_bound = constrains[compound]
                elif graph.has_predecessor(node[0], "D-glucose"):
                    lower_bound = -10
                else:
                    lower_bound = -INFINITY
                
            # if output variable, set to 0 -> infinity
            elif node[0].startswith("R_output_"):
                lower_bound = 0
            # if reversible reaction, set to -infinity -> infinity
            # no output or input reaction is reversible
            elif node[1]["reversible"]:
                lower_bound = -INFINITY
            # create pulp variable
            variables[node[0]] = pulp.LpVariable(
                node[0], lowBound=lower_bound, upBound=INFINITY, cat="Continuous"
            )
    return variables


def add_constraints(
    lp_problem: pulp.LpProblem, graph: nx.DiGraph, variables: dict[str, pulp.LpVariable]
):
    """adds constraints to the lp_problem given the graph

    Args:
        lp_problem (pulp.LpProblem): the problem to add constraints to
        graph (nx.DiGraph): the graph to retrieve the data from
        variables (dict[str, pulp.LpVariable]): the variable to use for the nodes
    """
    # for every compound node in the graph
    for node in graph.nodes(data=True):
        if node[1]["nodeType"] == 0:
            # create constraint
            constraint = []
            # for every predecessor reaction node with direction=1 label
            for predeccessor in graph.predecessors(node[0]):
                node_direction = graph.get_edge_data(predeccessor, node[0])['direction']
                if node_direction == 2:
                    continue
                coefficient = graph.get_edge_data(predeccessor, node[0])['multiplicity']
                # ignore variables with 0 coefficient
                if coefficient == 0:
                    continue
                # append variable to constraint
                constraint.append(coefficient * variables[predeccessor])
            # for every successor reaction node with direction=1 label
            for successor in graph.successors(node[0]):
                node_direction = graph.get_edge_data(node[0], successor)["direction"]
                if node_direction == 2:
                    continue
                # set coefficient to -1 * multiplicity
                coefficient = (
                    -1 * graph.get_edge_data(node[0], successor)["multiplicity"]
                )
                # ignore variables with 0 coefficient
                if coefficient == 0:
                    continue
                constraint.append(coefficient * variables[successor])
            lp_problem += pulp.lpSum(constraint) == 0


def run(
    datadir: str = "data/amino_reaction_cycle/",
    resultsdir: str = "data/flux_results/",
    proteomdir: str = "data/proteoms/"
):
    # proteom to be generated from acids
    # proteom = amino_acid_ratios.read_fasta(proteindir + proteinfile)
    # all available acids
    amino_acids = file_handler.read_json("input/amino_acids.json")
    # every essential compound that can be used as an input
    essential_compounds = file_handler.read_json("input/essential_compounds.json")

    # dictionary for the limitations of essential compounds (if compound listed in directory their lower of the import reaction is set to the value in the dic)
    essential_compounds_constrains = {es : DEFAULT_CONSTRAINT for es in essential_compounds}
    # set no limit for compounds that are definitly in the cell (and do not contain carbon)
    relevant_compounds = ["H2O", "phosphate", "H(+)"]
    for relevant_compound in relevant_compounds:
        essential_compounds_constrains[relevant_compound] = -INFINITY


    for entry in os.scandir(datadir):
        if entry.is_file() and entry.name.endswith(".pi"):

            # load proteom of species
            organism_name = entry.name.split("_")[0]
            proteom = amino_acid_ratios.read_fasta(f"{proteomdir}{organism_name}.faa")

            with open(entry.path, "rb") as reader:
                graph: nx.DiGraph = pickle.load(reader)
            model = pulp.LpProblem("Maximising_Problem", pulp.LpMaximize)

            # find all missing acids in the graph
            missing_acids = [aa for aa in amino_acids if graph.has_node(aa) is False]

            # calculate acid ratios used for the proteom
            acid_ratios = amino_acid_ratios.calculate_ratios(proteom, missing_acids)

            # extend graph
            add_acid_export_reactions(graph, acid_ratios)
            add_exchange_reactions(graph, essential_compounds + ["D-glucose"])

            variables: dict[str, pulp.LpVariable] = generate_variables(graph, essential_compounds_constrains)

            # objective function
            model += variables["R_output_BIOMASS"], "Profit"

            add_constraints(model, graph, variables)

            model.solve()

            relevant_counter = 0
            results: dict[str, any | int | None] = {}

            for name, variable in variables.items():
                # print(v, ":", variables[v].varValue)
                results[name] = variable.varValue
                if variable.varValue != 0:
                    relevant_counter += 1
            # write out results
            file_handler.write_json(results, resultsdir + entry.name.replace("_aa_cycle", "_flux").replace(".pi",".json"))


if __name__ == "__main__":
    run()
