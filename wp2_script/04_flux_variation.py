"""
Repeat the flux anaylsis over a list of constrains to visualize the biomass production at different constrain thresholds
 """

import pulp
import networkx as nx
import pickle
import os
import matplotlib.pyplot as plt

pulp_solve = __import__('02_pulp_solve')
amino_acid_ratios = __import__('01_amino_acid_ratios')

results_dir = "graphics/"
INFINITY = 1000
# define a list of constrains that should be analysed
constrains_list = [-1, -2, -5, -10, -20, -
                   50, -100, -200, -400, -600, -800, -1000]

def create_compound_constrains(essential_compounds: list[str], constraint: int) -> dict[str, int]:
    """Creates a dictionary with a constraint for all compounds except some that can be assumed as unlimited

    Args:
        essential_compounds list(str): list of essential compounds
        constraint (int): constraint for the essential compound

    Returns:
        dict[str, int]: constraint for each compound
    """
    # dictionary for the limitations of essential compounds (if compound listed in directory their lower of the import reaction is set to the value in the dic)
    essential_compounds_constrains = {
        es: constraint for es in essential_compounds}
    # set no limit for compounds that are definitly in the cell (and do not contain carbon)
    infinite_compounds = ["H2O", "phosphate", "H(+)"]
    for infinite_compound in infinite_compounds:
        essential_compounds_constrains[infinite_compound] = -INFINITY
    return essential_compounds_constrains


def main():

    # used for reading in graphs
    datadir = "data/amino_reaction_cycle_clean/"

    # proteom to be generated from acids
    # proteom = amino_acid_ratios.read_fasta("data/proteins/proteom_ecoli_uniprot.fasta")
    # all available acids
    amino_acids = pulp_solve.read_file("wp1_script/amino_acids.txt")
    # every essential compound that can be used as an input
    essential_compounds = pulp_solve.read_file(
        "wp1_script/essential_compounds.txt")

    results = {"constrain": [], "biomass": {},
               "activated_reactions": {}, "activation_level": {}}

    for constrain in constrains_list:
        results["constrain"].append(-1*constrain)
        results["activation_level"][-1*constrain] = dict()
        essential_compounds_constrains = create_compound_constrains(
            essential_compounds, constrain)
        for entry in os.scandir(datadir):
            if entry.is_file() and entry.name.endswith(".pi"):
                print(entry.path)

                # load proteom of species
                organism_name = entry.name.split("_")[0]
                proteom = amino_acid_ratios.read_fasta(
                    f"data/proteins/{organism_name}.faa")

                with open(entry.path, "rb") as reader:
                    graph: nx.DiGraph = pickle.load(reader)
                model = pulp.LpProblem("Maximising_Problem", pulp.LpMaximize)

                # find all missing acids in the graph
                missing_acids = [
                    aa for aa in amino_acids if graph.has_node(aa) is False]
                print("Missing Acids:")
                print(missing_acids)

                # calculate acid ratios used for the proteom
                acid_ratios = amino_acid_ratios.calculate_ratios(
                    proteom, missing_acids)

                # extend graph
                pulp_solve.add_acid_export_reactions(graph, acid_ratios)
                pulp_solve.add_exchange_reactions(
                    graph, essential_compounds + ['D-glucose'])

                variables: dict[str, pulp.LpVariable] = pulp_solve.generate_variables(
                    graph, essential_compounds_constrains)

                # objective function
                model += variables["R_output_BIOMASS"], "Profit"

                pulp_solve.add_constraints(model, graph, variables)

                model.solve()

                relevant_counter = 0
                # results: dict[str, any | int | None] = {}

                # collect activation level to visualize the proportion of highly activated reactions in a histogram
                activation_level = []
                for name, variable in variables.items():
                    #print(v, ":", variables[v].varValue)
                    # results[name] = variable.varValue
                    activation_level.append(abs(variable.varValue))
                    if variable.varValue != 0:
                        relevant_counter += 1

                # save biomass and number of activated reactions
                if organism_name not in results["biomass"]:
                    results["biomass"][organism_name] = []
                if organism_name not in results["activated_reactions"]:
                    results["activated_reactions"][organism_name] = []

                results["biomass"][organism_name].append(
                    pulp.value(model.objective))
                results["activated_reactions"][organism_name].append(
                    relevant_counter)
                results["activation_level"][-1 *
                                            constrain][organism_name] = activation_level

    # plot the biomass
    legend = []
    for oragnism in results["biomass"]:
        plt.plot(results["constrain"], results["biomass"][oragnism])
        legend.append(oragnism)
    plt.xlabel("compound threshold")
    plt.ylabel("biomass")
    plt.legend(legend)
    # plt.loglog()
    plt.savefig("graphics/biomass.png")
    plt.close()

    # plot the number of active reactions
    for organism in results["activated_reactions"]:
        plt.plot(results["constrain"], results["activated_reactions"][organism])
        plt.xlabel("compound threshold")
    plt.ylabel("number of activated reactions")
    plt.legend(legend, loc= "lower right")
    # plt.loglog()
    plt.savefig("graphics/activated_reactions.png")
    plt.close()

    # visualize activation level in a histogram
    for threshold in results["constrain"]:
        his_data = []
        legend = []
        for organism in results["biomass"]:
            his_data.append(results["activation_level"][threshold][organism])
            legend.append(organism)
        plt.hist(his_data, bins = 10, log=True)
        plt.legend(legend)
        plt.ylabel("number of reactions")
        plt.xlabel("reaction activation level")
        plt.title(f"activation level at essential compound threshold {threshold}")
        plt.savefig(f"{results_dir}activation_level_by_theshold_{threshold}.png")
        plt.close()

    # plot the proportion from unactivated to activated reactions
    for threshold in results["activation_level"]:
        activated = []
        not_activated = []
        legend = []
        for organism in results["activation_level"][threshold]:
            legend.append(organism)
            activation_levels = results["activation_level"][threshold][organism]
            active_percent = len([level for level in activation_levels if level != 0])/len(activation_levels) * 100
            activated.append(active_percent)
            not_activated.append(100-active_percent)
        
        plt.bar(legend, activated, width=0.5)
        plt.bar(legend, not_activated, width=0.5, bottom=activated)
        plt.legend(["active reactions", "inactive reactions"])
        plt.ylabel("percantage of total reactions")
        plt.savefig(f"{results_dir}acitive_reaction_proportions_compound_threshold_{threshold}.png")
        plt.close()

if __name__ == "__main__":
    main()
