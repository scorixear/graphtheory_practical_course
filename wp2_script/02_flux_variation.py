"""
Repeat the flux anaylsis over a list of constrains 
to visualize the biomass production at different constrain thresholds
 """

import os
import sys
import pickle
import pulp
import networkx as nx
import matplotlib.pyplot as plt
import amino_acid_ratios

sys.path.append("./library")
file_handler = __import__("file_handler")

pulp_solve = __import__('01_pulp_solve')

INFINITY = 1000
# define a list of constrains that should be analysed
DEFAULT_CONSTRAINTS = [-1, -2, -5, -10, -20, -
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


def run(
    datadir: str = "data/amino_reaction_cycle_clean/",
    proteomdir: str = "data/proteoms/",
    outputdir: str = "data/flux_variation/",
    inputdir: str = "input/"):

    # proteom to be generated from acids
    # proteom = amino_acid_ratios.read_fasta("data/proteins/proteom_ecoli_uniprot.fasta")
    # all available acids
    amino_acids = file_handler.read_json(inputdir+"amino_acids.json")
    # every essential compound that can be used as an input
    essential_compounds = file_handler.read_json(inputdir+"essential_compounds.json")

    results = {"constrain": [], "biomass": {},
               "activated_reactions": {}, "activation_level": {}}

    for constrain in DEFAULT_CONSTRAINTS:
        results["constrain"].append(-1*constrain)
        results["activation_level"][-1*constrain] = {}
        essential_compounds_constrains = create_compound_constrains(essential_compounds, constrain)
        for entry in os.scandir(datadir):
            if entry.is_file() and entry.name.endswith(".pi"):

                # load proteom of species
                organism_name = entry.name.split("_")[0]
                proteom = amino_acid_ratios.read_fasta(f"{proteomdir}{organism_name}.faa")

                with open(entry.path, "rb") as reader:
                    graph: nx.DiGraph = pickle.load(reader)
                model = pulp.LpProblem("Maximising_Problem", pulp.LpMaximize)

                # find all missing acids in the graph
                missing_acids = [
                    aa for aa in amino_acids if graph.has_node(aa) is False]

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

                # collect activation level to visualize the proportion of highly activated reactions in a histogram
                activation_level = []
                for _, variable in variables.items():
                    activation_level.append(abs(variable.varValue))
                    if variable.varValue != 0:
                        relevant_counter += 1

                # save biomass and number of activated reactions
                if organism_name not in results["biomass"]:
                    results["biomass"][organism_name] = []
                if organism_name not in results["activated_reactions"]:
                    results["activated_reactions"][organism_name] = []

                results["biomass"][organism_name].append(pulp.value(model.objective))
                results["activated_reactions"][organism_name].append(relevant_counter)
                results["activation_level"][-1 * constrain][organism_name] = activation_level

    # plot the biomass
    legend = []
    for oragnism in results["biomass"]:
        plt.plot(results["constrain"], results["biomass"][oragnism])
        legend.append(oragnism)
    plt.xlabel("compound threshold")
    plt.ylabel("biomass")
    plt.legend(legend)
    plt.savefig(f"{outputdir}biomass.png")
    plt.close()

    # plot the number of active reactions
    for organism in results["activated_reactions"]:
        plt.plot(results["constrain"], results["activated_reactions"][organism])
        plt.xlabel("compound threshold")
    plt.ylabel("number of activated reactions")
    plt.legend(legend, loc= "lower right")
    # plt.loglog()
    plt.savefig(f"{outputdir}activated_reactions.png")
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
        plt.savefig(f"{outputdir}activation_level_by_theshold_{threshold}.png")
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
        plt.savefig(f"{outputdir}acitive_reaction_proportions_compound_threshold_{threshold}.png")
        plt.close()
        
    file_handler.write_json(results, outputdir+"settings.json")

if __name__ == "__main__":
    run()
