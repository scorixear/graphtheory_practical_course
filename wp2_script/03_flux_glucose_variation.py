"""
Do the flux analysis for different levels on glucose
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
flux_variation = __import__('02_flux_variation')

INFINITY = 1000
ESSENTIAL_COMPOUND_THRESHOLD = -100 # -200 is a good mixture where all compounds are present but not to abundant to get an infinite biomass function
GLUCOSE_THRESHOLDS = [-1, -2, -5, -10, -20, -50, -100, -200, -400, -600, -800, -1000]

def run(
    datadir: str = "data/amino_reaction_cycle_clean/",
    proteomdir: str = "data/proteoms/",
    outputdir: str = "data/flux_glucose_variations/",
    inputdir: str = "input/"):

    # all available acids
    amino_acids = file_handler.read_json(inputdir+"amino_acids.json")
    # every essential compound that can be used as an input
    essential_compounds = file_handler.read_json(inputdir+"essential_compounds.json")

    results = {"constrain": [], "biomass": {}, "activated_reactions": {}, "activation_level": {}}

    essential_compounds_constrains = flux_variation.create_compound_constrains(
            essential_compounds, ESSENTIAL_COMPOUND_THRESHOLD)
    

    for glucose_threshold in GLUCOSE_THRESHOLDS:
        essential_compounds_constrains["D-glucose"] = glucose_threshold
        results["constrain"].append(-1 * glucose_threshold)
        results["activation_level"][-1 * glucose_threshold] = {}
        
        for entry in os.scandir(datadir):
            if entry.is_file() and entry.name.endswith(".pi"):
                # load proteom of species
                organism_name = entry.name.split("_")[0]
                proteom = amino_acid_ratios.read_fasta(
                    f"{proteomdir}{organism_name}.faa")

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

                results["biomass"][organism_name].append(
                    pulp.value(model.objective))
                results["activated_reactions"][organism_name].append(
                    relevant_counter)
                results["activation_level"][-1*glucose_threshold][organism_name] = activation_level

    # plot the biomass
    legend = []
    for oragnism in results["biomass"]:
        plt.plot(results["constrain"], results["biomass"][oragnism])
        legend.append(oragnism)
    plt.xlabel("glucose threshold")
    plt.ylabel("biomass")
    plt.legend(legend)
    # plt.loglog()
    plt.savefig(f"{outputdir}biomass_glucose.png")
    plt.close()

    # plot the number of active reactions
    for organism in results["activated_reactions"]:
        plt.plot(results["constrain"], results["activated_reactions"][organism])
        plt.xlabel("glucose threshold")
    plt.ylabel("number of activated reactions")
    plt.legend(legend, loc= "lower right")
    # plt.loglog()
    plt.savefig(f"{outputdir}activated_reactions_glucose.png")
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
        plt.ylabel("Anzahl Reaktionen")
        plt.xlabel("Aktivierungslevel")
        plt.title(f"Aktivierungslevel bei glucose threshold {threshold}")
        plt.savefig(f"{outputdir}activation_level_by_glucose_threshold_{threshold}.png")
        plt.close()
    file_handler.write_json(results, outputdir+"settings.json")

if __name__ == "__main__":
    run()
