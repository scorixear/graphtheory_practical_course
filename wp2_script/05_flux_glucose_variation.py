"""
Do the flux analysis for different levels on glucose
"""

import pulp
import networkx as nx
import pickle
import os
import matplotlib.pyplot as plt

pulp_solve = __import__('02_pulp_solve')
amino_acid_ratios = __import__('01_amino_acid_ratios')
flux_variation = __import__('04_flux_variation')

results_dir = "graphics/"
INFINITY = 1000
essential_compound_threshold = -100 #-200 is a good mixture where all compounds are present but not to abundant to get an infinite biomass function
glucose_thresholds = [-1, -2, -5, -10, -20, -50, -100, -200, -400, -600, -800, -1000]

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

    results = {"constrain": [], "biomass": {}, "activated_reactions": {}, "activation_level": {}}

    essential_compounds_constrains = flux_variation.create_compound_constrains(
            essential_compounds, essential_compound_threshold)
    

    for glucose_threshold in glucose_thresholds:
        essential_compounds_constrains["D-glucose"] = glucose_threshold
        results["constrain"].append(-1*glucose_threshold)
        results["activation_level"][-1*glucose_threshold] = dict()
        
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

                # collect activation level to visualize the proportion of highly activated reactions in a histogram
                activation_level = []
                for name, variable in variables.items():
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
    plt.savefig(f"{results_dir}biomass_glucose.png")
    plt.close()

    # plot the number of active reactions
    for organism in results["activated_reactions"]:
        plt.plot(results["constrain"], results["activated_reactions"][organism])
        plt.xlabel("glucose threshold")
    plt.ylabel("number of activated reactions")
    plt.legend(legend, loc= "lower right")
    # plt.loglog()
    plt.savefig(f"{results_dir}activated_reactions_glucose.png")
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
        plt.savefig(f"{results_dir}activation_level_by_glucose_threshold_{threshold}.png")
        plt.close()

if __name__ == "__main__":
    main()
