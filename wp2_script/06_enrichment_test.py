import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import pickle
import networkx as nx
import os
from scipy.stats import hypergeom
import matplotlib.cm as cm

bigg_dir = "data/bigg/"
kegg_dir = "data/kegg_pathways/"
organism = "ecoli_adam" #limit the enrichment analysis to ecoli
a = 0.05 # alpha niveau


def read_kegg(file:str) -> tuple([str, list[str]]):
    """Reads the reactions for a pathway

    Args:
        file (str): pathway file from Kegg

    Returns:
        str: name of the pathway, list of Kegg reaction IDs of the pathway
    """
    pathway_reactions = []
    pathway_name = ""
    with open(file, "r") as infile:
        for line in infile:
            if line.startswith("#"):
                pathway_name = line.replace("#", "").strip()
            else:
                pathway_reactions.append(line.strip())
    
    return pathway_name, pathway_reactions


def get_foreground_background(organism: str) -> tuple([list[str], list[str]]):
    """Reads out the foreground (reactions with flux) and the background (all reactions)
    and returns a list of the keggIDs. Reactions that can not be mapped are discarded.

    Args:
        organism (str): organism name and medium

    Returns:
        foreground_kegg_ids (list[str]), background_kegg_ids (list[str])
    """
    # create a translation dict from bigg to kegg id
    bigg_to_kegg = dict()

    bigg_data = pd.read_csv(f"{bigg_dir}bigg_models_reactions.txt", sep="\t")
    bigg_data = bigg_data[bigg_data["database_links"].notna()]

    for index, row in bigg_data.iterrows():
        bigg_id = row["bigg_id"]
        db_ids = row["database_links"]

        if "KEGG" in db_ids:
            kegg_id = db_ids.split("KEGG Reaction: ")[
                1].split(";")[0].split("/")[-1]
            bigg_to_kegg[bigg_id] = kegg_id
    
    # load the reactions that have flux as foreground
    data = pickle.load(
        open(f"data/flux_results_clean/{organism}_cleaned_flux.pi", "rb"))
    foreground_kegg_ids = []
    background_kegg_ids = []
    for reaction in data:
        reaction_id = reaction.replace("R_", "")
        if reaction_id in bigg_to_kegg:
            if data[reaction] != 0:
                foreground_kegg_ids.append(bigg_to_kegg[reaction_id])
            background_kegg_ids.append(bigg_to_kegg[reaction_id])
    return foreground_kegg_ids, background_kegg_ids

def main():

    foreground, background = get_foreground_background(organism)
    results = []

    for pathway_file in os.listdir(kegg_dir):
        pathway_name, pathway_reactions = read_kegg(kegg_dir + pathway_file)
        
        #do not test to short pathways
        if len(pathway_reactions) < 30:
            continue
            
        #calculate the variables for the hypergeometric test
        #background 
        M = len(background)
        #number of elements in background associated with the pathway
        n = len(set(pathway_reactions).intersection(set(background)))
        #foreground
        N = len(foreground)
        #significant reactions in foreground
        k = len(set(pathway_reactions).intersection(set(foreground)))

        p_value = hypergeom.sf(k-1,M,n,N)
        
        # if p_value < 0.05:
        #     print(f"found significant pathway: {pathway_name}, p-value: {p_value}, proportion of associated reactions in background: {round(n/M, 3)}, proportion of associated reactions in foreground: {round(k/N, 3)}")
        
        result = [pathway_name, p_value, k, N, n, M]
        results.append(result)

    #correct for multiple testing
    number_test = len(results)
    a_adj = a/number_test
    results = [result for result in results if result[1] < a_adj]
    
    #calculate fold enrichment
    for result in results[:]:
        k = result[2]
        N = result[3]
        n = result[4]
        M = result[5]
        fold_enrichment = (k/N - n/M)/(n/M)
        result.append(fold_enrichment)

    # results = [tuple(result) for result in results]
    
    #sort after p-value
    results = sorted(results, key=lambda result: result[6], reverse=True)
    
    pathway = [result[0].replace("Biosynthesis of ", "BS ").replace("biosynthesis", "BS") for result in results]
    p_value  = [result[1] for result in results]
    pathway[0] = "BS Phe, Tyr, Trp"
    reaction_number = [result[2] for result in results]
    fold_enrichment = [result[6] for result in results]
    y = range(len(pathway))
    
    # plt.style.use('_mpl-gallery')
    plt.rcParams["figure.figsize"] = (15,5)
    # matplotlib.rcParams['figure.subplot.left'] = 0
    # matplotlib.rcParams['figure.subplot.bottom'] = 0
    # matplotlib.rcParams['figure.subplot.right'] = 1
    # matplotlib.rcParams['figure.subplot.top'] = 1
    scatter = plt.scatter(x=fold_enrichment, y=y, s=reaction_number, c=p_value, cmap=cm.copper)
    plt.yticks(ticks= y,labels= pathway,  fontsize= 8)
    plt.ylabel("pathway")
    plt.xlabel("Fold Enrichment")
    plt.legend(*scatter.legend_elements(), title = "p-value")
    plt.savefig("graphics/fold_enrichment.png")
    plt.close()


if __name__ == "__main__":
    main()
