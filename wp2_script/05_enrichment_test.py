import sys
import pickle
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
from scipy.stats import hypergeom

sys.path.append("./library")

file_handler = __import__("file_handler")

ORGANISM = "ecoli_adam" #limit the enrichment analysis to ecoli
ALPHA_NIVEAU = 0.05 # alpha niveau


def get_foreground_background(organism: str, flux_results: str) -> tuple[list[str], list[str]]:
    """Reads out the foreground (reactions with flux) and the background (all reactions)
    and returns a list of the keggIDs. Reactions that can not be mapped are discarded.

    Args:
        organism (str): organism name and medium

    Returns:
        foreground_kegg_ids (list[str]), background_kegg_ids (list[str])
    """
    # create a translation dict from bigg to kegg id
    bigg_to_kegg = dict()

    bigg_data = pd.read_csv("input/bigg_models_reactions.txt", sep="\t")
    bigg_data = bigg_data[bigg_data["database_links"].notna()]

    for _, row in bigg_data.iterrows():
        bigg_id = row["bigg_id"]
        db_ids = row["database_links"]

        if "KEGG" in db_ids:
            kegg_id = db_ids.split("KEGG Reaction: ")[
                1].split(";")[0].split("/")[-1]
            bigg_to_kegg[bigg_id] = kegg_id
    
    # load the reactions that have flux as foreground
    with open(f"{flux_results}{organism}_flux.pi", "rb") as reader:
        data = pickle.load(reader)
    foreground_kegg_ids = []
    background_kegg_ids = []
    for reaction in data:
        reaction_id = reaction.replace("R_", "")
        if reaction_id in bigg_to_kegg:
            if data[reaction] != 0:
                foreground_kegg_ids.append(bigg_to_kegg[reaction_id])
            background_kegg_ids.append(bigg_to_kegg[reaction_id])
    return foreground_kegg_ids, background_kegg_ids

def run(
    flux_resultdir: str = "data/flux_results/",
    kegg_file: str = "data/kegg_pathways/pathways.json",
    outputdir: str = "data/enrichment/"
):

    foreground, background = get_foreground_background(ORGANISM, flux_resultdir)
    results = []
    pathways = file_handler.read_json(kegg_file)
    for pathway_name, pathway_reactions in pathways.items():
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
        
        result = [pathway_name, p_value, k, N, n, M]
        results.append(result)

    #correct for multiple testing
    number_test = len(results)
    a_adj = ALPHA_NIVEAU/number_test
    results = [result for result in results if result[1] < a_adj]
    
    #calculate fold enrichment
    for result in results[:]:
        k = result[2]
        N = result[3]
        n = result[4]
        M = result[5]
        fold_enrichment = (k/N - n/M)/(n/M)
        result.append(fold_enrichment)
    
    #sort after p-value
    results = sorted(results, key=lambda result: result[6], reverse=True)
    
    pathway = [result[0].replace("Biosynthesis of ", "BS ").replace("biosynthesis", "BS") for result in results]
    p_value  = [result[1] for result in results]
    pathway[0] = "BS Phe, Tyr, Trp"
    reaction_number = [result[2] for result in results]
    fold_enrichment = [result[6] for result in results]
    y = range(len(pathway))

    plt.rcParams["figure.figsize"] = (15,5)
    scatter = plt.scatter(x=fold_enrichment, y=y, s=reaction_number, c=p_value, cmap=cm.copper)
    plt.yticks(ticks= y,labels= pathway,  fontsize= 8)
    plt.ylabel("pathway")
    plt.xlabel("Fold Enrichment")
    plt.legend(*scatter.legend_elements(), title = "p-value")
    plt.savefig(f"{outputdir}fold_enrichment.png")
    plt.close()


if __name__ == "__main__":
    run()
