import Bio.KEGG.REST as kegg

outdir = "data/kegg_pathways/"

def main():
    pathway_list = kegg.kegg_list("pathway")

    for pathway in pathway_list:
        pathway_id = pathway.split("\t")[0].replace("path:", "")
        pathway_name = pathway.split("\t")[1].strip()
        pathway_reactions = kegg.kegg_link("reaction", pathway_id)
        with open(f"{outdir}{pathway_id}.txt", "w") as outfile:
             print("#" + pathway_name, file=outfile)
             for reaction in pathway_reactions:
                try:
                    reaction_id = reaction.split("\t")[1].replace("rn:", "")
                except:
                 continue
                print(reaction_id, file=outfile, end="")

if __name__ == "__main__":
    main()