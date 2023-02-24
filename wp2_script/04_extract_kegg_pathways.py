import sys
import Bio.KEGG.REST as kegg

sys.path.append("./library")
file_handler = __import__("file_handler")


def run(outputfile: str = "data/kegg_pathways/pathways.json"):
    pathway_list = kegg.kegg_list("pathway")

    result = {}
    for pathway in pathway_list:
        pathway_id = pathway.split("\t")[0].replace("path:", "")
        pathway_name = pathway.split("\t")[1].strip()
        pathway_reactions = kegg.kegg_link("reaction", pathway_id)
        result[pathway_name] = []
        for reaction in pathway_reactions:
            try:
                result[pathway_name].append(reaction.split("\t")[1].replace("rn:", ""))
            except: # pylint: disable=W0702
                pass
    file_handler.write_json(result, outputfile)

if __name__ == "__main__":
    run()