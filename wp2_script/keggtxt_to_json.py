import os
import sys
sys.path.append("./library")
file_handler = __import__("file_handler")

def main():
    datadir = "data/kegg_pathways"
    output_file = "input/03_kegg_pathways/pathways.json"
    result = {}
    for entry in os.scandir(datadir):
        if entry.is_file() and entry.name.endswith(".txt"):
            with open(entry.path, "r", encoding="UTF-8") as reader:
                lines = reader.readlines()
            pathway_name = ""
            for line in lines:
                if line.startswith("#"):
                    pathway_name = line.split("#")[1].strip()
                    result[pathway_name] = []
                else:
                    result[pathway_name].append(line.strip())
    file_handler.write_json(result, output_file)
            
  
if __name__ == "__main__":
  main()