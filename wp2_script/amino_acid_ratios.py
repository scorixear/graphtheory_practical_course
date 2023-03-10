def read_fasta(filepath: str) -> str:
    with open(filepath, "r", encoding="UTF-8") as fastaFile:
        proteomString = ""
        for line in fastaFile:
            if not line.startswith(">"):
                proteomString += line.strip()
        return proteomString


def calculate_ratios(proteomString: str, exceptions: list[str]) -> dict:
    """calculate the ratios for all aa that are not in the excetpion list"""

    translation_dir = {
        "L-cysteine": "C",
        "L-aspartate": "D",
        "L-serine": "S",
        "L-glutamine": "Q",
        "L-lysine": "K",
        "L-isoleucine": "I",
        "L-proline": "P",
        "L-threonine": "T",
        "L-phenylalanine": "F",
        "L-asparagine": "N",
        "glycine": "G",
        "L-histidine": "H",
        "L-leucine": "L",
        "L-arginine": "R",
        "L-tryptophan": "W",
        "L-alanine": "A",
        "L-valine": "V",
        "L-glutamate": "E",
        "L-tyrosine": "Y",
        "L-methionine": "M",
    }

    # remove unwanted aa (e.g. such that are not in the graph)
    for aa in exceptions:
        proteomString = proteomString.replace(translation_dir[aa], "")

    ratioDir = {}
    for [key, aa] in translation_dir.items():
        ratio = proteomString.count(aa) / len(proteomString)
        if ratio != 0:
            ratioDir[key] = ratio
    return ratioDir


def debug():
    datadir: str = "data/proteoms/"
    filepath = datadir + "proteom_ecoli_uniprot.fasta"
    proteomString = read_fasta(filepath)
    ratioDir = calculate_ratios(proteomString, ["L-cysteine"])
    print(ratioDir)


if __name__ == "__main__":
    debug()
