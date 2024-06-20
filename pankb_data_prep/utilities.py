import pandas as pd


def get_genome_list(args):
    with open(args.genomes, "r") as f:
        l = [genome.strip() for genome in f.readlines()]
    return l


COG_TABLE = pd.DataFrame.from_dict(
    {
        "A": "RNA processing and modification",
        "B": "Chromatin Structure and dynamics",
        "C": "Energy production and conversion",
        "D": "Cell cycle control and mitosis",
        "E": "Amino Acid metabolis and transport",
        "F": "Nucleotide metabolism and transport",
        "G": "Carbohydrate metabolism and transport",
        "H": "Coenzyme metabolis",
        "I": "Lipid metabolism",
        "J": "Tranlsation",
        "K": "Transcription",
        "L": "Replication and repair",
        "M": "Cell wall/membrane/envelop biogenesis",
        "N": "Cell motility",
        "O": "Post-translational modification, protein turnover, chaperone functions",
        "P": "Inorganic ion transport and metabolism",
        "Q": "Secondary Structure",
        "T": "Signal Transduction",
        "U": "Intracellular trafficing and secretion",
        "Y": "Nuclear structure",
        "Z": "Cytoskeleton",
        "R": "General Functional Prediction only",
        "S": "Function Unknown",
        "-": "Not found in COG",
    },
    orient="index",
    columns=["Function details"],
)
COG_TABLE.index.name = "Categories"
