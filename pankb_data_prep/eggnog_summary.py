import logging
from pathlib import Path
import math
import pandas as pd
from Bio import SeqIO

def get_cog_dict():
    """
    Get COG dict
    """

    cog_dict = {
        "A": "RNA processing and modification",
        "B": "Chromatin structure and dynamics",
        "C": "Energy production and conversion",
        "D": "Cell cycle control, cell division, chromosome partitioning",
        "E": "Amino acid transport and metabolism",
        "F": "Nucleotide transport and metabolism",
        "G": "Carbohydrate transport and metabolism",
        "H": "Coenzyme transport and metabolism",
        "I": "Lipid transport and metabolism",
        "J": "Translation, ribosomal structure and biogenesis",
        "K": "Transcription",
        "L": "Replication, recombination and repair",
        "M": "Cell wall/membrane/envelope biogenesis",
        "N": "Cell motility",
        "O": "Post-translational modification, protein turnover, and chaperones",
        "P": "Inorganic ion transport and metabolism",
        "Q": "Secondary metabolites biosynthesis, transport, and catabolism",
        "R": "General function prediction only",
        "S": "Function unknown",
        "T": "Signal transduction mechanisms",
        "U": "Intracellular trafficking, secretion, and vesicular transport",
        "V": "Defense mechanisms",
        "W": "Extracellular structures",
        "X": "Mobilome: prophages, transposons",
        "Y": "Nuclear structure",
        "Z": "Cytoskeleton",
        "-": "Not found in COG"
    }

    return cog_dict


def generate_eggnog_summary(gp_binary_path, summary_path, eggnog_table_path, reference_path, eggnog_summary_path):
    df_roary_binary = pd.read_csv(gp_binary_path, index_col = 0, low_memory = False)
    df_pangene_summary = pd.read_csv(summary_path, index_col=0)
    # eggnog_table_path = '/home/binhuan/data_pankb/bgcflow/data/interim/eggnog_roary/' + species + '/' + species + '.emapper.annotations'

    df_eggnog = pd.read_csv(eggnog_table_path, sep="\t", header=4, index_col="#query").iloc[:-3,:]
    df_eggnog.index.name = "locus_tag"

    # Sort the matrix by the sum of strains presence
    df = pd.DataFrame()
    df["data"] = df_roary_binary.sum(axis=1).sort_values(ascending=True)

    # Group the data by value and count the number of occurrences
    df = df.groupby('data').size().reset_index(name='frequency')
    x15 = math.floor(df['data'].quantile(0.15))
    x99 = math.floor(df['data'].quantile(0.99).round()) - 1

    df_accesory = df_roary_binary[(df_roary_binary.sum(1) >= x15) & (df_roary_binary.sum(1) < x99)]

    cog_dict = get_cog_dict()

    # fasta_file = '/home/binhuan/data_pankb/bgcflow/data/interim/roary/' + species + '/pan_genome_reference.fa'
    recs = SeqIO.parse(reference_path, format="fasta")
    for seq in recs:
        desc = seq.description
        locus_tag = desc.split(" ")[0]
        group_id = desc.split(" ")[1]
        if locus_tag in df_eggnog.index:
            df_eggnog.loc[locus_tag, "group_id"] = group_id

    for pan_gene_id in df_pangene_summary.index:
        if pan_gene_id in df_eggnog.group_id.tolist():
            locus_tag = df_eggnog[df_eggnog.group_id == pan_gene_id].index.tolist()[0]
            for col in df_eggnog.columns:
                df_pangene_summary.loc[pan_gene_id, col] = df_eggnog.loc[locus_tag, col]

            cog_cat = df_eggnog.loc[locus_tag, 'COG_category']
            if len(cog_cat) > 1:
                cog_name = ' | '.join([cog_dict[cog_letter] for cog_letter in cog_cat])
                df_pangene_summary.loc[pan_gene_id, 'COG_category_name'] = cog_name
                df_pangene_summary.loc[pan_gene_id, 'uniq_COG_category_name'] = 'multi_COGs_' + str(len(cog_cat))
            else:
                df_pangene_summary.loc[pan_gene_id, 'COG_category_name'] = cog_dict[cog_cat]
                df_pangene_summary.loc[pan_gene_id, 'uniq_COG_category_name'] = cog_dict[cog_cat]
        else:
            for col in df_eggnog.columns:
                df_pangene_summary.loc[pan_gene_id, col] = "-"
            df_pangene_summary.loc[pan_gene_id, 'COG_category_name'] = 'Not found in COG'
            df_pangene_summary.loc[pan_gene_id, 'uniq_COG_category_name'] = 'Not found in COG'
    
    df_pangene_summary.to_csv(eggnog_summary_path)
    # df_pangene_summary.to_csv('../data/genus/' + genus +'/'+ species + '/roary/df_pangene_eggnog_summary.csv')