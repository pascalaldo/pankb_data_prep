import pandas as pd
import json
import argparse


def initialize_parser(parser):
    parser.description = "Process data required for the landing page."
    parser.add_argument(
        "--species_summary",
        type=str,
        required=True,
        help="Species summary csv file.",
    )
    parser.add_argument(
        "--output_pankb_dimension",
        type=str,
        required=True,
        help="Gene count json file.",
    )
    parser.add_argument(
        "--output_json",
        type=str,
        required=True,
        help="Output in json format.",
    )


def generate_landing_page(
    species_summary_path, pankb_dimension_path, species_genome_gene_path
):
    df_species_summaries = pd.read_csv(species_summary_path, index_col=0)

    alleleome_done = ["Lactobacillaceae"]

    N_alleleome = 0
    Species_coding_alleleome = 0
    Species_coding_alleleome = []

    # for g in alleleome_done:
    #     Species_coding_alleleome += len(os.listdir('../data/genus/' + g))
    #     for s in os.listdir('../data/genus/' + g):
    #         gene_class = pd.read_csv('../data/genus/' + g + '/' + s + '/roary/df_pangene_summary_v2.csv')
    #         Core = len(gene_class.loc[gene_class['pangenome_class_2'] == 'Core',:])
    #         N_alleleome += Core

    gene_cluster_count = int(df_species_summaries["N_of_gene"].sum())
    genome_count = int(df_species_summaries["N_of_genome"].sum())
    species_count = int(df_species_summaries.shape[0])
    dimension = {
        "Genes": gene_cluster_count,
        "Coding alleleomes": N_alleleome,
        "Genomes": genome_count,
        "Species pangenomes": species_count,
        "Species coding alleleomes": Species_coding_alleleome,
    }
    with open(pankb_dimension_path, "w") as f:
        json.dump(dimension, f)

    family_genome_gene = {}
    for family in sorted(set(df_species_summaries["Family"].tolist())):
        family_genome_gene[family] = {}

        for species, row in df_species_summaries.loc[
            df_species_summaries["Family"] == family, :
        ].iterrows():
            family_genome_gene[family][species] = [
                int(row["N_of_genome"]),
                int(row["N_of_gene"]),
            ]

    with open(species_genome_gene_path, "w") as f:
        json.dump(family_genome_gene, f)


def run(args):
    generate_landing_page(
        args.species_summary,
        args.output_pankb_dimension,
        args.output_json,
    )


def main():
    parser = argparse.ArgumentParser()
    initialize_parser(parser)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
