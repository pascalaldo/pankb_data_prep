import pandas as pd
import json
import argparse


def initialize_parser(parser):
    parser.description = "Generate summary."
    parser.add_argument(
        "--species_summary",
        type=str,
        required=True,
        help="Species summary csv file.",
    )
    parser.add_argument(
        "--output_json",
        type=str,
        required=True,
        help="Output in json format.",
    )
    parser.add_argument(
        "--output_genome_count",
        type=str,
        required=True,
        help="Genome count json file.",
    )
    parser.add_argument(
        "--output_gene_count",
        type=str,
        required=True,
        help="Gene count json file.",
    )


def pangenome_summary(
    species_summary_path,
    species_list_out_path,
    genome_count_out_path,
    gene_count_out_path,
):
    df = pd.read_csv(species_summary_path)

    json_data = {
        "Family": df["Family"].tolist(),
        "Species": df["Species"].tolist(),
        "Openness": df["Openness"].tolist(),
        "N_of_genome": df["N_of_genome"].astype(int).tolist(),
        "Pangenome_analyses": df["Pangenome_analyses"].tolist(),
        "Gene_class": list(
            map(
                list,
                zip(
                    *[
                        df[x].astype(int).tolist()
                        for x in ["N_of_core", "N_of_accessory", "N_of_rare"]
                    ]
                ),
            )
        ),
    }
    with open(species_list_out_path, "w") as f:
        json.dump(json_data, f)

    family_group = df.groupby("Family")
    family_group["N_of_genome"].sum().sort_values(ascending=False).to_json(
        genome_count_out_path, orient="index"
    )

    family_group["N_of_gene"].sum().sort_values(ascending=False).to_json(
        gene_count_out_path, orient="index"
    )


def run(args):
    pangenome_summary(
        args.species_summary,
        args.output_json,
        args.output_genome_count,
        args.output_gene_count,
    )


def main():
    parser = argparse.ArgumentParser()
    initialize_parser(parser)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
