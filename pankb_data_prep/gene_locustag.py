import pandas as pd
from Bio import AlignIO
from Bio import SeqIO
import json
from pathlib import Path
import argparse


def initialize_parser(parser):
    parser.description = "Process gene-locustag relation info."
    parser.add_argument(
        "--gp_locustag",
        type=str,
        required=True,
        help="Gene presence locustag csv file.",
    )
    parser.add_argument(
        "--all_locustag",
        type=str,
        required=True,
        help="Locustag info df_all_locustag csv file.",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        required=True,
        help="Output file or directory.",
    )


def remove_slash(s):
    if "/" in str(s):
        return s.replace("/", "_")
    else:
        return s


def generate_locustag_data(gp_locustag_path, all_locustag_path, gene_locustag_dir):
    gene_locustag_dir = Path(gene_locustag_dir)
    gene_locustag_dir.mkdir(parents=True, exist_ok=True)
    df_gene_presence_locustag = pd.read_csv(
        gp_locustag_path, index_col="Gene", low_memory=False
    )
    df_gene_presence_locustag.index = [
        remove_slash(i) for i in list(df_gene_presence_locustag.index)
    ]
    all_locustag_df = pd.read_csv(all_locustag_path, index_col=0, low_memory=False)

    for gene_id in df_gene_presence_locustag.index.tolist():
        gene_locustag = []
        for genome_id, locus_tag in (
            df_gene_presence_locustag.loc[gene_id, :].dropna()
        ).items():
            if "\t" in locus_tag:
                locus_tag_list = locus_tag.split("\t")
                for locus_tag_id in locus_tag_list:
                    gene_locustag.append(f"{genome_id}@{locus_tag_id}")
            else:
                gene_locustag.append(f"{genome_id}@{locus_tag}")

        # Write the JSON object to a file
        with open(gene_locustag_dir / f"{gene_id}.json", "w") as f:
            json.dump(
                all_locustag_df.loc[gene_locustag, :].to_dict(orient="records"),
                f,
                separators=(",", ":"),
                ensure_ascii=False,
            )


def run(args):
    generate_locustag_data(
        args.gp_locustag,
        args.all_locustag,
        args.output,
    )


def main():
    parser = argparse.ArgumentParser()
    initialize_parser(parser)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
