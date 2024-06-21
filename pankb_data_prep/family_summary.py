import pandas as pd
import argparse


def initialize_parser(parser):
    parser.description = "Generate family summary."
    parser.add_argument(
        "name",
        type=str,
        help="Family or analysis name to output data for.",
    )
    parser.add_argument(
        "--gtdb_meta",
        type=str,
        required=True,
        help="GTDB meta csv file.",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        required=True,
        help="Output file or directory.",
    )


def family_summary_table(family, gtdb_meta_path, family_summary_path):
    if not family.startswith("f__"):
        family = f"f__{family}"
    df = pd.read_csv(gtdb_meta_path, low_memory=False, index_col=0).loc[
        :, ["Family", "gc_percentage", "genome_size"]
    ]
    df["source"] = "ncbi"
    df["gc_content"] = round((df["gc_percentage"]) * 0.01, 3)
    df["genome_len"] = df["genome_size"].astype(int)
    df = df.loc[df["Family"] == family, :]

    df.to_csv(family_summary_path)


def run(args):
    family_summary_table(args.name, args.gtdb_meta, args.output)


def main():
    parser = argparse.ArgumentParser()
    initialize_parser(parser)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
