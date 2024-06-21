import pandas as pd
import argparse


def initialize_parser(parser):
    parser.description = "Generate full summary table."
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


def full_summary_table(gtdb_meta_path, full_summary_path):
    df = pd.read_csv(gtdb_meta_path, low_memory=False, index_col=0).loc[
        :, ["Family", "gc_percentage", "genome_size"]
    ]
    df["source"] = "ncbi"
    df["gc_content"] = round((df["gc_percentage"]) * 0.01, 3)
    df["genome_len"] = df["genome_size"].astype(int)

    df.to_csv(full_summary_path)


def run(args):
    full_summary_table(args.gtdb_meta, args.output)


def main():
    parser = argparse.ArgumentParser()
    initialize_parser(parser)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
