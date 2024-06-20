import argparse
import pandas as pd
import os


def initialize_parser(parser):
    parser.description = "Process COG data."
    parser.add_argument(
        "--gp_binary",
        type=str,
        required=True,
        help="Gene presence binary csv file.",
    )
    parser.add_argument(
        "--eggnog_summary",
        type=str,
        required=True,
        help="Eggnog summary file.",
    )
    parser.add_argument(
        "--output_json",
        type=str,
        required=True,
        help="Output in json format.",
    )


def remove_slash(s):
    if "/" in str(s):
        return s.replace("/", "_")
    else:
        return s


def generate_cog_data(eggnog_summary_path, gp_binary_path, cog_data_path):
    eggnog_data = pd.read_csv(eggnog_summary_path, low_memory=False)
    # df_summary = pd.read_csv(summary_v2_path, low_memory=False)
    gp_binary = pd.read_csv(gp_binary_path, low_memory=False)
    gp_binary["Occurency"] = gp_binary.iloc[:, 1:].sum(axis=1)

    cog_class = pd.merge(eggnog_data, gp_binary.iloc[:, [0, -1]], on=["Gene"])

    selected_col_list = [
        "Gene",
        "COG_category",
        "COG_category_name",
        "Description",
        "Annotation",
        "PFAMs",
        "Occurency",
        "pangenome_class_2",
    ]
    df_all = cog_class.loc[:, selected_col_list]
    # df_core = cog_class.loc[cog_class.loc[:,"pangenome_class_2"] == "Core", selected_col_list]
    # df_accessory = cog_class.loc[cog_class.loc[:,"pangenome_class_2"] == "Accessory", selected_col_list]
    # df_rare = cog_class.loc[cog_class.loc[:,"pangenome_class_2"] == "Rare", selected_col_list]

    df_all.Gene = df_all.Gene.apply(remove_slash)

    df_all.to_json(cog_data_path, orient="values")


def run(args):
    generate_cog_data(
        args.eggnog_summary,
        args.gp_binary,
        args.output_json,
    )


def main():
    parser = argparse.ArgumentParser()
    initialize_parser(parser)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
