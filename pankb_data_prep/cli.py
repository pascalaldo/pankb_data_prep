import argparse
import datetime
import logging
import os

from . import (
    isolation_source,
    eggnog_summary,
    summary_table,
    mash_list,
    heatmap,
    cog_data,
    heaps_law,
    gene_locustag,
    genome_page,
    landing_page,
)

logging.basicConfig(
    # filename=log_filename,
    level=logging.DEBUG,
    format="%(asctime)s - %(levelname)s - %(message)s",
)

def get_genome_list(args):
    with open(args.genomes, "r") as f:
        l = [genome.strip() for genome in f.readlines()]
    return l

def main_isosource(args):
    genomes = get_genome_list(args)
    isolation_source.find_isolation_source(genomes, args.output)

def main_eggnog(args):
    eggnog_summary.generate_eggnog_summary(
        args.gp_binary,
        args.summary,
        args.eggnog_table,
        args.reference,
        args.output,
    )

def main_family(args):
    summary_table.family_summary_table(args.name, args.gtdb_meta, args.output)

def main_species(args):
    summary_table.species_pangenome_summary(args.name, args.gp_binary, args.summary, args.gtdb_meta, args.output, args.output_json)

def main_all(args):
    summary_table.pangenome_summary(args.species_summary, args.output_json, args.output_genome_count, args.output_gene_count)

def main_mash(args):
    genomes = get_genome_list(args)
    mash_list.generate_mash_list(genomes, args.gtdb_meta, args.mash_file, args.output)

def main_heatmap(args):
    heatmap.generate_heatmap(
        args.gp_binary,
        args.gp_locustag,
        args.summary,
        args.eggnog_summary,
        args.mash_list,
        args.isosource,
        args.species_info,
        args.output_json,
    )

def main_cog(args):
    cog_data.generate_cog_data(
        args.eggnog_summary,
        args.gp_binary,
        args.output_json,
    )

def main_heaps(args):
    heaps_law.heaps_law(
        args.gp_binary,
        args.output_json,
    )

def main_locustag(args):
    gene_locustag.generate_locustag_data(
        args.gp_locustag,
        args.all_locustag,
        args.output,
    )

def main_genome(args):
    genome_page.generate_genome_page(
        args.species_summary,
        args.isosource,
        args.species_info,
        args.gp_binary,
        args.eggnog_summary,
        args.output,
    )

def main_landing(args):
    landing_page.generate_landing_page(
        args.species_summary,
        args.output_pankb_dimension,
        args.output_json,
    )
def ask_select_mode(args):
    logging.error("Please select a mode, see --help for more info.")


def main():
    logging.info("Application started")
    parser = argparse.ArgumentParser(
        description=(
            "PanKB data preparation."
        )
    )
    parser.set_defaults(func=ask_select_mode)
    subparsers = parser.add_subparsers()

    modes = {
        "isosource": main_isosource,
        "eggnog": main_eggnog,
        "family": main_family,
        "species": main_species,
        "all": main_all,
        "mash": main_mash,
        "heatmap": main_heatmap,
        "cog": main_cog,
        "heaps": main_heaps,
        "locustag": main_locustag,
        "genome": main_genome,
        "landing": main_landing,
    }
    parsers = {}
    for x, f in modes.items():
        parsers[x] = subparsers.add_parser(x)
        parsers[x].set_defaults(func=f)

    for x in ["family", "species"]:
        parsers[x].add_argument(
            "name",
            type=str,
            help="Family or analysis name to output data for.",
        )

    for x in ["isosource", "mash"]:
        parsers[x].add_argument(
            "--genomes", "-g",
            type=str,
            required=True,
            help="Genome IDs to process.",
        )
    for x in ["heatmap", "locustag"]:
        parsers[x].add_argument(
            "--gp_locustag",
            type=str,
            required=True,
            help="Gene presence locustag csv file.",
        )
    for x in ["eggnog", "species", "heatmap", "cog", "heaps", "genome"]:
        parsers[x].add_argument(
            "--gp_binary",
            type=str,
            required=True,
            help="Gene presence binary csv file.",
        )
    for x in ["eggnog", "species", "heatmap", "genome"]:
        parsers[x].add_argument(
            "--summary",
            type=str,
            required=True,
            help="Pangene summary csv file.",
        )
    parsers["eggnog"].add_argument(
        "--eggnog_table",
        type=str,
        required=True,
        help="Eggnog table file (.emapper.annotations).",
    )
    for x in ["heatmap", "cog", "genome"]:
        parsers[x].add_argument(
            "--eggnog_summary",
            type=str,
            required=True,
            help="Eggnog summary file.",
        )
    parsers["eggnog"].add_argument(
        "--reference",
        type=str,
        required=True,
        help="Pangenome reference fasta file.",
    )
    for x in ["isosource", "eggnog", "family", "species", "mash", "genome", "locustag"]:
        parsers[x].add_argument(
            "--output", "-o",
            type=str,
            required=True,
            help="Output file or directory.",
        )
    for x in ["species", "all", "heatmap", "cog", "heaps", "landing"]:
        parsers[x].add_argument(
            "--output_json",
            type=str,
            required=True,
            help="Output in json format.",
        )
    for x in ["family", "species", "mash"]:
        parsers[x].add_argument(
            "--gtdb_meta",
            type=str,
            required=True,
            help="GTDB meta csv file.",
        )
    for x in ["all", "genome", "landing"]:
        parsers[x].add_argument(
            "--species_summary",
            type=str,
            required=True,
            help="Species summary csv file.",
        )
    parsers["all"].add_argument(
        "--output_genome_count",
        type=str,
        required=True,
        help="Genome count json file.",
    )
    parsers["all"].add_argument(
        "--output_gene_count",
        type=str,
        required=True,
        help="Gene count json file.",
    )
    parsers["landing"].add_argument(
        "--output_pankb_dimension",
        type=str,
        required=True,
        help="Gene count json file.",
    )
    parsers["mash"].add_argument(
        "--mash_file",
        type=str,
        required=True,
        help="Mash file df_mash.csv.",
    )
    parsers["heatmap"].add_argument(
        "--mash_list",
        type=str,
        required=True,
        help="Mash list file.",
    )
    for x in ["heatmap", "genome"]:
        parsers[x].add_argument(
            "--isosource",
            type=str,
            required=True,
            help="Isolation source file.",
        )
        parsers[x].add_argument(
            "--species_info",
            type=str,
            required=True,
            help="Species info df_ncbi_meta csv file.",
        )
    parsers["locustag"].add_argument(
        "--all_locustag",
        type=str,
        required=True,
        help="Locustag info df_all_locustag csv file.",
    )
    # parser.add_argument(
    #     "--version", action="version", version="%(prog)s " + __version__
    # )

    args = parser.parse_args()
    args.func(args)
