import argparse
import datetime
import logging
import os

from . import (
    isolation_source,
    eggnog_summary,
    summary_table,
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
    summary_table.family_summary_table(args.family, args.gtdb_meta, args.output)

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
    }
    parsers = {}
    for x, f in modes.items():
        parsers[x] = subparsers.add_parser(x)
        parsers[x].set_defaults(func=f)

    parsers["family"].add_argument(
        "family",
        type=str,
        help="Family name to output data for.",
    )

    parsers["isosource"].add_argument(
        "--genomes", "-g",
        type=str,
        required=True,
        help="Genome IDs to process.",
    )
    parsers["eggnog"].add_argument(
        "--gp_binary",
        type=str,
        required=True,
        help="Gene presence binary csv file.",
    )
    parsers["eggnog"].add_argument(
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
    parsers["eggnog"].add_argument(
        "--reference",
        type=str,
        required=True,
        help="Pangenome reference fasta file.",
    )
    for x in ["isosource", "eggnog"]:
        parsers[x].add_argument(
            "--output", "-o",
            type=str,
            required=True,
            help="Output file.",
        )
    parsers["family"].add_argument(
        "--gtdb_meta",
        type=str,
        required=True,
        help="GTDB meta csv file.",
    )
    # parser.add_argument(
    #     "--version", action="version", version="%(prog)s " + __version__
    # )

    args = parser.parse_args()
    args.func(args)
