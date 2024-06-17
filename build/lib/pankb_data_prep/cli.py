import argparse
import datetime
import logging
import os

from . import (
    isolation_source,
)

logging.basicConfig(
    # filename=log_filename,
    level=logging.DEBUG,
    format="%(asctime)s - %(levelname)s - %(message)s",
)

def main_isosource(args):
    isolation_source.find_isolation_source(args.genome_ids, args.output)

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
    }
    parsers = {}
    for x, f in modes.items():
        parsers[x] = subparsers.add_parser(x)
        parsers[x].set_defaults(func=f)

    parsers["isosource"].add_argument(
        "genome_ids",
        type=str,
        nargs="+",
        required=True,
        help="Genome IDs to process.",
    )
    parsers["isosource"].add_argument(
        "--output", "-o",
        type=str,
        required=True,
        help="Output file.",
    )
    # parser.add_argument(
    #     "--version", action="version", version="%(prog)s " + __version__
    # )

    args = parser.parse_args()
    args.func(args)
