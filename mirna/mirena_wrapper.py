#!/usr/bin/env python

"""
Initial documentation
"""

import sys
import re
import subprocess
import argparse

from Bio import SeqIO

def parse_arguments():
    """Parse sys.argv for arguments"""

    prog_description = (
        "This script is a multiprocess whole-genome wrapper for MIReNA, specifically for the"
        " functionality which searches for miRNAs within a text string. Output pre-miRNA fasta"
        " sequence headers are modified from the MIReNA default to include the genomic"
        " contig/scaffold ID and position of the found pre-miRNA."
        )
    parser = argparse.ArgumentParser(description=prog_description, add_help=False)
    args_required = parser.add_argument_group("Required")
    args_required.add_argument("-g", "--genome", required=True, dest="genome_fn", metavar="FASTA",
                               help="fasta file to search for miRNAs")
    args_required.add_argument("-m", "--mirna", required=True, dest="mirna_fn", metavar="FASTA",
                               help="fasta file containing miRNA sequences")
    
    args_optional = parser.add_argument_group("Optional")
    args_optional.add_argument("-h", "--help", help="show this help message and exit",
                               action="help")
    args_optional.add_argument("-p", "--n_procs", help="number of processors (default:" \
                               " %(default)s)", type=int, default=1)
    args_optional.add_argument("-t", "--temp_dir", help="directory to use for temporary files" \
                               " (default: %(default)s)", default="/tmp")
    args = parser.parse_args()
    return args

def main(args):
    do_something()
    # Read in genome; concatenate into individual long genome seq delimited by blocks of X's
    # Read in miRNAs; break into individual files of appropriate length
    # 

    # Future work:
    # Implement checkpointing; if some files are already done.

if __name__ == "__main__":
    args = parse_arguments()
    main(args)


# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4 textwidth=100 wrap
