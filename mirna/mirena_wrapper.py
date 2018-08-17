#!/usr/bin/env python

"""
Initial documentation
"""

import sys
import re
import subprocess
import argparse

from Bio import SeqIO


def main():
    prog_description = (
        "This script is a multiprocess whole-genome wrapper for MIReNA,"
        " specifically for the functionality which searches for miRNAs within"
        " a text string. Output pre-miRNA fasta sequence headers are modified"
        " from the MIReNA default to include the genomic contig/scaffold ID"
        " and position of the found pre-miRNA."
        )
    
    parser = argparse.ArgumentParser(description=prog_description)
    parser.add_argument("genome_file", metavar="METAVAR_TEST", help="fasta file to search for miRNAs")
    parser.add_argument("mirna_file",
                        help="fasta file containing miRNA sequences")
    parser.add_argument("--n_procs", "-p", help="number of processors",
                        type=int, default=1)
    parser.add_argument("--temp_dir", "-t", help="directory to use for " \
                        "temporary files", default="/tmp")
    args = parser.parse_args()

if __name__ == "__main__":
    main()


# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4 textwidth=80 wrap

