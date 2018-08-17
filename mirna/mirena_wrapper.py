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
    parser = argparse.ArgumentParser()
    parser.add_argument("genome_file", help="fasta file to search for miRNAs")
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

