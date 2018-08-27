#!/usr/bin/env python

"""
Initial documentation
"""

import sys
import re
import subprocess
import argparse

from Bio import SeqIO

SPACER_WIDTH=1000
SPACER_CHAR='X'

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

def generate_single_genome_str(genome_list):
    """ Generate a single string from SeqRecord sequences in which each constituent sequence is
    separated from the others by a number (SPACER_WIDTH) of SPACER_CHAR bases"""
    
    spacer = SPACER_CHAR*SPACER_WIDTH
    return spacer.join(str(seq_record.seq) for seq_record in genome_list)


def create_temp_files():
    # Content here

def run_parallel_jobs():
    # Content here

def merge_results():
    # Content here

# Create miRNA group class?  Could track temp associated temp file, the miRNAs themselves, etc

def main(args):
    
    # Read in genome
    # Need a function to create all temporary files?
    genome_list = [seq_record for seq_record in SeqIO.parse(args.genome, "fasta")]
    genome_text = generate_single_genome_str(genome_list)
   
    # Write file could have a check to see if the file exists?

    # Create temporary files
        # Single genome file
        # Multiple miRNA files

    # Perform processing
        

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
