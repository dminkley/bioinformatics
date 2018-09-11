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

class SeqCollection:

    def __init__(self, file_id_obj):
        """ Given a file_id_obj of a type that can be used as input to Bio.SeqIO.parse(), create a
        genome object that can return a list of seq_ids in the order originally present in the file,
        as well as a dict object which accesses SeqRecord objects by id. """
        
        self.seq_dict = {}
        self.seq_list = []
        for seq_record in SeqIO.parse(file_id_obj, "fasta"):
            self.seq_dict[seq_record.id] = seq_record
            self.seq_list.append(seq_record)

    def __iter__(self):
        for seq_record in self.seq_list:
            yield seq_record

    @property
    def id_list(self):
        return [seq_record.id for seq_record in self.seq_list]

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
                               " (default: %(default)s)", dest="temp_dir", default="/tmp")
    args = parser.parse_args()
    return args

def create_genome_text_line(genome):
    """ Generate a single string from SeqRecord sequences in which each constituent sequence is
    separated from the others by a number (SPACER_WIDTH) of SPACER_CHAR bases"""
    
    spacer = SPACER_CHAR*SPACER_WIDTH
    return spacer.join(str(seq_record.seq) for seq_record in genome.seq_list)

def run_MIReNA_tasks(tmp_text_fn, miRNA_group_filenames, args.n_procs):
    mirena_args = [
                    "MIReNA.sh",
                    "--file "]


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
    genome = SeqCollection(args.genome_fn)
    miRNAs = SeqCollection(args.mirna_fn)

    # Create genome text file
    genome_as_text_line = create_genome_text_line(genome)
    tmp_text_fn = os.path.join(args.temp_dir, "temp_genome.txt")
    tmp_text_fh = open(tmp_text_fn, 'w')
    tmp_text_fh.write(genome_as_text_line)
    tmp_text_fh.close()
   
    # Create individual temp miRNA files
    miRNA_groups = [ [] for i in range(args.n_procs) ]
    for i, seq_record in enumerate(miRNAs):
        miRNA_groups[i % args.n_procs].append(seq_record)

    max_n_group_digits = len(str(len(miRNA_groups)))
    miRNA_group_filenames = []
    for group_num, seq_group in enumerate(miRNA_groups):
        base_fn = "temp_miRNA_group_{0:0>{width}}.fa".format(group_num, width=max_n_group_digits)
        miRNA_group_fn = os.path.join(args.temp_dir, base_fn)
        miRNA_group_filenames.append(miRNA_group_fn)
        SeqIO.write(seq_group, miRNA_group_fn, "fasta")

    # Now need to run the jobs themselves
    run_MIReNA_tasks(tmp_text_fn, miRNA_group_filenames, args.n_procs)

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
