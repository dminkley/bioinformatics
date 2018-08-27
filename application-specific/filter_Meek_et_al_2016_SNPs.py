#!/usr/bin/env python

import sys
import argparse
import itertools

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class VariantLocus:
    # If I was really doing this, and wanted to make something useful for all my future work or to
    # put into a library, I would make a VariantLocus that is extensible.  I would then make a more
    # specific SNPLocus or something like that inherits from VariantLocus.

    # Object variables
    # variant_id
    # samples
    # allele1_region_seq
    # allele2_region_seq
    # variant_pos_in_seq
    # allele1_base
    # allele2_base

    def __init__(self, variant_id):
        """ A class to represent single-nucleotide variants and their properties across multiple
        individuals.  Assumes bi-allelic variants only. """

        self.variant_id = variant_id
        self.samples = {}

    def add_sample_read_counts(self, sample_id, allele_read_counts):
        """ Given a string of the format 'num_allele1_reads,num_allele2_reads', add the individual
        (and its read counts) to this variant """

        num_allele1_reads, num_allele2_reads = allele_read_counts.split(',')
        self.samples[sample_id] = map(int, [num_allele1_reads, num_allele2_reads])

    def read_coverage_by_sample(self):
        """ Return a dict that reports sample_id:total_reads for each sample at this variant
        locus """

        return {sample_id: sum(allele_counts) for sample_id, allele_counts in self.samples.items()}

    def add_seq_variants(self, allele1_seq, allele2_seq):
        """ Add two sequences to this locus that correspond to the variant itself and the
        surrounding bases. Also determine where the variant is located in the sequence and set
        relevent parameters. """

        self.allele1_region_seq = allele1_seq
        self.allele2_region_seq = allele2_seq

        bases_different = [base1!=base2 for base1, base2 in zip(allele1_seq.seq, allele2_seq.seq)]
        variant_base_positions = [pos for pos, different in enumerate(bases_different) if
                different]

        if len(variant_base_positions) != 1:
            print str(variant_base_positions)
            print str(allele1_seq.seq)
            print str(allele2_seq.seq)
            print "ERROR: Two variant sequences do not have only a single difference"
            print "The sequences are: " + allele1_seq.id + " and " + allele2_seq.id
            sys.exit(-1)

        self.variant_pos_in_seq = variant_base_positions[0]
        self.allele1_base = allele1_seq[self.variant_pos_in_seq]
        self.allele2_base = allele2_seq[self.variant_pos_in_seq]

        template_seq = "{}{}{}".format(allele1_seq.seq[:self.variant_pos_in_seq], 'N',
                allele1_seq.seq[self.variant_pos_in_seq+1:])
        self.seq_region = SeqRecord(Seq(template_seq), id="{}".format(self.variant_id), description="")

def filter_variants(variant, min_read_depth=20, num_good_samples=50):
    """ Given a VariantLocus object, return True if it has at least num_good_samples samples that
    have at least a combined (allele1 + allele2) read depth of min_read_depth. """

    # Get number of samples which have the minimum read depth
    good_samples = [sample_id for sample_id, total_reads in
            variant.read_coverage_by_sample().items() if
            total_reads >= min_read_depth]
    return len(good_samples) >= num_good_samples

def parse_arguments():
    """Parse sys.argv for arguments"""

    prog_description = (
        "Requires a table of RADseq-like variant markers which includes columns corresponding to "
        "marker loci and rows corresponding to samples/individuals, with cell values consisting of "
        "<allele1>,<allele2> read depth counts, as well as a fasta file with allele sequences "
        "(with neighbouting regions).  Filters the variants to include only those which reach some "
        "total variant read depth in some number of samples.  A fasta file containing the retained "
        "variants is produced, with one sequence per variant and an ambiguous 'N' base at the "
        "position of the variant base itself."
        )
    parser = argparse.ArgumentParser(description=prog_description, add_help=False)
    args_required = parser.add_argument_group("Required")
    args_required.add_argument("-t", "--table", required=True, dest="SNP_counts_fn", metavar="TABLE",
                               help="allele reads depths by variant and sample")
    args_required.add_argument("-s", "--seq", required=True, dest="SNP_seq_fn", metavar="FASTA",
                               help="fasta file with variant and neighbouring sequence, for each \
                               allele")
    args_required.add_argument("-o", "--out", required=True, dest="out_fn", metavar="FASTA",
                               help="output fasta file with retained variants.")
    
    args_optional = parser.add_argument_group("Optional")
    args_optional.add_argument("-h", "--help", help="show this help message and exit",
                               action="help")
    args_optional.add_argument("-d", "--depth", help="read depth required at a variant locus in a \
                               sample for that sample to contribute to the total number of good \
                               samples at a variant locus (default: %(default)s)", type=int,
                               default=10, dest="min_read_depth")
    args_optional.add_argument("-n", "--n_samples", help="minimum number of samples with \
                               sufficient read depth required for a variant locus to be \
                               retained (default: %(default)s)", type=int, default=50,
                               dest="num_good_samples")
    args = parser.parse_args()
    return args

def main(args):

    # Process variant/allele read count file
    fh_SNP_counts = open(args.SNP_counts_fn)
    total_snp_coverage = {}
    first_line = fh_SNP_counts.readline().strip()
    variants = []
    variant_lookup = {}
    for marker_id in first_line.split(' ')[3:]:
        variant = VariantLocus(marker_id)
        variants.append(variant)
        variant_lookup[marker_id] = variant

    for line in fh_SNP_counts:
        line_list = line.rstrip().split(' ')
        sample_id = line_list[0]
        for allele_read_counts, variant in itertools.izip(line_list[3:], variants):
            variant.add_sample_read_counts(sample_id, allele_read_counts)
        
    # Add seq information to variants.  Assumes that allele1 and allele2 sequences for a given
    # marker are adjacent in the fasta file.
    variant_allele_seqs = list(SeqIO.parse(args.SNP_seq_fn, "fasta"))
    for allele1_seq, allele2_seq in itertools.izip(variant_allele_seqs[::2],
            variant_allele_seqs[1::2]):
        allele1_marker_id = allele1_seq.id.split('_')[0]
        allele2_marker_id = allele2_seq.id.split('_')[0]
        if allele1_marker_id != allele2_marker_id:
            print "ERROR! Two successive alleles should have the same marker ID, but they dont!"
            print "The marker ids are: " + allele1_marker_id + " and " + allele2_marker_id
        variant_lookup[allele1_marker_id].add_seq_variants(allele1_seq, allele2_seq)

    # Identify high-quality variants
    hq_variants = filter(lambda x: filter_variants(x, args.min_read_depth, args.num_good_samples),
                         variants)

    # TEMP: Output hq_variants
    hq_variant_seqs = [variant.seq_region for variant in hq_variants]
    SeqIO.write(hq_variant_seqs, args.out_fn, "fasta")

    print "Done! :-)"


if __name__ == "__main__":
    #fh_SNP_counts = open("Meek_et_al_2016_SNP_counts_in_individuals.txt")
    #fh_SNP_seq = open("Meek_et_al_2016_SNP_regions.fa")
    #fh_outfile = open("high_quality_variants.fa", 'w')
    args = parse_arguments()
    main(args)

# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4 textwidth=100 wrap
