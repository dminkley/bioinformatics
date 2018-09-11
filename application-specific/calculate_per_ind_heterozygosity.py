#!/usr/bin/env python

import argparse
import sys

VALID_GENOME_BASES = ['A', 'T', 'G', 'C']

def parse_arguments():
    """Parse sys.argv for arguments"""

    prog_description = "Given a VCF file, display genotype counts and, if a genome is specified, \
                       heterozygosity for each sample.  Heterozygosity is calculated as \
                       (# het sites)/(# valid genome bases); a valid genome base is a \
                       non-ambiguous (ie ATGCatgc) base.  Indels, and sites which are not \
                       bi-allelic are ignored."
    parser = argparse.ArgumentParser(description=prog_description, add_help=False)
    args_required = parser.add_argument_group("Required")
    args_required.add_argument("-i", "--in", required=True, dest="in_vcf_fn", metavar="INPUT_VCF",
                               help="VCF file for which stats are to be calculated")
    args_required.add_argument("-o", "--out", required=True, dest="out_info_fn",
                               metavar="OUTPUT_INFO", help="Output file with stats")
    
    args_optional = parser.add_argument_group("Optional")
    args_optional.add_argument("-h", "--help", help="show this help message and exit",
                               action="help")
    args_optional.add_argument("-g", "--genome", dest="genome_fn", metavar="FASTA", 
                               help="Fasta-format reference genome for this VCF file")
    args = parser.parse_args()
    return args


def main(args):

    if args.genome_fn is not None:
        sys.stderr.write("# A genome Fasta file was specified.  Calculating number of valid "
                         "bases...\n".format(args.in_vcf_fn))
        calculate_het = True
        n_valid_genome_bases = 0
        line_count = 0
        for line in open(args.genome_fn):
            # Output progress
            line_count += 1
            if line_count % 1000000 == 0:
                sys.stderr.write("# Processing fasta line: {}\n".format(line_count))
            if line[0] == ">":
                continue
            for base in line.rstrip():
                if base.upper() in VALID_GENOME_BASES:
                    n_valid_genome_bases += 1
        sys.stderr.write("# Number of non-ambiguous genome bases: "
                         "{}\n".format(n_valid_genome_bases))
    else:
        calculate_het = False
        
    sys.stderr.write("# Opening VCF file: {}\n".format(args.in_vcf_fn))
    fh_in_vcf = open(args.in_vcf_fn)
    
    line_count = 0
    for line in fh_in_vcf:

        # Output progress
        line_count += 1
        if line_count % 1000000 == 0:
            sys.stderr.write("# Processing VCF line: {}\n".format(line_count))

        # Skip header lines (except the table header)
        if line[:2] == "##":
            continue
        
        line_list = line.rstrip().split('\t')

        if line[0] == "#":
            sample_dict = {sample_id:{"0/0":0, "1/1":0, "0/1":0, "./.":0} for sample_id in
                           line_list[9:]}
            sample_list = [sample_id for sample_id in line_list[9:]]
            sys.stderr.write("# Number of samples detected: {}\n".format(len(sample_list)))
            continue


        # Check if the line is valid
        alt_allele_list = line_list[4].split(',')
        if len(alt_allele_list) > 1:
            #sys.stderr.write("# Position {} on scaffold {} is being skipped as it has more than "
            #                 "one alternate allele.\n".format(line_list[1], line_list[0]))
            continue
      
        if len(line_list[3]) > 1 or len(line_list[4]) > 1:
            #sys.stderr.write("# Position {} on scaffold {} is being skipped as it appears to be an "
            #                 "indel.\n".format(line_list[1], line_list[0]))
            continue

        # Add this line's genotypes to sample counts
        genotype_list = [ind_info.split(':')[0] for ind_info in line_list[9:]]
        for sample_id, genotype in zip(sample_list, genotype_list):
            sample_dict[sample_id][genotype] += 1

    # Output stats
    fh_out_info = open(args.out_info_fn, 'w')
    if calculate_het:
        fh_out_info.write("# Number of non-ambiguous genome bases: {}\n".format(n_valid_genome_bases))
        fh_out_info.write("# sample_id\tnum_hom_ref\tnum_hom_alt\tnum_het\tnum_uncalled\theterozygosity\n")
    else:
        fh_out_info.write("# sample_id\tnum_hom_ref\tnum_hom_alt\tnum_het\tnum_uncalled\n")

    for sample_id in sample_list:
        counts = sample_dict[sample_id]
        if calculate_het:
            heterozygosity = float(counts["0/1"])/float(n_valid_genome_bases)
            fh_out_info.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(sample_id, counts["0/0"],
                              counts["1/1"], counts["0/1"], counts["./."], heterozygosity))
        else:
            fh_out_info.write("{}\t{}\t{}\t{}\t{}\n".format(sample_id, counts["0/0"],
                              counts["1/1"], counts["0/1"], counts["./."]))

    # Print a friendly message :-)
    print "Done! :-)"


if __name__ == "__main__":
    args = parse_arguments()
    main(args)

# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4 textwidth=100 wrap
