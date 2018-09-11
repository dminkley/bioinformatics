#!/usr/bin/env python

import argparse
import re

from Bio import SeqIO

def parse_arguments():
    """Parse sys.argv for arguments"""

    prog_description = "Given a VCF file filter based on two INFO fields: MAF and CAT."
    parser = argparse.ArgumentParser(description=prog_description, add_help=False)
    args_required = parser.add_argument_group("Required")
    args_required.add_argument("-i", "--in", required=True, dest="in_vcf_fn", metavar="INPUT_VCF",
                               help="VCF file to be filtered")
    args_required.add_argument("-o", "--out", required=True, dest="out_vcf_fn", metavar="OUTPUT_VCF",
                               help="output, filtered VCF file")
    
    args_optional = parser.add_argument_group("Optional")
    args_optional.add_argument("-h", "--help", help="show this help message and exit",
                               action="help")
    args_optional.add_argument("-m", "--min_maf", help="minimum minor allele frequency required in \
                               order for a variant to be retained (default: %(default)s)", type=float,
                               default=0, dest="min_maf")
    args_optional.add_argument("-c", "--categories", help="comma-delimited list of categories \
                               (eg \"PolyHighResolution\" or \"NoMinoHom\"), variants from which \
                               are to be removed (default: \"%(default)s\")", default="",
                               dest="cats_to_remove")
    args = parser.parse_args()
    return args


def main(args):

    fh_in_vcf = open(args.in_vcf_fn)
    fh_out_vcf = open(args.out_vcf_fn, 'w')
    
    vcf_data_re = re.compile(r"(?P<chrom>[\w\.]+)\t(?P<pos>\d+)\t(?P<id>\S+)\t(?P<ref>\S+)\t"
                             r"(?P<alt>\S+)\t(?P<qual>\S+)\t(?P<filter>\S+)\t(?P<info>\S*)")
    info_re = re.compile(r"MAF=(?P<maf>\d+\.*\d*);CAT=(?P<cat>[A-Za-z]+).*")

    # Process list of categories to remove
    cats_to_remove = {cat:True for cat in args.cats_to_remove.split(',')}

    for line in fh_in_vcf:
        if line[0] == "#":
            fh_out_vcf.write(line)
        else:
            vcf_line_m = vcf_data_re.match(line)
            info = vcf_line_m.group("info")
            info_field_m = info_re.match(info)
            if (float(info_field_m.group("maf")) < args.min_maf or 
                info_field_m.group("cat") in cats_to_remove):
                continue
            else:
                fh_out_vcf.write(line)

    print "Done! :-)"

            
if __name__ == "__main__":
    args = parse_arguments()
    main(args)

# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4 textwidth=100 wrap
