#!/usr/bin/env python

import argparse
import re
import datetime

from Bio import SeqIO

def parse_arguments():
    """Parse sys.argv for arguments"""

    prog_description = (
        "Given a scaffold index file produced by ScaffoldStitcher and a VCF with coordinates "
        "based on the original (un-stitched) genome, output a VCF file with coordinates which "
        "refer to the same loci in the scaffold-stitched genome." 
        )
    parser = argparse.ArgumentParser(description=prog_description, add_help=False)
    args_required = parser.add_argument_group("Required")
    args_required.add_argument("-i", "--index", required=True, dest="scaffold_index_fn",
                               metavar="INDEX_FILE", help="scaffold index file produced by \
                               ScaffoldStitcher")
    args_required.add_argument("-v", "--vcf", required=True, dest="in_vcf_fn", metavar="INPUT_VCF",
                               help="VCF file with coordinates based on the original (unstitched) \
                               genome sequence")
    args_required.add_argument("-g", "--genome", required=True, dest="genome_fn", metavar="FASTA",
                               help="original (unstitched) genome sequence fasta file")
    args_required.add_argument("-r", "--ref_name", required=True, dest="stitched_ref_fn", 
                               metavar="FILENAME", help="file name of the stitched reference file \
                               to be used in the VCF header")
    args_required.add_argument("-o", "--out", required=True, dest="out_vcf_fn", metavar="OUTPUT_VCF",
                               help="VCF file with coordinates based on the stitched genome \
                               sequence")
    
    args_optional = parser.add_argument_group("Optional")
    args_optional.add_argument("-h", "--help", help="show this help message and exit",
                               action="help")
    args_optional.add_argument("-s", "--spacer", help="number of 'N' spacer bases between \
                               original scaffolds in a superscaffold (default: %(default)s)",
                               type=int, default=1000, dest="spacer_len")
    
    args = parser.parse_args()
    return args

class ConcatScaffoldMap:

    def __init__(self, fh_scaffold_index, genome_seq_dict, spacer_len):
        """ Create an object which coordinates between stitched and unstitched sequences from a
        ScaffoldStitcher scaffold index file. """

        curr_superscaff = ""
        self.map = {}
        used_scaffs = {}
        for line in fh_scaffold_index:
            m = re.match(r">(\S+)\t(\d+)\t>(\S+).*$", line)
            superscaff_id = m.group(1)
            order_in_superscaff = int(m.group(2))
            orig_scaff_id = m.group(3)
            
            if superscaff_id != curr_superscaff:
                curr_pos = 0
                curr_superscaff = superscaff_id
            
            self.map[orig_scaff_id] = [superscaff_id, curr_pos]
            curr_pos = curr_pos + len(genome_seq_dict[orig_scaff_id]) + spacer_len
            used_scaffs[orig_scaff_id] = True
        
        # Add any old scaffolds/chromosomes which were not added to superscaffolds
        for orig_scaff_id in genome_seq_dict:
            if orig_scaff_id not in used_scaffs:
                self.map[orig_scaff_id] = [orig_scaff_id, 0]

    def map_old_to_new(self, old_scaff_id, old_pos):
        """ Given an old scaffold id and position, return the corresponding position in the sequence
        create by ScaffoldStitcher """

        new_superscaff_id, start_pos_in_superscaff = self.map[old_scaff_id]

        return new_superscaff_id, start_pos_in_superscaff + old_pos



def main(args):
    
    # Set up the coordinate map
    fh_scaffold_index = open(args.scaffold_index_fn)
    genome_seq_dict = SeqIO.to_dict(SeqIO.parse(args.genome_fn, "fasta"))

    coord_map = ConcatScaffoldMap(fh_scaffold_index, genome_seq_dict, args.spacer_len)

    # Read the old vcf, output the new vcf with mapped coordinates.
    fh_in_vcf = open(args.in_vcf_fn)
    
    header_re = re.compile(r"##(\w+)=(.+)$")
    vcf_data_re = re.compile(r"(?P<chrom>[\w\.]+)\t(?P<pos>\d+)\t(?P<id>\S+)\t(?P<ref>\S+)\t"
                             r"(?P<alt>\S+)\t(?P<qual>\S+)\t(?P<filter>\S+)\t(?P<info>\S*)")

    now = datetime.datetime.now()
    fh_outfile = open(args.out_vcf_fn, 'w')
    for line in fh_in_vcf:
        re_match = header_re.match(line)
        if re_match:
            key = re_match.group(1)
            value = re_match.group(2)
            if key == "fileDate":
                out_str = "##fileDate={}\n".format(now.strftime("%Y%m%d"))
            elif key == "reference":
                out_str = "##reference={}\n".format(args.stitched_ref_fn)
            else:
                out_str = line
        else:
            if line[:6] == "#CHROM":
                out_str = line
            else:
                vcf_data = vcf_data_re.match(line)
                new_chrom, new_pos = coord_map.map_old_to_new(vcf_data.group("chrom"),
                                                              int(vcf_data.group("pos"))-1)
                out_str = "{}\t{}\t{}\n".format(new_chrom, new_pos+1,'\t'.join(vcf_data.groups()[2:]))
        fh_outfile.write(out_str)

    print "Done! :-)"

            
if __name__ == "__main__":
    args = parse_arguments()
    main(args)

# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4 textwidth=100 wrap
