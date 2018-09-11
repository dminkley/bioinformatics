#!/usr/bin/env python

class BlastHSP:
    
    def __init__(self, qry, sbj, perc_id, length, mismatch, gapopen, qry_start, qry_end, sbj_start,
            sbj_end, evalue, bitscore, strand):

        # Assumes values are already python-indexed and the correct type
        self.qry = qry
        self.sbj = sbj
        self.perc_id = perc_id
        self.length = length
        self.mismatch = mismatch
        self.gapopen = gapopen
        self.qry_start = qry_start
        self.qry_end = qry_end
        self.sbj_start = sbj_start
        self.sbj_end = sbj_end
        self.evalue = evalue
        self.bitscore = bitscore
        self.strand = strand

    @property
    def qry_len(self):
        return self.qry_end - self.qry_start

    @property
    def sbj_len(self):
        return self.sbj_end - self.sbj_start

    @property
    def blast_tbl_line(self):
        if self.strand == '+':
            line = '\t'.join(['{}']*12).format(self.qry, self.sbj, self.perc_id, self.length,
                self.mismatch, self.gapopen, self.qry_start+1, self.qry_end, self.sbj_start+1,
                self.sbj_end, self.evalue, self.bitscore)
        else:
            line = '\t'.join(['{}']*12).format(self.qry, self.sbj, self.perc_id, self.length,
                self.mismatch, self.gapopen, self.qry_start+1, self.qry_end, self.sbj_end+1,
                self.sbj_start, self.evalue, self.bitscore)
        return line
        
    def perc_qry_cov(self, qry_seq_len):
        return float(self.qry_len) / qry_seq_len * 100.0

    def perc_sbj_cov(self, sbj_seq_len):
        return float(self.sbj_len) / sbj_seq_len * 100.0

def read_table(fh_blast):
    # Generator of BlastHSPs
    for line in fh_blast:
        if line[0] == "#":
            continue
        line_list = line.rstrip().split('\t')
        qry, sbj = line_list[:2]
        perc_id = float(line_list[2])
        length, mismatch, gapopen = map(int, line_list[3:6])
        qry_start = min(int(line_list[6]), int(line_list[7])) - 1
        qry_end = max(int(line_list[6]), int(line_list[7]))
        sbj_start = min(int(line_list[8]), int(line_list[9])) - 1
        sbj_end = max(int(line_list[8]), int(line_list[9]))
        evalue = float(line_list[10])
        bitscore = float(line_list[11])
        strand = '+' if (int(line_list[9]) - int(line_list[8])) > 0 else '-'
        
        yield BlastHSP(qry, sbj, perc_id, length, mismatch, gapopen, qry_start, qry_end, sbj_start,
                sbj_end, evalue, bitscore, strand)

def main():
    print "In main"
    fh_blast = open("temp.blast")
    for i, hsp in enumerate(read_table(fh_blast)):
        print str(i) + "\t" + hsp.qry
    print "Done!"

if __name__=="__main__":
    main()



# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4 textwidth=100 wrap
