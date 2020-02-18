#!/usr/bin/env python

import pysam
import argparse
from collections import deque

GT_TYPES = ["homR", "homA", "het", "uncalled"]

def parse_arguments():
    return args

class Window:

    def __init__(self, contig, start, end, win_num):
        self.contig = contig
        self.start = start
        self.end = end
        self.win_num = win_num


class GTSurveyWindow(Window):

    def __init__(self, contig, start, end, win_num, sample_ids):
        super().__init__(contig, start, end, win_num)
        self.homR = {sample_id: 0 for sample_id in sample_ids}
        self.homA = {sample_id: 0 for sample_id in sample_ids}
        self.het = {sample_id: 0 for sample_id in sample_ids}
        self.uncalled = {sample_id: 0 for sample_id in sample_ids}

def gt_survey_windows(vcf_in, size, step):

    win_nums_in_q_start = 0
    win_nums_in_q_end = 0    # Note that this is half-open
    sample_ids = list(vcf_in.header.samples)
    open_windows = deque()


    curr_contig = ""
    for var in vcf_in:
        
        if var.config != curr_contig:
            # New contig (or first contig)
            curr_contig = var.contig

            # Dump (yield) all open windows (from previous contig) to begin reset for next contig
            for win_num in range(win_nums_in_q_start, win_nums_in_q_end:
                win_start = win_num * step
                win_end = win_start + size
                finished_win = open_windows.popleft()
                assert win_start == finished_win.start
                assert win_end == finished_win.end
                yield finished_win
            
            # Double check that window queue is empty and reset trackers
            assert len(open_windows) == 0
            win_nums_in_q_start = 0
            win_nums_in_q_end = 0


        # First, update which windows are active.  Return any that are no longer active.
        # Which window numbers include this variant?
        this_var_win_nums_start = ((var.pos - size) // step) + 1
        this_var_win_nums_end = (var.pos // step) + 1

        # Add newly opened windows as necessary
        for win_num in range(win_nums_in_q_end, this_var_win_nums_end):
            win_start = win_num * step
            win_end = win_start + size      # Note these produce 0-based half-open coords
            new_window = GTSurveyWindow(contig, win_start, win_end, win_num, sample_ids)
            open_windows.append(new_window)
        win_nums_in_q_end = this_var_win_nums_end

        # Remove windows that we are no longer using.  Yield their contents
        for win_num in range(win_nums_in_q_start, this_var_win_nums_start):
            win_start = win_num * step
            win_end = win_start + size
            finished_win = open_windows.popleft()

            # The window we dequeue should match with our calculated window number info
            assert win_start == finished_win.start
            assert win_end == finished_win.end

            yield finished_win
        win_nums_in_q_start = this_var_win_nums_start

        # Add the stats for this variant to all of the currently open windows
        for window in open_windows:
            for sample_id in sample_ids:
                gt = var.samples[sample_id]["GT"]
                if gt[0] == None:   # Uncalled
                    window.uncalled[sample_id] += 1
                elif gt[0] == gt[1] == 0:
                    window.homR[sample_id] += 1
                elif gt[0] == gt[1] == 1:
                    window.homA[sample_id] += 1
                elif gt[0] == 0 and gt[1] == 1:
                    window.het[sample_id] += 1
                else:
                    raise ValueError("Unexpected genotype ({}, {}) for sample_id: {} on \
                                     chromosome {} at position {}".format(gt[0], gt[1], 
                                         sample_id, var.contig, var.pos))
    
    # Dump (yield) all open windows and finish
    for win_num in range(win_nums_in_q_start, win_nums_in_q_end):
        win_start = win_num * step
        win_end = win_start + size
        finished_win = open_windows.popleft()
        assert win_start == finished_win.start
        assert win_end == finished_win.end
        yield finished_win

    # Finally check q is empty and head home
    assert len(open_windows) == 0
    return


def main(args):
    vcf_in = pysam.VariantFile(args.fn_in_vcf)
    sample_ids = list(vcf_in.header.samples)
    

    fh_homR_out = open("{}.homR.txt", 'w')
    fh_homA_out = open("{}.homA.txt", 'w')
    fh_het_out = open("{}.het.txt", 'w')
    fh_uncalled_out = open("{}.uncalled.txt", 'w')

    header = "contig\tstart\tend\t{}\n".format('\t'.join(sample_ids))
    fh_homR_out.write(header)
    fh_homA_out.write(header)
    fh_het_out.write(header)
    fh_uncalled_out.write(header)

    for gt_window in gt_survey_windows(vcf_in, size=args.win_size, step=args.win_step):
        contig = gt_window.contig
        start = gt_window.start
        end = gt_window.end

        fh_homR.out.write("{}\t{}\t{}\t{}\n".format(contig, start, end, '\t'.join([gt_window.homR[sample_id] for sample_id in sample_ids])))
        fh_homA.out.write("{}\t{}\t{}\t{}\n".format(contig, start, end, '\t'.join([gt_window.homA[sample_id] for sample_id in sample_ids])))
        fh_het.out.write("{}\t{}\t{}\t{}\n".format(contig, start, end, '\t'.join([gt_window.het[sample_id] for sample_id in sample_ids])))
        fh_uncalled.out.write("{}\t{}\t{}\t{}\n".format(contig, start, end, '\t'.join([gt_window.uncalled[sample_id] for sample_id in sample_ids])))

    echo "Done! :-)"

if __name__ == "__main__":
    args = parse_arguments()
    main(args)

