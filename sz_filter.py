import os
import sys
import argparse

import sz_utils
from colortext import ColorText

def run_filter(args):
	sz_utils.check_if_files_exist(args.ac_file)
	fOUT = None
	if args.out == sys.stdout:
		fOUT = sys.stdout
	else:
		sz_utils.make_dirs_if_necessary(args.out)
		fOUT = open(args.out, 'w')
	before, after = 0, 0
	with open(args.ac_file, 'r') as fAC:
		for line in fAC:
			tmp_line = line.strip().split("\t")
			before += 1
			ref_base = tmp_line[2]
			alt_base  = tmp_line[3]
			counts = map(int, tmp_line[4].split(':'))
			if (counts[1] >= args.min_alt_ac and
				counts[0] >= args.min_ref_ac and
				counts[2] >= args.min_ref_ac and
				counts[3] >= args.min_alt_ac and
				counts[0]+counts[1] >= args.min_cov and
				counts[2]+counts[3] >= args.min_cov):
				fOUT.write(line)
				after += 1
	fOUT.close()
	ColorText().info("Number of SNPs before filtering: %d\n" %(before), "stderr")
	ColorText().info("Number of SNPs after filtering: %d\n" %(after), "stderr")
