'''
	poolseq_tk.py filter
	Description: filter SNPs that do not satisfy the user-specified conditions
	Author: Simo Zhang

	Input: allele counts file generated from step 3

	Output: filtered set of SNPs with their allele counts
			(1) chr
			(2) pos
			(3) ref base
			(4) alt base
			(5) allele counts

'''

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
	before, after = 0, 0			# number of SNPs before and after filteration
	with open(args.ac_file, 'r') as fAC:
		for line in fAC:
			tmp_line = line.strip().split("\t")
			before += 1
			ref_base = tmp_line[2]
			alt_base  = tmp_line[3]
			counts = map(int, tmp_line[4].split(':'))
			if len(counts) < 2:
				ColorText().error("At least two counts (separated by colon) required"
								  "for column five\n")
				sys.exit(1)
			fail = 0
			for i in range(len(counts)):
				if i == 0:
					if counts[i] < args.min_ref_ac:
						fail = 1
						break
				elif i == 1:
					if (counts[i] < args.min_alt_ac or
					    counts[i] + counts[i-1] < args.min_cov):
						fail = 1
						break
				elif i == 2:
					if counts[i] < args.min_ref_ac:
						fail = 1
						break
				elif i == 3:
					if (counts[i] < args.min_alt_ac or
					   counts[i] + counts[i-1] < args.min_cov):
						fail = 1
						break
			if fail == 0:
				fOUT.write(line)
				after += 1
	fOUT.close()
	ColorText().info("Number of SNPs before filtering: %d\n" %(before), "stderr")
	ColorText().info("Number of SNPs after filtering: %d\n" %(after), "stderr")
