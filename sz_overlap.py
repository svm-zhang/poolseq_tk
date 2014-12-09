'''
	python poolseq_tk.py overlap
	Description: Getting SNPs identified from both pools/samples
	Author: Simo V. Zhang

	Input: statistical tests generated from step 5
	Output: set of SNPs in both
'''

import os
import sys
import argparse
import collections

import sz_utils
from colortext import ColorText

def run_overlap(args):
	''' getting SNPs identified from both pools '''
	sz_utils.check_if_files_exist(args.file_a, args.file_b)

	snp_a = collections.defaultdict(list)
	with open(args.file_a, 'r') as fA:
		for line in fA:
			tmp_line = line.strip().split("\t")
			snp_a[int(tmp_line[1])] = tmp_line
	ColorText().info("[poolseq_tk]: %d SNPs parsed from %s\n" %(len(snp_a),
					 os.path.basename(args.file_a)), "stderr")

	sz_utils.make_dirs_if_necessary(args.out)

	num_overlapion = 0
	with open(args.out, 'w') as fOUT:
		with open(args.file_b, 'r') as fB:
			for line in fB:
				tmp_line = line.strip().split("\t")
				if int(tmp_line[1]) in snp_a:
					num_overlapion += 1
					fOUT.write("%s\t%s\n" %("\t".join(snp_a[int(tmp_line[1])]), "\t".join(tmp_line[-4:])))
	ColorText().info("[poolseq_tk]: %d SNPs identified from both pools\n" %(num_overlapion),
					 "stderr")

