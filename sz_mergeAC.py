'''
	poolseq_tk.py mergeAC
	Description: combining allele counts across samples/replicates
	Author: Simo V. Zhang

	Input: allele counts files

	Output: counts-merged file
			(1) chr
			(2) pos
			(3) ref base
			(4) alt base
			(5) merged allele counts
'''

import os
import sys
import collections
import argparse

import sz_utils
from colortext import ColorText

def run_merge(args):
	''' combine allele counts across replicates '''
	allele_counts = collections.defaultdict(list)
	data = collections.defaultdict(list)
	for ac_file in args.acs:
		sz_utils.check_if_files_exist(ac_file)
		ColorText().info("[poolseq_tk] reading and updating allele counts from %s ..."
						 %(ac_file), "stderr")
		with open(ac_file) as fAC:
			for line in fAC:
				tmp_line = line.strip().split()
				pos = int(tmp_line[1])
				if not pos in data:
					data[pos] = tmp_line[0:4]
				if not pos in allele_counts:
					allele_counts[pos] = map(int, tmp_line[4].split(':'))
				else:
					allele_counts[pos] = map(sum, zip(allele_counts[pos], map(int, tmp_line[4].split(':'))))
		ColorText().info(" [done]\n", "stderr")

	# output to file
	fOUT = None
	if args.out == sys.stdout:
		fOUT = sys.stdout
	else:
		sz_utils.make_dirs_if_necessary(args.out)
		fOUT = open(args.out, 'w')
	ColorText().info("[poolseq_tk] outputting to %s ..."
					 %(fOUT.name), "stderr")
	for pos in sorted(allele_counts.iterkeys()):
		fOUT.write("%s\t%s\n" %("\t".join(data[pos]), ":".join(map(str, allele_counts[pos]))))
	ColorText().info(" [done]\n", "stderr")
