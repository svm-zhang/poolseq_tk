'''
	python poolseq_tk.py count

	Description: Count alleles at each SNP give the pileups
	Author: Simo V. Zhang

	Input: pileup file with reads bases converted to corresponding alleles

	Output: pielup file with allele counts
		   (1) chr
		   (2) pos
		   (3) ref base
		   (4) alt base
		   (5) allele counts in the order of ref and alt, separated by colon

'''

import collections
import sys
import os

import sz_utils
from colortext import ColorText

def run_count(args):
	''' Counting alleles at each SNP in the given pileup files '''
	ac = collections.defaultdict(tuple)
	for pileup in args.pileups:
		sz_utils.check_if_files_exist(pileup)
		nsnps = 0
		ColorText().info("[poolseq_tk] counting alleles in %s:" %(os.path.basename(pileup)), "stderr")
		with open(pileup, 'r') as fMPILEUP:
			for line in fMPILEUP:
				nsnps += 1
				tmp_line = line.strip().split("\t")
				pos = int(tmp_line[1])
				ref_base = tmp_line[2]
				alt_base = tmp_line[3]
				ref_count, alt_count = 0, 0
				if tmp_line[-1] != "N/A":
					ref_count = tmp_line[-1].count(ref_base) + \
								tmp_line[-1].count(ref_base.lower())
					alt_count = tmp_line[-1].count(alt_base) + \
								tmp_line[-1].count(alt_base.lower())
				if pos not in ac:
					ac[pos] = [tmp_line[0], ref_base, alt_base, str(ref_count), str(alt_count)]
				else:
					ac[pos] += [str(ref_count), str(alt_count)]
		ColorText().info(" %d SNPs parsed\n" %(nsnps), "stderr")

	fOUT = None
	if args.out == sys.stdout:
		fOUT = sys.stdout
	else:
		sz_utils.make_dirs_if_necessary(args.out)
		fOUT = open(args.out, 'w')
	ColorText().info("[poolseq_tk] outputting allele counts to table ...", "stderr")
	for pos in sorted(ac.iterkeys()):
		fOUT.write("%s\t%d\t%s" %(ac[pos][0], pos, "\t".join(ac[pos][1:3])))
		i = 3
		if len(args.pileups) != 1:
			while i <= len(ac[pos])-4:
				fOUT.write("\t%s" %(":".join(ac[pos][i:i+4])))
				i += 4
		else:			# case where only one pileup file provided
			while i <= len(ac[pos])-2:
				fOUT.write("\t%s" %(":".join(ac[pos][i:i+2])))
				i += 2
		fOUT.write("\n")
	ColorText().info(" [done]\n", "stderr")
	fOUT.close()
