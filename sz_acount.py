import argparse
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
		outdir = os.path.dirname(os.path.realpath(args.out))
		sz_utils.make_dirs_if_necessary(outdir)
		fOUT = open(args.out, 'w')
	ColorText().info("[poolseq_tk] outputting allele counts to table ...", "stderr")
	for pos in sorted(ac.iterkeys()):
		fOUT.write("%s\t%d\t%s" %(ac[pos][0], pos, "\t".join(ac[pos][1:3])))
		i = 3
		while i <= len(ac[pos])-4:
			fOUT.write("\t%s" %(":".join(ac[pos][i:i+4])))
			i += 4
		fOUT.write("\n")
	ColorText().info(" [done]\n", "stderr")
	fOUT.close()
