import os
import sys
import collections

import sz_utils

def run_biallelic(args):
	dPileups = syncPileups(args.pileups)

	sz_utils.make_dirs_if_necessary(args.out)
	fOUT = open(args.out, 'w')
	nZeroCov = 0
	nSNPsKept = 0
	nMulti = 0
	for k in sorted(dPileups.iterkeys()):
		chr = k[0]
		pos = k[1]
		ref_base = dPileups[k][0]
		alt_base = dPileups[k][1]
		reads_bases = dPileups[k][2]
		if len(reads_bases) > 0:
			ref_count = reads_bases.count(ref_base) + \
						reads_bases.count(ref_base.lower())
			alt_count = reads_bases.count(alt_base) + \
						reads_bases.count(alt_base.lower())
			other_count = len(reads_bases) - ref_count - alt_count
			if float(other_count)/len(reads_bases) < 0.05:
				fOUT.write("%s\t%d\t%s\t%s\n" %(chr, pos, ref_base, alt_base))
				nSNPsKept += 1
			else:
				nMulti += 1
		else:
			nZeroCov += 1
	fOUT.close()

	print nSNPsKept
	print nMulti
	print nZeroCov

def syncPileups(pileups):
	dPileups = collections.defaultdict(list)
	for pileup in pileups:
		with open(pileup, 'r') as fIN:
			for line in fIN:
				tmp_line = line.strip().split("\t")
				chr = tmp_line[0]
				pos = int(tmp_line[1])
				ref_base = tmp_line[2]
				alt_base = tmp_line[3]
				reads_bases = ""
				if len(tmp_line) == 5:
					reads_bases = tmp_line[4]
				else:
					reads_bases = ""
				if (chr, pos) not in dPileups:
					dPileups[chr, pos] = [ref_base, alt_base, reads_bases]
				else:
					dPileups[chr, pos][2] += reads_bases

	return dPileups

