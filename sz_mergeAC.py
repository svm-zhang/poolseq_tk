import os
import sys
import collections
import argparse

def getopts():
	usage = "Merging allele counts from multiple replicates"
	mergeAC_parser = argparse.ArgumentParser(description=usage)
	mergeAC_parser.add_argument("-o", metavar="FILE", dest="out", help="output file of combined counts at each SNP across replicates")
	mergeAC_parser.add_argument("acs", metavar="ac_file", nargs='+', help="allele counts files")
	return mergeAC_parser.parse_args()

def run_merge():
	''' combine allele counts across replicates '''
	args = getopts()
	allele_counts = collections.defaultdict(list)
	data = collections.defaultdict(list)
	for ac_file in args.acs:
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

	# output to file
	with open(args.out, 'w') as fOUT:
		for pos in sorted(allele_counts.iterkeys()):
			fOUT.write("%s\t%s\n" %("\t".join(data[pos]), ":".join(map(str, allele_counts[pos]))))


