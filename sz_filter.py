import os
import sys
import argparse

def getopts():
	usage = "Filter SNPs that are not satisfied specified conditions"
	filter_parser = argparse.ArgumentParser(description=usage)
	filter_parser.add_argument("-ac",
								metavar="FILE",
								dest="ac_file",
								help="allele counts file")
	filter_parser.add_argument("-o",
								metavar="FILE",
								dest="out",
								help="output file without filtered SNPs")
	filter_parser.add_argument("-min_ref_ac",
								metavar="INT",
								dest="min_ref_ac",
								type=int,
								default=5,
								help="minimum number of the ref allele (3rd column) per sample/pool")
	filter_parser.add_argument("-min_alt_ac",
								metavar="INT",
								dest="min_alt_ac",
								type=int,
								default=5,
								help="minimum number of the alt allele (4th column) per sample/pool")
	filter_parser.add_argument("-min_cov",
								metavar="INT",
								dest="min_cov",
								type=int,
								default=10,
								help="specify minimum coverage per site per sample/pool")
	return filter_parser.parse_args()

def run_filter():
	args = getopts()
	with open(args.out, 'w') as fOUT:
		with open(args.ac_file, 'r') as fAC:
			for line in fAC:
				tmp_line = line.strip().split("\t")
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


