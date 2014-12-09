import os
import sys
import argparse
import collections

def getopts():
	usage = "Get overlaps of significant SNPs between replicates/pools"
	overlap_parser = argparse.ArgumentParser(description=usage)
	overlap_parser.add_argument("-a", metavar="FILE", dest="file_a", help="significant SNPs identified from pool A")
	overlap_parser.add_argument("-b", metavar="FILE", dest="file_b", help="significant SNPs identified from pool B")
	overlap_parser.add_argument("-o", metavar="FILE", dest="out", help="output file of overlapion of significant SNPs identified from both pools")
	return overlap_parser.parse_args()

def overlap(args):
	''' getting SNPs identified from both pools '''
	snp_a = collections.defaultdict(list)
	with open(args.file_a, 'r') as fA:
		for line in fA:
			tmp_line = line.strip().split("\t")
			snp_a[int(tmp_line[1])] = tmp_line
	sys.stdout.write("[pool_gwas]: %d SNPs parsed from %s\n" %(len(snp_a), os.path.basename(args.file_a)))

	num_overlapion = 0
	with open(args.out, 'w') as fOUT:
		with open(args.file_b, 'r') as fB:
			for line in fB:
				tmp_line = line.strip().split("\t")
				if int(tmp_line[1]) in snp_a:
					num_overlapion += 1
					fOUT.write("%s\t%s\n" %("\t".join(snp_a[int(tmp_line[1])]), "\t".join(tmp_line[-4:])))
	sys.stdout.write("[pool_gwas]: %d SNPs identified from both pools\n" %(num_overlapion))

