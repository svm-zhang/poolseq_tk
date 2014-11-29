import argparse
import collections
import sys
import os

def parse_cmd():
	count_parser = argparse.ArgumentParser(description="counting alleles")
	count_parser.add_argument("-o", metavar="FILE", dest="out_count", help="output file of allele counts at each SNP")
	count_parser.add_argument("mpileups", metavar="MILEUP FILE", nargs='+', help="mpileup files")
	return count_parser.parse_args()

def count_alleles():
	args = parse_cmd()

	ac = collections.defaultdict(tuple)
	for mpileup in args.mpileups:
		nsnps = 0
		sys.stdout.write("[pool_gwas] counting alleles in %s:" %(os.path.basename(mpileup)))
		with open(mpileup, 'r') as fMPILEUP:
			for line in fMPILEUP:
				nsnps += 1
				tmp_line = line.strip().split("\t")
				pos = int(tmp_line[1])
				ref_base = tmp_line[2]
				alt_base = tmp_line[3]
				ref_count, alt_count = 0, 0
				if tmp_line[-1] != "N/A":
					ref_count = tmp_line[-1].count(ref_base) + tmp_line[-1].count(ref_base.lower())
					alt_count = tmp_line[-1].count(alt_base) + tmp_line[-1].count(alt_base.lower())
				if pos not in ac:
					ac[pos] = [tmp_line[0], ref_base, alt_base, str(ref_count), str(alt_count)]
				else:
					ac[pos] += [str(ref_count), str(alt_count)]
		sys.stdout.write(" %d SNPs parsed\n" %(nsnps))

	sys.stdout.write("[pool_gwas] outputting allele counts to table ...")
	with open(args.out_count, 'w') as fCOUNT:
		for pos in sorted(ac.iterkeys()):
			fCOUNT.write("%s\t%d\t%s" %(ac[pos][0], pos, "\t".join(ac[pos][1:3])))
			i = 3
			while i <= len(ac[pos])-4:
				fCOUNT.write("\t%s" %(":".join(ac[pos][i:i+4])))
				i += 4
			fCOUNT.write("\n")
	sys.stdout.write(" [done]\n")
