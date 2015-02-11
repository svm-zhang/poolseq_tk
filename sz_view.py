import os
import sys
import collections

from colortext import ColorText
from sz_utils import parseReadsBases, getSNPs

def viewPileup(ipileup, out, dSNPs):
	nRemoved = 0
	with open(ipileup, 'r') as fIN:
		for line in fIN:
			tmp_line = line.strip().split("\t")
			chr = tmp_line[0]
			pos = int(tmp_line[1])
			if (chr, pos) in dSNPs:
				cov = int(tmp_line[3])
				refBase = tmp_line[2].upper()
				altBase = dSNPs[chr, pos][1]
				if refBase == dSNPs[chr, pos][0]:
					if cov == 0:
						continue
					reads_bases = tmp_line[4]
					reads_bases_viewed, nReadsBases, nRefBases, dMultiBases, dIndels = parseReadsBases(reads_bases, refBase, altBase)
					out.write("%s\t%d\t%s\t%s\t%s\n" %(chr, pos, refBase, altBase,
													   reads_bases_viewed))

					# the following is a checkup on other alleles
					# at this moment this checkup is inactive
					# number of alleles (SNPs, Indels) other than the alternative allele
#					nMultiBases = sum(dMultiBases.values()) + sum(dIndels.values())
#					if (nReadsBases == nRefBases or
#						(nMultiBases)/float(nReadsBases) <= 0.05):
#						out.write("%s\t%d\t%s\t%s\t%s\n" %(chr, pos, refBase, altBase,
#														   reads_bases_viewed))
#					else:
#						nRemoved += 1
#						print pos, refBase, altBase, reads_bases, reads_bases_viewed
#						print dMultiBases
#						print dIndels
#						print
				else:
					sys.stderr.write("reference base not consistent\n")
					sys.stderr.write(line)
					sys.exit()
				del dSNPs[chr, pos]

	ColorText().info("[poolseq_tk]: There are %d sites where more than 2 alleles occur\n"
					 %(nRemoved), "stderr")
	ColorText().info("[poolseq_tk]: There are %d SNPs having no reads coverage for this dataset\n"
					 %(len(dSNPs)), "stderr")

def main():
	ipileup = sys.argv[1]
	isnps = sys.argv[2]
	out = sys.stdout

	dSNPs = getSNPs(isnps)
	viewPileup(ipileup, out, dSNPs)

if __name__ == "__main__":
	main()
