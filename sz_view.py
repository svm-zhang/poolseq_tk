import os
import sys
import collections

from colortext import ColorText
from sz_utils import parseReadsBases
from sz_utils import check_if_files_exist, make_dirs_if_necessary

def run_view(args):
	check_if_files_exist(args.ipileup)
	make_dirs_if_necessary(args.out)

	dSNPs = getSNPs(args.isnp)
	fOUT = open(args.out, 'w')
#	nRemoved = 0
	with open(args.ipileup, 'r') as fIN:
		for line in fIN:
			tmp_line = line.strip().split("\t")
			chr = tmp_line[0]
			pos = int(tmp_line[1])
			if (chr, pos) in dSNPs:
				cov = int(tmp_line[3])
				ref_base = tmp_line[2].upper()
				alt_base = dSNPs[chr, pos][1]
				if ref_base == dSNPs[chr, pos][0]:
					if cov > 0:
						reads_bases = tmp_line[4]
						reads_bases_parsed = parseReadsBases(ref_base, alt_base, reads_bases)
						fOUT.write("%s\t%d\t%s\t%s\t%s\n" %(chr, pos, ref_base, alt_base,
														    reads_bases_parsed))
#						reads_bases_parsed, nReadsBases, nRefBases, dMultiBases, dIndels = parseReadsBases(reads_bases, ref_base, alt_base)

					# the following is a checkup on other alleles
					# at this moment this checkup is inactive
					# number of alleles (SNPs, Indels) other than the alternative allele
#					nMultiBases = sum(dMultiBases.values()) + sum(dIndels.values())
#					if (nReadsBases == nRefBases or
#						(nMultiBases)/float(nReadsBases) <= 0.05):
#						out.write("%s\t%d\t%s\t%s\t%s\n" %(chr, pos, ref_base, alt_base,
#														   reads_bases_parsed))
#					else:
#						nRemoved += 1
#						print pos, ref_base, alt_base, reads_bases, reads_bases_parsed
#						print dMultiBases
#						print dIndels
#						print
				else:
					sys.stderr.write("reference base not consistent\n")
					sys.stderr.write(line)
					sys.exit()
				del dSNPs[chr, pos]
	fOUT.close()

	ColorText().info("[poolseq_tk]: There are %d sites where more than 2 alleles occur\n"
					 %(nRemoved), "stderr")
	ColorText().info("[poolseq_tk]: There are %d SNPs having no reads coverage for this dataset\n"
					 %(len(dSNPs)), "stderr")

def getSNPs(isnp):
	'''
		getting polymorphic sties from file
		file format:
					1. chr name
					2. pos
					3. ref allele
					4. alt allele
	'''
	check_if_files_exist(isnp)
	dSNPs = collections.defaultdict(tuple)
	with open(isnp, 'r') as fSNPS:
		for line in fSNPS:
			tmp_line = line.strip().split("\t")
			chr = tmp_line[0]
			pos = int(tmp_line[1])
			refBase = tmp_line[2]
			altBase = tmp_line[3]
			if not (chr, pos) in dSNPs:
				dSNPs[chr, pos] = (refBase, altBase)

	return dSNPs
