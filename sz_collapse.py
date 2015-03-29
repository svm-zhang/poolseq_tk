import os
import sys
import argparse
import collections
import re

import sz_utils
from colortext import ColorText

def run_collapse(args):
	'''
		Given two pileup files of the same region, like 2l+ and 2la,
		collapse the pileups at each corresponding SNP
		Some SNPs are not reported in one or the other pileup file.
		A full list of SNP positions are required
	'''
	m1_base = os.path.basename(args.m1)
	m2_base = os.path.basename(args.m2)

	# first, getting the full list of SNPs
	dSNPs = get_SNPs(args.snps)

	# second, reading each of the pileup files
	chr1, dM1 = read_mpileup(args.m1, args.offset1)
	chr2, dM2 = read_mpileup(args.m2, args.offset2)

	ColorText().info("[poolseq_tk] %s: %d SNPs parsed\n" %(m1_base, len(dM1)), "stderr")
	ColorText().info("[poolseq_tk] %s: %d SNPs parsed\n" %(m2_base, len(dM2)), "stderr")

	fOUT = None
	if args.out != sys.stdout:
		outdir = os.path.dirname(os.path.realpath(args.out))
		sz_utils.make_dirs_if_necessary(outdir)
		fOUT = open(args.out, 'w')
	else:
		fOUT = args.out
	ColorText().info("[poolseq_tk]: collapsing mpileups %s and %s ..."
					 %(m1_base, m2_base), "stderr")
	for pos in sorted(dSNPs.iterkeys()):
		reads_bases_collapsed = ""
		if pos in dM1 and pos in dM2:
			'''
				dSNPs[pos][0]: ref base of m1 pileup
				dSNPs[pos][1]: ref base of m2 pileup
				dM1[pos][0]: ref base of m1 pileup
				dM2[pos][1]: ref base of m2 pileup
			'''
			if dSNPs[pos][0] == dM1[pos][0] and dSNPs[pos][1] == dM2[pos][0]:
				reads_bases_collapsed = parseReadsBases(dM1[pos][0],
															 dM2[pos][0],
															 dM1[pos][1])
				reads_bases_collapsed += parseReadsBases(dM2[pos][0],
															 dM1[pos][0],
															 dM2[pos][1])
				fOUT.write("%s/%s\t%d\t%s\t%s\t%s\n"
								%(chr1, chr2, pos,
								  dSNPs[pos][0], dSNPs[pos][1], reads_bases_collapsed))
			else:
				# this should bark if the same sites having different states
				ColorText().error("SNP position: %d %s %s\t\tMpileup position: %d %s %s\n"
								 %(pos, dSNPs[pos][0], dSNPs[pos][1],
								   pos, dM1[pos][0], dM2[pos][0]),
								   "stderr")
		# SNPS missed in both pileup file
		elif pos not in dM1 and pos not in dM2:
			fOUT.write("%s/%s\t%d\t%s\t%s\n"
							%(chr1, chr2, pos,
							  dSNPs[pos][0], dSNPs[pos][1]))
		# SNPs in m1 pileup file but not in m2
		elif pos in dM1 and pos not in dM2:
			reads_bases_collapsed = parseReadsBases(dM1[pos][0],
														 dSNPs[pos][1],
														 dM1[pos][1])
			fOUT.write("%s/%s\t%d\t%s\t%s\t%s\n"
							%(chr1, chr2, pos,
							  dSNPs[pos][0], dSNPs[pos][1], reads_bases_collapsed))
		# SNPs in m2 pileup file but not in m1
		elif pos not in dM1 and pos in dM2:
			reads_bases_collapsed = parseReadsBases(dM2[pos][0],
														 dSNPs[pos][0],
														 dM2[pos][1])
			fOUT.write("%s/%s\t%d\t%s\t%s\t%s\n"
							%(chr1, chr2, pos,
							  dSNPs[pos][0], dSNPs[pos][1], reads_bases_collapsed))
	ColorText().info(" [done]\n", "stderr")
	fOUT.close()

def get_SNPs(snps_file):
	''' read the SNP positions into a dictionary of tuple '''
	dSNPs = collections.defaultdict(tuple)
	sz_utils.check_if_files_exist(snps_file)
	with open(snps_file, 'r') as fSNP:
		for line in fSNP:
			if not line.startswith('#'):
				tmp_line = line.strip().split("\t")
				pos = int(tmp_line[0])
				ref_base = tmp_line[1]
				alt_base = tmp_line[2]
				if pos not in dSNPs:
					dSNPs[pos] = (ref_base, alt_base)
	return dSNPs

def read_mpileup(mpileup_file, offset):
	''' read certain columns in a pileup file into a dictionary of tuple '''
	ColorText().info("[poolseq_tk]: reading %s ..." %(mpileup_file), "stderr")
	dMpileups = collections.defaultdict(tuple)
	chr = ""
	sz_utils.check_if_files_exist(mpileup_file)
	with open(mpileup_file, 'r') as fMPILEUP:
		for line in fMPILEUP:
			tmp_line = line.strip().split("\t")
			chr = tmp_line[0]
			pos = int(tmp_line[1])
			cov = int(tmp_line[3])
			ref_base = tmp_line[2].upper()
			if cov > 0:
				reads_bases = tmp_line[4]
				dMpileups[pos+offset] = (ref_base, reads_bases)
	ColorText().info(" [done]\n", "stderr")
	return chr, dMpileups
#			if cov == 0:		# the forth column: coverage at a position
#				dMpileups[pos+offset] = (tmp_line[2].upper(), "N/A")
#			else:
#				dMpileups[pos+offset] = (tmp_line[2].upper(), tmp_line[4])

def parseReadsBases(ref_base, alt_base, reads_bases):
	''' collapse reads bases of the two pileup files '''
	reads_bases_collapsed = ""
	if reads_bases != "":
		i = 0
		while i <= len(reads_bases)-1:
			if reads_bases[i] == '.':
				reads_bases_collapsed += ref_base
				i += 1
			elif reads_bases[i] == ',':
				reads_bases_collapsed += ref_base.lower()
				i += 1
			elif reads_bases[i] == alt_base:
				reads_bases_collapsed += alt_base
				i += 1
			elif reads_bases[i] ==  alt_base.lower():
				reads_bases_collapsed += alt_base.lower()
				i += 1
			elif reads_bases[i] in ['+', '-']:
				len_indel = int(re.search(r'\d+', reads_bases[i+1:i+3]).group())
				i += len_indel + len(str(len_indel)) + 1
			elif reads_bases[i] in ['N', 'n', '$', '*']:
				i += 1
			elif reads_bases[i] == '^':
				i += 2
			else:
				reads_bases_collapsed += reads_bases[i]
				i += 1
	return reads_bases_collapsed
