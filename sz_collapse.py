import os
import sys
import argparse
import collections

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
	snp_pos = get_SNPs(args.snps)

	# second, reading each of the pileup files
	chr1, m1_info = read_mpileup(args.m1, args.offset1)
	chr2, m2_info = read_mpileup(args.m2, args.offset2)

	ColorText().info("[poolseq_tk] %s: %d SNPs parsed\n" %(m1_base, len(m1_info)), "stderr")
	ColorText().info("[poolseq_tk] %s: %d SNPs parsed\n" %(m2_base, len(m2_info)), "stderr")

	fOUT = None
	if args.out != sys.stdout:
		outdir = os.path.dirname(os.path.realpath(args.out))
		sz_utils.make_dirs_if_necessary(outdir)
		fOUT = open(args.out, 'w')
	else:
		fOUT = args.out
	ColorText().info("[poolseq_tk]: collapsing mpileups %s and %s ..."
					 %(m1_base, m2_base), "stderr")
	for pos in sorted(snp_pos.iterkeys()):
		reads_bases_collapsed = ""
		if pos in m1_info and pos in m2_info:
			'''
				snp_pos[pos][0]: ref base of m1 pileup
				snp_pos[pos][1]: ref base of m2 pileup
				m1_info[pos][0]: ref base of m1 pileup
				m2_info[pos][1]: ref base of m2 pileup
			'''
			if snp_pos[pos][0] == m1_info[pos][0] and snp_pos[pos][1] == m2_info[pos][0]:
				reads_bases_collapsed = collapse_reads_bases(m1_info[pos][0],
															 m2_info[pos][0],
															 m1_info[pos][1],
															 reads_bases_collapsed)
				reads_bases_collapsed = collapse_reads_bases(m2_info[pos][0],
															 m1_info[pos][0],
															 m2_info[pos][1],
															 reads_bases_collapsed)
				if reads_bases_collapsed == "":
					fOUT.write("%s/%s\t%d\t%s\t%s\tN/A\n"
									%(chr1, chr2, pos,
									  snp_pos[pos][0], snp_pos[pos][1]))
				else:
					fOUT.write("%s/%s\t%d\t%s\t%s\t%s\n"
									%(chr1, chr2, pos,
									  snp_pos[pos][0], snp_pos[pos][1], reads_bases_collapsed))
			else:
				# this should bark if the same sites having different states
				ColorText().error("SNP position: %d %s %s\t\tMpileup position: %d %s %s\n"
								 %(pos, snp_pos[pos][0], snp_pos[pos][1],
								   pos, m1_info[pos][0], m2_info[pos][0]),
								   "stderr")
		# SNPS missed in both pileup file
		elif pos not in m1_info and pos not in m2_info:
			fOUT.write("%s/%s\t%d\t%s\t%s\t%s\n"
							%(chr1, chr2, pos,
							  snp_pos[pos][0], snp_pos[pos][1], "N/A"))
		# SNPs in m1 pileup file but not in m2
		elif pos in m1_info and pos not in m2_info:
			reads_bases_collapsed = collapse_reads_bases(m1_info[pos][0],
														 snp_pos[pos][1],
														 m1_info[pos][1],
														 reads_bases_collapsed)
			if reads_bases_collapsed == "":
				fOUT.write("%s/%s\t%d\t%s\t%s\tN/A\n"
								%(chr1, chr2, pos,
								  snp_pos[pos][0], snp_pos[pos][1]))
			else:
				fOUT.write("%s/%s\t%d\t%s\t%s\t%s\n"
								%(chr1, chr2, pos,
								  snp_pos[pos][0], snp_pos[pos][1], reads_bases_collapsed))
		# SNPs in m2 pileup file but not in m1
		elif pos not in m1_info and pos in m2_info:
			reads_bases_collapsed = collapse_reads_bases(m2_info[pos][0],
														 snp_pos[pos][0],
														 m2_info[pos][1],
														 reads_bases_collapsed)
			if reads_bases_collapsed == "":
				fOUT.write("%s/%s\t%d\t%s\t%s\tN/A\n"
								%(chr1, chr2, pos,
								  snp_pos[pos][0], snp_pos[pos][1]))
			else:
				fOUT.write("%s/%s\t%d\t%s\t%s\t%s\n"
								%(chr1, chr2, pos,
								  snp_pos[pos][0], snp_pos[pos][1], reads_bases_collapsed))
	ColorText().info(" [done]\n", "stderr")
	fOUT.close()

def get_SNPs(snps_file):
	''' read the SNP positions into a dictionary of tuple '''
	snp_pos = collections.defaultdict(tuple)
	sz_utils.check_if_files_exist(snps_file)
	with open(snps_file, 'r') as fSNP:
		for line in fSNP:
			if not line.startswith('#'):
				tmp_line = line.strip().split("\t")
				'''
					key: snp position in integer
					value: a tuple with two elements
						   1) one allele corresponding to the ref base of first pileup
						   2) second allele correspoding to the ref base of second pileup
				'''
				snp_pos[int(tmp_line[0])] = (tmp_line[1], tmp_line[2])
	return snp_pos

def read_mpileup(mpileup_file, offset):
	''' read certain columns in a pileup file into a dictionary of tuple '''
	ColorText().info("[poolseq_tk]: reading %s ..." %(mpileup_file), "stderr")
	mpileup_info = collections.defaultdict(tuple)
	chr = ""
	sz_utils.check_if_files_exist(mpileup_file)
	with open(mpileup_file, 'r') as fMPILEUP:
		for line in fMPILEUP:
			tmp_line = line.strip().split("\t")
			chr = tmp_line[0]
			'''
				key: SNP position in integer
				value: a tuple with two elements
					   1) ref base at the position
					   2) reads base covering that position
						  if the coverage at that position > 0
						  else N/A
			'''
			if int(tmp_line[3]) == 0:		# the forth column: coverage at a position
				mpileup_info[int(tmp_line[1])+offset] = (tmp_line[2].upper(), "N/A")
			else:
				mpileup_info[int(tmp_line[1])+offset] = (tmp_line[2].upper(), tmp_line[4])
	ColorText().info(" [done]\n", "stderr")
	return chr, mpileup_info

def collapse_reads_bases(ref_base, alt_base, reads_bases, reads_bases_collapsed):
	''' collapse reads bases of the two pileup files '''
	if reads_bases != "N/A":
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
				len_indel = int(reads_bases[i+1])
				i += len_indel + 1 + 1
			elif reads_bases[i] == '^':
				i += 2
			else:
				i += 1
	return reads_bases_collapsed
