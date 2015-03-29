'''
	python poolseq_tk.py count

	Description: Count alleles at each SNP give the pileups
	Author: Simo V. Zhang

	Input: pileup file with reads bases converted to corresponding alleles

	Output: pielup file with allele counts
		   (1) chr
		   (2) pos
		   (3) ref base
		   (4) alt base
		   (5) allele counts in the order of ref and alt, separated by colon

'''

import collections
import sys
import os

import sz_utils
from colortext import ColorText

def run_count(args):
	''' Counting alleles at each SNP in the given pileup files '''

	dPos = {}
	if args.pos:
		ColorText().info("[poolseq_tk] reading SNPs positions:", "stderr")
		with open(args.pos, 'r') as fPOS:
			for line in fPOS:
				tmp_line = line.strip().split("\t")
				chr = tmp_line[0]
				pos = int(tmp_line[1])
				if (chr, pos) not in dPos:
					dPos[chr, pos] = 1
		ColorText().info(" %d\n" %(len(dPos)), "stderr")
	else:
		ColorText().info("[poolseq_tk] no SNP positions provided ... [skipped]\n", "stderr")

	ac = collections.defaultdict(tuple)
	fOUT = open(args.out, 'w')
	for pileup in args.pileups:
		sz_utils.check_if_files_exist(pileup)
		nsnps = 0
		ColorText().info("[poolseq_tk] counting alleles in %s:" %(os.path.basename(pileup)), "stderr")
		with open(pileup, 'r') as fMPILEUP:
			for line in fMPILEUP:
				nsnps += 1
				tmp_line = line.strip().split("\t")
				chr = tmp_line[0]
				pos = int(tmp_line[1])
				if (chr, pos) in dPos:
					ref_base = tmp_line[2]
					alt_base = tmp_line[3]
					nRefAlleles, nAltAlleles = 0, 0
					if len(tmp_line) == 5:
						nRefAlleles = tmp_line[-1].count(ref_base) + \
									  tmp_line[-1].count(ref_base.lower())
						nAltAlleles = tmp_line[-1].count(alt_base) + \
									  tmp_line[-1].count(alt_base.lower())
					if (chr, pos) not in ac:
						ac[chr, pos] = [ref_base, alt_base, str(nRefAlleles), str(nAltAlleles)]
					else:
						ac[chr, pos] += [str(nRefAlleles), str(nAltAlleles)]
		ColorText().info(" %d SNPs parsed\n" %(nsnps), "stderr")

	fOUT = None
	if args.out == sys.stdout:
		fOUT = sys.stdout
	else:
		sz_utils.make_dirs_if_necessary(args.out)
		fOUT = open(args.out, 'w')
	ColorText().info("[poolseq_tk] outputting allele counts to table ...", "stderr")
	for k in sorted(ac.iterkeys()):
		chr = k[0]
		pos = k[1]
		i = 2
		if len(ac[k][i:]) == 2*len(args.pileups):
			fOUT.write("%s\t%d\t%s" %(chr, pos, "\t".join(ac[k][0:2])))
			while i <= len(ac[k])-4:
				fOUT.write("\t%s" %(":".join(ac[k][i:i+4])))
				i += 4
			fOUT.write("\n")
	ColorText().info(" [done]\n", "stderr")
	fOUT.close()

def parseReadsBases(reads_bases, refBase, altBase):
	i = 0
	nRefAlleles, nAltAlleles = 0, 0
	nOtherAlleles = 0
	cov = 0
	while i < len(reads_bases):
		if reads_bases[i] == '.':
			nRefAlleles += 1
			i += 1
		elif reads_bases[i] == ',':
			nRefAlleles += 1
			i += 1
		elif reads_bases[i] == altBase:
			nAltAlleles += 1
			i += 1
		elif reads_bases[i] ==  altBase.lower():
			nAltAlleles += 1
			i += 1
		elif reads_bases[i] in ['+', '-', '*']:
			if reads_bases[i] == '*':
				i += 1
			else:
				len_indel = int(re.search(r'\d+', reads_bases[i+1:i+3]).group())
				i += len_indel + len(str(len_indel)) + 1
		elif reads_bases[i] == '^':
			i += 2
		elif reads_bases[i] in ['N', 'n', '$']:
			i += 1
		else:
			nOtherAlleles += 1
			i += 1
		cov += 1
	return cov, nRefAlleles, nAltAlleles, nOtherAlleles

#def run_count(args):
#	''' Counting alleles at each SNP in the given pileup files '''
#	ac = collections.defaultdict(tuple)
#	fOUT = open(args.out, 'w')
#	for pileup in args.pileups:
#		sz_utils.check_if_files_exist(pileup)
#		nsnps = 0
#		ColorText().info("[poolseq_tk] counting alleles in %s:" %(os.path.basename(pileup)), "stderr")
#		with open(pileup, 'r') as fMPILEUP:
#			for line in fMPILEUP:
#				nsnps += 1
#				tmp_line = line.strip().split("\t")
#				chr = tmp_line[0]
#				pos = int(tmp_line[1])
#				ref_base = tmp_line[2]
#				alt_base = tmp_line[3]
#				ref_count, alt_count = 0, 0
#				baseCounts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
#				if tmp_line[-1] != "N/A":
#					ref_count = tmp_line[-1].count(ref_base) + \
#								tmp_line[-1].count(ref_base.lower())
#					alt_count = tmp_line[-1].count(alt_base) + \
#								tmp_line[-1].count(alt_base.lower())
#				if (chr, pos) not in ac:
#					ac[chr, pos] = [ref_base, alt_base, str(ref_count), str(alt_count)]
#				else:
#					ac[chr, pos] += [str(ref_count), str(alt_count)]
#				for k in baseCounts.iterkeys():
#					baseCounts[k] = tmp_line[-1].count(k) + \
#									tmp_line[-1].count(k.lower())
#				fOUT.write("%s\t%s\t%d:%d:%d:%d:%d\n" %(chr, pos, baseCounts['A'], baseCounts['T'], baseCounts['C'], baseCounts['G'], baseCounts['N']))
#				if pos not in ac:
#					ac[pos] = [chr, ref_base, alt_base, str(ref_count), str(alt_count)]
#				else:
#					ac[pos] += [str(ref_count), str(alt_count)]
#		ColorText().info(" %d SNPs parsed\n" %(nsnps), "stderr")
#	ColorText().info(" [done]\n", "stderr")
#	fOUT.close()

#		if len(args.pileups) != 1:
#			# this deals with case where some pools have no data at this site
#			if len(ac[pos][i:]) >= 2*len(args.pileups):
#				fOUT.write("%s\t%d\t%s" %(ac[pos][0], pos, "\t".join(ac[pos][1:3])))
#				while i <= len(ac[pos])-4:
#					fOUT.write("\t%s" %(":".join(ac[pos][i:i+4])))
#					i += 4
#				fOUT.write("\n")
#		else:			# case where only one pileup file provided
#			fOUT.write("%s\t%d\t%s" %(ac[pos][0], pos, "\t".join(ac[pos][1:3])))
#			while i <= len(ac[pos])-2:
#				fOUT.write("\t%s" %(":".join(ac[pos][i:i+2])))
#				i += 2
#			fOUT.write("\n")
