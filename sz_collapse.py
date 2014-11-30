import os
import sys
import argparse

def getopts():
	usage = "Collapse two pileup files at corresponding SNPs"
	collapse_parser = argparse.ArgumentParser(description=usage)
	collapse_parser.add_argument("-m1",
								 metavar="FILE",
								 dest="m1",
								 help="one of the two mpileup files")
	collapse_parser.add_argument("-m2",
								 metavar="FILE",
								 dest="m2",
								 help="one of the two mpileup files")
	collapse_parser.add_argument("-offset1",
								 metavar="INT",
								 dest="offset1",
								 type=int,
								 default=0,
								 help="offset add in for the first mpileup file specified by -m1")
	collapse_parser.add_argument("-offset2",
								 metavar="INT",
								 dest="offset2",
								 type=int,
								 default=0,
								 help="offset add in for the second mpileup file specified by -m2")
	collapse_parser.add_argument("-snps",
								 metavar="FILE",
								 dest="snps",
								 help="a list of SNP positions. e.g. chr	pos")
	collapse_parser.add_argument("-p",
								 metavar="PREFIX",
								 dest="outp",
								 help="prefix for output file")
	return collapse_parser.parse_args()

def run_collapse():
	args = getopts()

	snp_pos = get_SNPs(args.snps)

	m1_info = read_mpileup(args.m1, args.offset1)
	m2_info = read_mpileup(args.m2, args.offset2)

	sys.stdout.write("[pool_gwas] %s: %d SNPs parsed\n" %(os.path.basename(args.m1), len(m1_info)))
	sys.stdout.write("[pool_gwas] %s: %d SNPs parsed\n" %(os.path.basename(args.m2), len(m2_info)))

	sys.stdout.write("[pool_gwas]: collapsing mpileups %s and %s ..." %(os.path.basename(args.m1), os.path.basename(args.m2)))
	merged_mpileup = args.outp + ".collapse.snps.mpileup"
	with open(merged_mpileup, 'w') as fCOLLAPSE:
		for pos in sorted(snp_pos.iterkeys()):
			cover_base_str = ""
			if pos in m1_info and pos in m2_info:
				if snp_pos[pos][0] == m1_info[pos][0] and snp_pos[pos][1] == m2_info[pos][0]:
					cover_base_str = collapse_cover_bases(m1_info[pos][0], m2_info[pos][0], m1_info[pos][1], cover_base_str)
					cover_base_str = collapse_cover_bases(m2_info[pos][0], m1_info[pos][0], m2_info[pos][1], cover_base_str)
					if cover_base_str == "":
						fCOLLAPSE.write("2L\t%d\t%s\t%s\tN/A\n" %(pos, snp_pos[pos][0], snp_pos[pos][1]))
					else:
						fCOLLAPSE.write("2L\t%d\t%s\t%s\t%s\n" %(pos, snp_pos[pos][0], snp_pos[pos][1], cover_base_str))
				else:
					# this should bark if the same sites having different states
					sys.stderr.write("SNP position: %d %s %s\t\tMpileup position: %d %s %s\n" %(pos, snp_pos[pos][0], snp_pos[pos][1], pos, m1_info[pos][0], m2_info[pos][0]))
			elif pos not in m1_info and pos not in m2_info:
					fCOLLAPSE.write("2L\t%d\t%s\t%s\t%s\n" %(pos, snp_pos[pos][0], snp_pos[pos][1], "N/A"))
			elif pos in m1_info and pos not in m2_info:
					cover_base_str = collapse_cover_bases(m1_info[pos][0], snp_pos[pos][1], m1_info[pos][1], cover_base_str)
					if cover_base_str == "":
						fCOLLAPSE.write("2L\t%d\t%s\t%s\tN/A\n" %(pos, snp_pos[pos][0], snp_pos[pos][1]))
					else:
						fCOLLAPSE.write("2L\t%d\t%s\t%s\t%s\n" %(pos, snp_pos[pos][0], snp_pos[pos][1], cover_base_str))
			elif pos not in m1_info and pos in m2_info:
					cover_base_str = collapse_cover_bases(m2_info[pos][0], snp_pos[pos][0], m2_info[pos][1], cover_base_str)
					if cover_base_str == "":
						fCOLLAPSE.write("2L\t%d\t%s\t%s\tN/A\n" %(pos, snp_pos[pos][0], snp_pos[pos][1]))
					else:
						fCOLLAPSE.write("2L\t%d\t%s\t%s\t%s\n" %(pos, snp_pos[pos][0], snp_pos[pos][1], cover_base_str))

	sys.stdout.write(" [done]\n")

def get_SNPs(snps_file):
	snp_pos = collections.defaultdict(tuple)
	with open(snps_file, 'r') as fSNP:
		for line in fSNP:
			if not line.startswith('#'):
				tmp_line = line.strip().split("\t")
				snp_pos[int(tmp_line[0])] = (tmp_line[1], tmp_line[2])
	return snp_pos

def read_mpileup(mpileup_file, offset):
	sys.stdout.write("[pool_gwas]: reading %s ..." %(mpileup_file))
	mpileup_info = collections.defaultdict(tuple)
	with open(mpileup_file, 'r') as fMPILEUP:
		for line in fMPILEUP:
			tmp_line = line.strip().split("\t")
			if int(tmp_line[3]) == 0:
				mpileup_info[int(tmp_line[1])+offset] = (tmp_line[2].upper(), "N/A")
			else:
				mpileup_info[int(tmp_line[1])+offset] = (tmp_line[2].upper(), tmp_line[4])
	sys.stdout.write(" [done]\n")
	return mpileup_info

def collapse_cover_bases(ref_base, alt_base, cover_bases, collapsed_cover_base_str):
	if cover_bases != "N/A":
		i = 0
		while i <= len(cover_bases)-1:
			if cover_bases[i] == '.':
				collapsed_cover_base_str += ref_base
				i += 1
			elif cover_bases[i] == ',':
				collapsed_cover_base_str += ref_base.lower()
				i += 1
			elif cover_bases[i] == alt_base:
				collapsed_cover_base_str += alt_base
				i += 1
			elif cover_bases[i] ==  alt_base.lower():
				collapsed_cover_base_str += alt_base.lower()
				i += 1
			elif cover_bases[i] in ['+', '-']:
				len_indel = int(cover_bases[i+1])
				i += len_indel + 1 + 1
			elif cover_bases[i] == '^':
				i += 2
			else:
				i += 1
	return collapsed_cover_base_str

