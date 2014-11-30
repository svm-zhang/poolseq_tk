import os
import sys
import argparse
import collections
import multiprocessing as mp
import glob
import subprocess
import shlex
import re
import shutil

import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
import rpy2.rlike.container as rlc

class SubcommandHelpFormatter(argparse.RawDescriptionHelpFormatter):
	def _format_action(self, action):
		flag = 0
		parts = super(argparse.RawDescriptionHelpFormatter, self)._format_action(action)
		if action.nargs == argparse.PARSER:
			sub_cmd = "\n"
			for i, part in enumerate(parts.split("\n")):
				if i == 0:
					continue
				else:
					if flag == 1:
						sub_cmd += 4*" "+ " ".join(filter(None, part.split(" "))) + "\n"
						flag = 0
						continue
					if len(part.split(" ")) > 4:
						if len(part.split(" ")[4]) > 7:
							sub_cmd += " ".join(part.split(" ")[0:5])
							flag = 1
						else:
							sub_cmd += " ".join(part.split(" ")[0:5]) + (9-len(part.split(" ")[4])+4)*" " + " ".join(filter(None, part.split(" "))) + "\n"
			return sub_cmd
		else:
			return parts

def show_subcommand_overviews():
	parser = argparse.ArgumentParser(description="Toolkits for Genome-wide Association Mapping using Pooled Sequencing")
	sub_parsers = parser.add_subparsers(title="Commands", metavar="", dest="command")

	# merging two mpileup files
	collapse_parser = sub_parsers.add_parser("collapse", help="collapse two mpileupis that covering the same set of SNPs")
	filter_parser = sub_parsers.add_parser("filter", help="filter out SNPs in allele counts file that do not satisfy given conditions")
	count_parser = sub_parsers.add_parser("count", help="counting alleles given mpileup file")
	mergeAC_parser = sub_parsers.add_parser("mergeAC", help="combining allele counts from replicates")
	fisher_parser = sub_parsers.add_parser("fisher", help="run Fisher's Exact Test at each SNP")
	cmh_parser = sub_parsers.add_parser("cmh", help="run Cochran-Mantel-Haenszel test with multi-testing adjustment")
	plot_parser = sub_parsers.add_parser("plot", help="making Manhattan and Q-Q plots")
	adjust_parser = sub_parsers.add_parser("adjust", help="getting significant SNPs with FDR correction")
	overlap_parser = sub_parsers.add_parser("overlap", help="overlapion of significant SNPs identified from two pools")
	diff_parser = sub_parsers.add_parser("diff", help="get SNPs that significant in one replicate but not in the other")


#	collapse_parser.add_argument("-m1", metavar="FILE", dest="m1", help="one of the two mpileup files")
#	collapse_parser.add_argument("-m2", metavar="FILE", dest="m2", help="one of the two mpileup files")
#	collapse_parser.add_argument("-offset1", metavar="INT", dest="offset1", type=int, default=0, help="offset add in for the first mpileup file specified by -m1")
#	collapse_parser.add_argument("-offset2", metavar="INT", dest="offset2", type=int, default=0, help="offset add in for the second mpileup file specified by -m2")
#	collapse_parser.add_argument("-snps", metavar="FILE", dest="snps", help="a list of SNP positions. e.g. chr	pos")
#	collapse_parser.add_argument("-p", metavar="PREFIX", dest="outp", help="prefix for output file")
#	collapse_parser.set_defaults(func=collapse_mpileups)

	# Arguments for filtering unwanted SNPs
#	filter_parser.add_argument("-ac", metavar="FILE", dest="ac_file", help="allele counts file")
#	filter_parser.add_argument("-o", metavar="FILE", dest="out", help="output file without filtered SNPs")
#	filter_parser.add_argument("-min_ref_ac", metavar="INT", dest="min_ref_ac", type=int, default=5, help="minimum number of the ref allele (3rd column) per sample/pool")
#	filter_parser.add_argument("-min_alt_ac", metavar="INT", dest="min_alt_ac", type=int, default=5, help="minimum number of the alt allele (4th column) per sample/pool")
#	filter_parser.add_argument("-min_cov", metavar="INT", dest="min_cov", type=int, default=10, help="specify minimum coverage per site per sample/pool")
#	filter_parser.set_defaults(func=filter_SNPs)

	# counting alleles
#	count_parser.add_argument("-i", metavar="FILE", dest="out_count", help="output file of allele counts at each SNP")
#	count_parser.add_argument("mpileups", metavar="MILEUP FILE", nargs='+', help="mpileup files")
#	count_parser.set_defaults(func=run_count)

	# arguments for combining allele counts from replicates
#	mergeAC_parser.add_argument("-o", metavar="FILE", dest="out", help="output file of combined counts at each SNP across replicates")
#	mergeAC_parser.add_argument("acs", metavar="ac_file", nargs='+', help="allele counts files")
#	mergeAC_parser.set_defaults(func=mergeAC)

	# arguments for running Fisher's Exact test
#	fisher_parser.add_argument("-ac", metavar="FILE", dest="ac_file", help="allele counts for one pool")
#	fisher_parser.add_argument("-outp", metavar="PREFIX", dest="outp", help="output file for Fisher's Exact tests")
#	fisher_parser.add_argument("-t", metavar="INT", dest="nproc", type=int, default=1, help="Specify number of processes running simultaneously")
#	fisher_parser.add_argument("-adj_cutoff", metavar="FLOAT", dest="adj_cutoff", type=float, default=0.05, help="specify the cutoff below which adjusted p-values will be considered as significant")
#	fisher_parser.add_argument("-adj_method", metavar="STR", dest="adj_method", default="BH", choices=["holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"], help="specify the adjustment methods")
#	fisher_parser.add_argument("-direction", metavar="STR", dest="oddsr_direction", choices=["greater", "less"], help="specify whether odds ration greater, or less, than 1")
#	fisher_parser.set_defaults(func=run_fisher)

	# arguments for running Cochran-Mantel-Haenszel test
#	cmh_parser.add_argument("-table", metavar="FILE", dest="table_file", required=True, help="output file with the table that CMH test run on")
#	cmh_parser.add_argument("-outp", metavar="PREFIX", dest="outp", required=True, help="output file with CMH test results")
#	cmh_parser.add_argument("-t", metavar="INT", dest="nproc", type=int, default=1, help="Specify number of processes running simultaneously")
#	cmh_parser.add_argument("-adj_cutoff", metavar="FLOAT", dest="adj_cutoff", type=float, default=0.05, help="specify the cutoff below which adjusted p-values will be considered as significant")
#	cmh_parser.add_argument("-adj_method", metavar="STR", dest="adj_method", default="BH", choices=["holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"], help="specify the adjustment methods")
#	cmh_parser.add_argument("-direction", metavar="STR", dest="oddsr_direction", choices=["greater", "less"], help="specify whether odds ration greater, or less, than 1")
#	cmh_parser.set_defaults(func=run_cmh)

	# arguments for making Q-Q plot and Manhattan plot
#	plot_parser.add_argument("-i", metavar="FILE", dest="input", required=True, help="input file of test results with all SNPs, e.g. *.fisher.all, *.cmh.all")
#	plot_parser.add_argument("-interests", metavar="FILE", dest="interests_snps", help="file of a list of SNPs of interests. These SNPs will be highlighted in the Manhattan plot")
#	plot_parser.add_argument("-outp", metavar="PREFIX", dest="outp", help="prefix of output file")
#	plot_parser.add_argument("-adj_cutoff", metavar="FLOAT", dest="adj_cutoff", type=float, default=0.05, help="specify the cutoff below which adjusted p-values will be considered as significant")
#	plot_mutual_group = plot_parser.add_mutually_exclusive_group(required=True)
#	plot_mutual_group.add_argument("-pdf", dest="pdf", action="store_true", help="output qqplot in pdf format")
#	plot_mutual_group.add_argument("-png", dest="png", action="store_true", help="output qqplot in pdf format")
#	plot_parser.set_defaults(func=making_plot)

	# arguments for adjusting p-values
#	adjust_parser.add_argument("-fisher", metavar="FILE", dest="power_file", help="*.fisher file with a set of SNPs at which statistical test ran with power")
#	adjust_parser.add_argument("-outp", metavar="PREFIX", dest="outp", help="prefix of output file")
#	adjust_parser.add_argument("-cutoff", metavar="FLOAT", dest="adj_cutoff", type=float, default=0.05, help="specify the FDR rate cutoff")
#	adjust_parser.set_defaults(func=multi_testing_correction)

	# arguments for get overlapion of significant SNPs between pools, or across replicates
#	overlap_parser.add_argument("-a", metavar="FILE", dest="file_a", help="significant SNPs identified from pool A")
#	overlap_parser.add_argument("-b", metavar="FILE", dest="file_b", help="significant SNPs identified from pool B")
#	overlap_parser.add_argument("-o", metavar="FILE", dest="out", help="output file of overlapion of significant SNPs identified from both pools")
#	overlap_parser.set_defaults(func=overlap)

	# arugments for get significant SNPs in one pool/replicate, but not in the other
#	diff_parser.add_argument("-sigrep1", metavar="FILE", dest="sigrep1_file", help="file with a set of significant SNPs identified from one replicated, e.g. *.b.fdr, *.b.fdr.expect")
#	diff_parser.add_argument("-rep2", metavar="FILE", dest="rep2_file", help="file with all SNPs identified from another replicate, e.g. *.a.power")
#	diff_parser.add_argument("-o", metavar="FILE", dest="out", help="output file with SNPs significant in -sigrep1 file but not in -rep2 file")
#	diff_parser.set_defaults(func=call_diff)

	return parser.parse_args()

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

def collapse_mpileups(args):
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

def filter_SNPs(args):
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

def mergeAC(args):
	''' combine allele counts across replicates '''
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

def cat_split_files(file_list, out_file):
	''' concatenate split files '''
	with open(out_file, 'w') as fOUT:
		for file in sorted(file_list):
			shutil.copyfileobj(open(file, 'r'), fOUT)

def main():
	subcommand = sys.argv[1]
	if subcommand == "count":
		import sz_acount
		sz_acount.count_alleles()
	elif subcommand == "fisher":
		import sz_fisher
		sz_fisher.run_fisher()
	elif subcommand == "cmh":
		import sz_cmh
		sz_cmh.run_cmh()
	elif subcommand in ["-h", "--help"]:
		show_subcommand_overviews()

if __name__ == "__main__":
	main()
