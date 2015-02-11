import os
import sys
import argparse
import collections
import multiprocessing as mp
import glob
import subprocess
import shlex
import re

import sz_collapse
import sz_acount
import sz_mergeAC
import sz_filter
import sz_fisher
import sz_cmh
import sz_plotting
import sz_overlap
import sz_prepVCF
import sz_view

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

def getopts():
	parser = argparse.ArgumentParser(description="Toolkits for Genome-wide Association Mapping using Pooled Sequencing")
	sub_parsers = parser.add_subparsers(title="Commands", metavar="", dest="command")

	usage = "Preparing VCF file from tests result file for snpEff"
	prepVCF_parser = sub_parsers.add_parser("vcf", help=usage)
	prepVCF_parser.add_argument("-i",
								metavar="FILE",
								dest="infile",
								required="True",
								help="test result file generated from poolseq_tk.py fisher or poolseq_tk.py cmh")
	prepVCF_parser.add_argument("-o",
								metavar="FILE",
								dest="out",
								type=argparse.FileType('w'),
								default=sys.stdout,
								help="output in VCF format. Default: STDOUT")
	prepVCF_parser.add_argument("-samples",
								metavar="LIST",
								dest="samples",
								default="table1,table2,table3,table4",
								help="a list of sample names separated by comma")
	prepVCF_parser.add_argument("-filter",
								metavar="EXPR",
								nargs='*',
								dest="filters",
								default=list(),
								help="a set of filters to apply. Only support INFO field ratio, e.g. ratio>1")
	prepVCF_parser.add_argument("-fst",
								metavar="FILE",
								dest="ifst",
								help="a file of Fst values")
	prepVCF_parser.set_defaults(func=sz_prepVCF.run_prepVCF)

	usage = "Viewing mpileup file (transforming 5th column in mpileup into human readable letters)"
	view_parser = sub_parsers.add_parser("view", help=usage)
	view_parser.add_argument("-mpileup",
							 metavar="FILE",
							 dest="ipileup",
							 required=True,
							 help="mpileup file")
	view_parser.add_argument("-snp",
							 metavar="FILE",
							 dest="isnp",
							 required=True,
							 help="tab-delimited snp file with four columns: chr, pos, ref, alt")
	view_parser.add_argument("-o",
							 metavar="FILE",
							 dest="out",
							 help="tab-delimited file with five columns: chr, pos, ref, alt, transformed reads bases ")
	view_parser.set_defaults(func=sz_view.viewPileup)

	# Collapsing two mpileup files
	usage = "Collapsing two pileup files at corresponding SNPs"
	collapse_parser = sub_parsers.add_parser("collapse", help=usage)
	collapse_parser.add_argument("-m1",
								 metavar="FILE",
								 dest="m1",
								 required="True",
								 help="one of the two mpileup files")
	collapse_parser.add_argument("-m2",
								 metavar="FILE",
								 dest="m2",
								 required="True",
								 help="one of the two mpileup files")
	collapse_parser.add_argument("-snps",
								 metavar="FILE",
								 dest="snps",
								 required="True",
								 help="a list of SNP positions. e.g. chr\\tpos")
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
	collapse_parser.add_argument("-o",
								 metavar="FILE",
								 dest="out",
								 default=sys.stdout,
								 help="output file. Default: STDOUT")
	collapse_parser.set_defaults(func=sz_collapse.run_collapse)

	usage = "Counting number of alleles given number of pileup files"
	count_parser = sub_parsers.add_parser("count", help=usage)
	count_parser.add_argument("-o",
							  metavar="FILE",
							  dest="out",
							  default=sys.stdout,
							  help="output file of allele counts at each SNP. Default: STDOUT")
	count_parser.add_argument("pileups",
							  metavar="PILEUP",
							  nargs='+',
							  help="pileup files")
	count_parser.set_defaults(func=sz_acount.run_count)

	usage = "Filter SNPs that are not satisfied specified conditions"
	filter_parser = sub_parsers.add_parser("filter", help=usage)
	filter_parser.add_argument("-ac",
								metavar="FILE",
								dest="ac_file",
								required=True,
								help="allele counts file")
	filter_parser.add_argument("-o",
								metavar="FILE",
								dest="out",
								default=sys.stdout,
								help="output file without filtered SNPs. Default: STDOUT")
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
	filter_parser.set_defaults(func=sz_filter.run_filter)

	usage = "Merging allele counts from multiple replicates"
	mergeAC_parser = sub_parsers.add_parser("mergeAC", help=usage)
	mergeAC_parser.add_argument("-o",
								metavar="FILE",
								dest="out",
								default=sys.stdout,
								help="output file of combined counts at each SNP across replicates")
	mergeAC_parser.add_argument("acs",
								metavar="ac_file",
								nargs='+',
								help="allele counts files")
	mergeAC_parser.set_defaults(func=sz_mergeAC.run_merge)

	usage = "Run Fisher's Exact Test at each SNP"
	fisher_parser = sub_parsers.add_parser("fisher", help=usage)
	fisher_parser.add_argument("-ac",
								metavar="FILE",
								dest="ac_file",
								help="allele counts for one pool")
	fisher_parser.add_argument("-outp",
								metavar="PREFIX",
								dest="outp",
								default="poolseq_tk.fisher",
								help="output file for Fisher's Exact tests")
	fisher_parser.add_argument("-t",
								metavar="INT",
								dest="nproc",
								type=int,
								default=1,
								help="Specify number of processes running simultaneously")
	fisher_parser.add_argument("-adj_cutoff",
								metavar="FLOAT",
								dest="adj_cutoff",
								type=float,
								default=0.05,
								help="specify the cutoff below which adjusted p-values will be considered as significant")
	fisher_parser.add_argument("-adj_method",
								metavar="STR",
								dest="adj_method",
								default="fdr",
#								choices=["holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"],
								help="specify the adjustment methods. Only BH procedure supported")
	fisher_parser.add_argument("-direction",
							   metavar="STR",
							   dest="oddsr_direction",
							   choices=["greater", "less"],
							   required=True,
							   help="specify whether odds ration greater, or less, than 1")
	fisher_parser.set_defaults(func=sz_fisher.run_fisher)

	usage="run Cochran-Mantel-Haenszel test with multi-testing adjustment"
	cmh_parser = sub_parsers.add_parser("cmh", help=usage)
	cmh_parser.add_argument("-table",
							metavar="FILE",
							dest="table_file",
							required=True,
							help="output file with the table that CMH test run on")
	cmh_parser.add_argument("-outp",
							metavar="PREFIX",
							dest="outp",
							default="poolseq_tk.cmh",
							required=True, help="output file with CMH test results")
	cmh_parser.add_argument("-t",
							metavar="INT",
							dest="nproc",
							type=int,
							default=1,
							help="Specify number of processes running simultaneously")
	cmh_parser.add_argument("-adj_cutoff",
							metavar="FLOAT",
							dest="adj_cutoff",
							type=float,
							default=0.05,
							help="specify the cutoff below which adjusted p-values will be considered as significant")
	cmh_parser.add_argument("-adj_method",
							metavar="STR",
							dest="adj_method",
							default="BH",
							choices=["holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"],
							help="specify the adjustment methods")
	cmh_parser.add_argument("-direction",
							metavar="STR",
							dest="oddsr_direction",
							choices=["greater", "less"],
							   required=True,
							help="specify whether odds ration greater, or less, than 1")
	cmh_parser.set_defaults(func=sz_cmh.run_cmh)

	usage = "Making Manhattan plot and QQ plot"
	plot_parser = sub_parsers.add_parser("plot", help=usage)
	plot_parser.add_argument("-i",
							 metavar="FILE",
							 dest="input",
							 required=True,
							 help="input file of test results with all SNPs")
	plot_parser.add_argument("-highlight",
							 metavar="FILE",
							 dest="highlight_snps",
							 help="file of a list of SNPs to be highlighted in the Manhattan plot")
	plot_parser.add_argument("-outp",
							 metavar="PREFIX",
							 dest="outp",
							 help="prefix of output file")
	plot_parser.add_argument("-pcutoff",
							 metavar="FLOAT",
							 dest="pcutoff",
							 type=float,
							 help="specify the p value cutoff to draw on the Mahhatan plot")
	plot_parser.add_argument("-fdrlevel",
							  metavar="FLOAT",
							  dest="fdrlevel",
							  type=float,
							  default=0.05,
							  help="specify at which level FDR will be applied")
	plot_parser.add_argument("-qqtitle",
							  metavar="STR",
							  dest="qqtitle",
							  help="specify the title for QQ plot")
	plot_parser.add_argument("-manx",
							  metavar="STR",
							  dest="manx",
							  help="specify xlab for manhattan plot")
	plot_parser.add_argument("-manxlim",
							  metavar="STR",
							  dest="manxlim",
							  default="-",
							  help="an interval defined by min and max, sperated by comma, e.g. 19,27. Default=\"-\"")
	plot_parser.add_argument("-mantitle",
							  metavar="STR",
							  dest="mantitle",
							  help="specify the title for Manhattan plot")
	plot_mutual_group = plot_parser.add_mutually_exclusive_group(required=True)
	plot_mutual_group.add_argument("-pdf",
								   dest="pdf",
								   action="store_true",
								   help="output qqplot in pdf format")
	plot_mutual_group.add_argument("-png",
								   dest="png",
								   action="store_true",
								   help="output qqplot in pdf format. Probably not working!")
	plot_parser.set_defaults(func=sz_plotting.making_plot)

	adjust_parser = sub_parsers.add_parser("adjust", help="getting significant SNPs with FDR correction")
	diff_parser = sub_parsers.add_parser("diff", help="get SNPs that significant in one replicate but not in the other")

	usage = "Get overlaps of significant SNPs between replicates/pools"
	overlap_parser = sub_parsers.add_parser("overlap", help=usage)
	overlap_parser.add_argument("-a",
								metavar="FILE",
								dest="file_a",
								help="significant SNPs identified from pool A")
	overlap_parser.add_argument("-b",
								metavar="FILE",
								dest="file_b",
								help="significant SNPs identified from pool B")
	overlap_parser.add_argument("-o",
								metavar="FILE",
								dest="out",
								help="output file of overlapion of significant SNPs identified from both pools")
	overlap_parser.set_defaults(func=sz_overlap.run_overlap)

#	adjust_parser.set_defaults(func=multi_testing_correction)
#	diff_parser.set_defaults(func=call_diff)

	return parser.parse_args()

def main():
	args = getopts()
	args.func(args)

if __name__ == "__main__":
	main()
