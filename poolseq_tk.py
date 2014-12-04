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

	fisher_parser = sub_parsers.add_parser("fisher", help="run Fisher's Exact Test at each SNP")
	cmh_parser = sub_parsers.add_parser("cmh", help="run Cochran-Mantel-Haenszel test with multi-testing adjustment")
	plot_parser = sub_parsers.add_parser("plot", help="making Manhattan and Q-Q plots")
	adjust_parser = sub_parsers.add_parser("adjust", help="getting significant SNPs with FDR correction")
	overlap_parser = sub_parsers.add_parser("overlap", help="overlapion of significant SNPs identified from two pools")
	diff_parser = sub_parsers.add_parser("diff", help="get SNPs that significant in one replicate but not in the other")
	return parser.parse_args()

#	fisher_parser.set_defaults(func=run_fisher)
#	cmh_parser.set_defaults(func=run_cmh)
#	plot_parser.set_defaults(func=making_plot)
#	adjust_parser.set_defaults(func=multi_testing_correction)
#	overlap_parser.set_defaults(func=overlap)
#	diff_parser.set_defaults(func=call_diff)

def main():
	args = getopts()
	args.func(args)
	#	subcommand = sys.argv[1]
#	if subcommand == "collapse":
#		import sz_collapse
#		sz_collapse.run_collapse()
#	elif subcommand == "filter":
#		import sz_filter
#		sz_filter.run_filter()
#	elif subcommand == "mergeAC":
#		import sz_mergeAC
#		sz_mergeAC.run_merge()
#	elif subcommand == "count":
#		import sz_acount
#		sz_acount.count_alleles()
#	elif subcommand == "fisher":
#		import sz_fisher
#		sz_fisher.run_fisher()
#	elif subcommand == "cmh":
#		import sz_cmh
#		sz_cmh.run_cmh()
#	elif subcommand == "adjust":
#		import sz_multitests_corr
#		sz_multitests_corr.run_correction()
#	elif subcommand == "overlap":
#		import sz_overlap
#		sz_overlap.overlap()
#	elif subcommand == "diff":
#		import sz_diff
#		sz_diff.run_diff()
#	elif subcommand == "plot":
#		import sz_plotting
#		sz_plotting.making_plot()
#	else:
#		show_subcommand_overviews()

if __name__ == "__main__":
	main()
