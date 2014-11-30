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

#	collapse_parser.set_defaults(func=collapse_mpileups)

	# Arguments for filtering unwanted SNPs
#	filter_parser.set_defaults(func=filter_SNPs)

	# counting alleles
#	count_parser.set_defaults(func=run_count)

	# arguments for combining allele counts from replicates
#	mergeAC_parser.set_defaults(func=mergeAC)

#	fisher_parser.set_defaults(func=run_fisher)

#	cmh_parser.set_defaults(func=run_cmh)

#	plot_parser.set_defaults(func=making_plot)

#	adjust_parser.set_defaults(func=multi_testing_correction)

#	overlap_parser.set_defaults(func=overlap)
#	diff_parser.set_defaults(func=call_diff)

	return parser.parse_args()



def cat_split_files(file_list, out_file):
	''' concatenate split files '''
	with open(out_file, 'w') as fOUT:
		for file in sorted(file_list):
			shutil.copyfileobj(open(file, 'r'), fOUT)

def main():
	subcommand = sys.argv[1]
	if subcommand == "collapse":
		import sz_collapse
		sz_collapse.run_collapse()
	elif subcommand == "filter":
		import sz_filter
		sz_filter.run_filter()
	elif subcommand == "mergeAC":
		import sz_mergeAC
		sz_mergeAC.run_merge()
	elif subcommand == "count":
		import sz_acount
		sz_acount.count_alleles()
	elif subcommand == "fisher":
		import sz_fisher
		sz_fisher.run_fisher()
	elif subcommand == "cmh":
		import sz_cmh
		sz_cmh.run_cmh()
	elif subcommand == "adjust":
		import sz_multitests_corr
		sz_multitests_corr.run_correction()
	elif subcommand == "overlap":
		import sz_overlap
		sz_overlap.overlap()
	elif subcommand == "diff":
		import sz_diff
		sz_diff.run_diff()
	elif subcommand == "plot":
		import sz_plotting
		sz_plotting.making_plot()
#	elif subcommand in ["-h", "--help"]:
	else:
		show_subcommand_overviews()

if __name__ == "__main__":
	main()
