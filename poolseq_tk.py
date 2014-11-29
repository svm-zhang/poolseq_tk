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

def parse_cmd():
	parser = argparse.ArgumentParser(description="Association mapping using pooled samples")
	sub_parsers = parser.add_subparsers(title="Commands", metavar="[command]", dest="command")

	# merging two mpileup files
	collapse_parser = sub_parsers.add_parser("collapse", help="collapse two mpileupis that covering the same set of SNPs")
	collapse_parser.add_argument("-m1", metavar="FILE", dest="m1", help="one of the two mpileup files")
	collapse_parser.add_argument("-m2", metavar="FILE", dest="m2", help="one of the two mpileup files")
	collapse_parser.add_argument("-offset1", metavar="INT", dest="offset1", type=int, default=0, help="offset add in for the first mpileup file specified by -m1")
	collapse_parser.add_argument("-offset2", metavar="INT", dest="offset2", type=int, default=0, help="offset add in for the second mpileup file specified by -m2")
	collapse_parser.add_argument("-snps", metavar="FILE", dest="snps", help="a list of SNP positions. e.g. chr	pos")
	collapse_parser.add_argument("-p", metavar="PREFIX", dest="outp", help="prefix for output file")
	collapse_parser.set_defaults(func=collapse_mpileups)

	# Arguments for filtering unwanted SNPs
	filter_parser = sub_parsers.add_parser("filter", help="filter out SNPs in allele counts file that do not satisfy given conditions")
	filter_parser.add_argument("-ac", metavar="FILE", dest="ac_file", help="allele counts file")
	filter_parser.add_argument("-o", metavar="FILE", dest="out", help="output file without filtered SNPs")
	filter_parser.add_argument("-min_ref_ac", metavar="INT", dest="min_ref_ac", type=int, default=5, help="minimum number of the ref allele (3rd column) per sample/pool")
	filter_parser.add_argument("-min_alt_ac", metavar="INT", dest="min_alt_ac", type=int, default=5, help="minimum number of the alt allele (4th column) per sample/pool")
	filter_parser.add_argument("-min_cov", metavar="INT", dest="min_cov", type=int, default=10, help="specify minimum coverage per site per sample/pool")
	filter_parser.set_defaults(func=filter_SNPs)

	# counting alleles
	count_parser = sub_parsers.add_parser("count", help="counting alleles given mpileup file")
	count_parser.add_argument("-o", metavar="FILE", dest="out_count", help="output file of allele counts at each SNP")
	count_parser.add_argument("mpileups", metavar="MILEUP FILE", nargs='+', help="mpileup files")
	count_parser.set_defaults(func=count_alleles)

	# arguments for combining allele counts from replicates
	combineAC_parser = sub_parsers.add_parser("combineAC", help="combining allele counts from replicates")
	combineAC_parser.add_argument("-o", metavar="FILE", dest="out", help="output file of combined counts at each SNP across replicates")
	combineAC_parser.add_argument("acs", metavar="ac_file", nargs='+', help="allele counts files")
	combineAC_parser.set_defaults(func=combineAC)

	# arguments for running Fisher's Exact test
	fisher_parser = sub_parsers.add_parser("fisher", help="run Fisher's Exact Test at each SNP")
	fisher_parser.add_argument("-ac", metavar="FILE", dest="ac_file", help="allele counts for one pool")
	fisher_parser.add_argument("-outp", metavar="PREFIX", dest="outp", help="output file for Fisher's Exact tests")
	fisher_parser.add_argument("-t", metavar="INT", dest="nproc", type=int, default=1, help="Specify number of processes running simultaneously")
	fisher_parser.add_argument("-adj_cutoff", metavar="FLOAT", dest="adj_cutoff", type=float, default=0.05, help="specify the cutoff below which adjusted p-values will be considered as significant")
	fisher_parser.add_argument("-adj_method", metavar="STR", dest="adj_method", default="BH", choices=["holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"], help="specify the adjustment methods")
	fisher_parser.add_argument("-direction", metavar="STR", dest="oddsr_direction", choices=["greater", "less"], help="specify whether odds ration greater, or less, than 1")
	fisher_parser.set_defaults(func=run_fisher)

	# arguments for running Cochran-Mantel-Haenszel test
	cmh_parser = sub_parsers.add_parser("cmh", help="run Cochran-Mantel-Haenszel test with multi-testing adjustment")
	cmh_parser.add_argument("-table", metavar="FILE", dest="table_file", required=True, help="output file with the table that CMH test run on")
	cmh_parser.add_argument("-outp", metavar="PREFIX", dest="outp", required=True, help="output file with CMH test results")
	cmh_parser.add_argument("-t", metavar="INT", dest="nproc", type=int, default=1, help="Specify number of processes running simultaneously")
	cmh_parser.add_argument("-adj_cutoff", metavar="FLOAT", dest="adj_cutoff", type=float, default=0.05, help="specify the cutoff below which adjusted p-values will be considered as significant")
	cmh_parser.add_argument("-adj_method", metavar="STR", dest="adj_method", default="BH", choices=["holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"], help="specify the adjustment methods")
	cmh_parser.add_argument("-direction", metavar="STR", dest="oddsr_direction", choices=["greater", "less"], help="specify whether odds ration greater, or less, than 1")
	cmh_parser.set_defaults(func=run_cmh)

	# arguments for making Q-Q plot and Manhattan plot
	plot_parser = sub_parsers.add_parser("plot", help="making Manhattan and Q-Q plots")
	plot_parser.add_argument("-i", metavar="FILE", dest="input", required=True, help="input file of test results with all SNPs, e.g. *.fisher.all, *.cmh.all")
	plot_parser.add_argument("-interests", metavar="FILE", dest="interests_snps", help="file of a list of SNPs of interests. These SNPs will be highlighted in the Manhattan plot")
	plot_parser.add_argument("-outp", metavar="PREFIX", dest="outp", help="prefix of output file")
	plot_parser.add_argument("-adj_cutoff", metavar="FLOAT", dest="adj_cutoff", type=float, default=0.05, help="specify the cutoff below which adjusted p-values will be considered as significant")
	plot_mutual_group = plot_parser.add_mutually_exclusive_group(required=True)
	plot_mutual_group.add_argument("-pdf", dest="pdf", action="store_true", help="output qqplot in pdf format")
	plot_mutual_group.add_argument("-png", dest="png", action="store_true", help="output qqplot in pdf format")
	plot_parser.set_defaults(func=making_plot)

	# arguments for adjusting p-values
	adjust_parser = sub_parsers.add_parser("adjust", help="getting significant SNPs with FDR correction")
	adjust_parser.add_argument("-fisher", metavar="FILE", dest="power_file", help="*.fisher file with a set of SNPs at which statistical test ran with power")
	adjust_parser.add_argument("-outp", metavar="PREFIX", dest="outp", help="prefix of output file")
	adjust_parser.add_argument("-cutoff", metavar="FLOAT", dest="adj_cutoff", type=float, default=0.05, help="specify the FDR rate cutoff")
	adjust_parser.set_defaults(func=multi_testing_correction)

	# arguments for get intersection of significant SNPs between pools, or across replicates
	intersect_parser = sub_parsers.add_parser("intersect", help="intersection of significant SNPs identified from two pools")
	intersect_parser.add_argument("-a", metavar="FILE", dest="file_a", help="significant SNPs identified from pool A")
	intersect_parser.add_argument("-b", metavar="FILE", dest="file_b", help="significant SNPs identified from pool B")
	intersect_parser.add_argument("-o", metavar="FILE", dest="out", help="output file of intersection of significant SNPs identified from both pools")
	intersect_parser.set_defaults(func=intersect)

	# arugments for get significant SNPs in one pool/replicate, but not in the other
	diff_parser = sub_parsers.add_parser("diff", help="get SNPs that significant in one replicate but not in the other")
	diff_parser.add_argument("-sigrep1", metavar="FILE", dest="sigrep1_file", help="file with a set of significant SNPs identified from one replicated, e.g. *.b.fdr, *.b.fdr.expect")
	diff_parser.add_argument("-rep2", metavar="FILE", dest="rep2_file", help="file with all SNPs identified from another replicate, e.g. *.a.power")
	diff_parser.add_argument("-o", metavar="FILE", dest="out", help="output file with SNPs significant in -sigrep1 file but not in -rep2 file")
	diff_parser.set_defaults(func=call_diff)

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

def count_alleles(args):
	ac = collections.defaultdict(tuple)
	for mpileup in args.mpileups:
		nsnps = 0
		sys.stdout.write("[pool_gwas] counting alleles in %s:" %(os.path.basename(mpileup)))
		with open(mpileup, 'r') as fMPILEUP:
			for line in fMPILEUP:
				nsnps += 1
				tmp_line = line.strip().split("\t")
				pos = int(tmp_line[1])
				ref_base = tmp_line[2]
				alt_base = tmp_line[3]
				ref_count, alt_count = 0, 0
				if tmp_line[-1] != "N/A":
					ref_count = tmp_line[-1].count(ref_base) + tmp_line[-1].count(ref_base.lower())
					alt_count = tmp_line[-1].count(alt_base) + tmp_line[-1].count(alt_base.lower())
				if pos not in ac:
					ac[pos] = [tmp_line[0], ref_base, alt_base, str(ref_count), str(alt_count)]
				else:
					ac[pos] += [str(ref_count), str(alt_count)]
		sys.stdout.write(" %d SNPs parsed\n" %(nsnps))

	sys.stdout.write("[pool_gwas] outputting allele counts to table ...")
	with open(args.out_count, 'w') as fCOUNT:
		for pos in sorted(ac.iterkeys()):
			fCOUNT.write("%s\t%d\t%s" %(ac[pos][0], pos, "\t".join(ac[pos][1:3])))
			i = 3
			while i <= len(ac[pos])-4:
				fCOUNT.write("\t%s" %(":".join(ac[pos][i:i+4])))
				i += 4
			fCOUNT.write("\n")
	sys.stdout.write(" [done]\n")

def combineAC(args):
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

def calculate_min_power(alpha):
	'''
		calculate the minimum row/column totals
		to make sure I have enough power
	'''
	power = 1
	while True:
		data_vector = robjects.IntVector([power, 0, 0, power])
		rfisher = robjects.r['fisher.test']
		test = rfisher(robjects.r['matrix'](data_vector, ncol=2), alternative='t')
		pval = test[0][0]
		if pval <= alpha:
			break
		else:
			power += 1
	return power

def run_fisher_worker(task_q, result_q):
	while True:
		try:
			tables, nth_job = task_q.get()
			sys.stdout.write("[pool_gwas]: %s running Fisher's Exact test on %d tables ...\n" %(mp.current_process().name, len(tables)))
			pvals_split, odds_ratios_split = {}, {}
			for pos in sorted(tables.iterkeys()):
				oddsr = 0.0
				alt_base = tables[pos][2]
				ref_base = tables[pos][1]
				ref_ac1 = int(tables[pos][3])
				alt_ac1 = int(tables[pos][4])
				ref_ac2 = int(tables[pos][5])
				alt_ac2 = int(tables[pos][6])
				if (sum(map(int, tables[pos][3:7])) >= 10 and
					alt_ac1 + ref_ac1 >= 5 and			# row subtotals
					alt_ac2 + ref_ac2 >= 5 and
					alt_ac1 + alt_ac2 >= 5 and			# column subtotals
					ref_ac1 + ref_ac2 >= 5):
					data_vector = robjects.IntVector([ref_ac1, alt_ac1, ref_ac2, alt_ac2])
					table = robjects.r['matrix'](data_vector, ncol=2)
					rfisher = robjects.r['fisher.test'](table, alternative='t')
					pvals_split[pos] = float(rfisher[0][0])
					if (ref_ac1 == 0 or ref_ac2 == 0 or
						alt_ac1 == 0 or alt_ac2 == 0):
						oddsr = (float(ref_ac1+1)/(alt_ac1+1))/(float(ref_ac2+1)/(alt_ac2+1))
					else:
						oddsr = rfisher[2][0]
					odds_ratios_split[pos] = oddsr
			result_q.put((pvals_split, odds_ratios_split))
		finally:
			task_q.task_done()

def cat_split_files(file_list, out_file):
	''' concatenate split files '''
	with open(out_file, 'w') as fOUT:
		for file in sorted(file_list):
			shutil.copyfileobj(open(file, 'r'), fOUT)

def count2table(ac_file):
	sys.stdout.write("[pool_gwas]: reading counts and preparing 2*2 tables ...")
	tables = collections.defaultdict(list)
	ntables_per_snp = 0
	with open(ac_file, 'r') as fAC:
		for line in fAC:
			tmp_line = line.strip().split("\t")
			if ntables_per_snp == 0:
				ntables_per_snp = len(tmp_line[4:])
			pos = int(tmp_line[1])
			base1 = tmp_line[2]
			base2 = tmp_line[3]
			tables[pos] = [tmp_line[0], base1, base2]		# chr, allele1, allele2
			for counts in tmp_line[4:]:
				tables[pos] += counts.split(':')			# counts
	sys.stdout.write(" [done]\n")
	return tables, ntables_per_snp

def creat_fisher_procs(nproc, task_q, result_q):
	''' initialize processes '''
	sys.stdout.write("[pool_gwas]: Initializing processes ...")
	for _ in range(nproc):
		p = mp.Process(target=run_fisher_worker, args=(task_q, result_q))
		p.daemon = True
		p.start()
	sys.stdout.write(" [done]\n")

def make_dirs_if_necessary(*dirs):
	for dir in dirs:
		if not os.path.exists(dir):
			os.makedirs(dir)

def run_fisher(args):
	''' run Fisher's Exact test '''
	make_dirs_if_necessary(os.path.dirname(args.outp))
	tables = count2table(args.ac_file)[0]

	task_q = mp.JoinableQueue()
	result_q = mp.Queue()
	creat_fisher_procs(args.nproc,task_q, result_q)
	assign_tables(tables, task_q, args.nproc)

	try:
		task_q.join()
	except KeyboardInterrupt:
		sys.stderr.write("[pool_gwas]: Terminated unexpectedly by keyboard\n")
		sys.exit()
	else:
		pvals, odds_ratios = {}, {}
		while args.nproc:
			pvals_split, odds_ratios_split = result_q.get()
			pvals.update(pvals_split)
			odds_ratios.update(odds_ratios_split)
			args.nproc -= 1

		# correcting raw p-values and make QQ plots
		sys.stdout.write("[pool_gwas]: multi-testing correction using %s method at %d%% level ..." %(args.adj_method, args.adj_cutoff*100))
		raw_pvals = [pvals[k] for k in sorted(pvals.iterkeys())]
		raw_pvals_vector = robjects.FloatVector(raw_pvals)
		padjust = robjects.r['p.adjust'](raw_pvals_vector, method=args.adj_method)
		sys.stdout.write(" [done]\n")

		# output p-values
		sys.stdout.write("[pool_gwas]: output to files ...")
		out_all = args.outp + ".snps.fisher.all"
		out_fdr = args.outp + ".snps.fisher.fdr%d" %(args.adj_cutoff*100)
		out_expect = args.outp + ".snps.fisher.fdr%d.expect" %(args.adj_cutoff*100)
		with open(out_all, 'w') as fALL, open(out_fdr, 'w') as fFDR, open(out_expect, 'w') as fEXPECT:
			for i, pos in enumerate(sorted(pvals.iterkeys())):
				if padjust[i] <= args.adj_cutoff:
					results_outputter(fFDR, pos, tables[pos][0], "\t".join(tables[pos][1:3]), tables[pos][3:], pvals[pos], padjust[i], odds_ratios[pos])
					if ((args.oddsr_direction == "less" and odds_ratios[pos] < 1) or
						 args.oddsr_direction == "greater" and odds_ratios[pos] > 1):
						results_outputter(fEXPECT, pos, tables[pos][0], "\t".join(tables[pos][1:3]), tables[pos][3:], pvals[pos], padjust[i], odds_ratios[pos])
				results_outputter(fALL, pos, tables[pos][0], "\t".join(tables[pos][1:3]), tables[pos][3:], pvals[pos], padjust[i], odds_ratios[pos])
		sys.stdout.write(" [done]\n")
		sys.stdout.write("[pool_gwas]: Program finishes successfully\n")

def making_plot(args):
	''' making Q-Q plot and Manhattan plot '''

	# install qqman package if not installed
	if not rpackages.isinstalled("qqman"):
		utils = rpackages.importr('utils')
		utils.chooseCRANmirror(ind=84)
		utils.install_packages("qqman")

	# get pvalues
	data = collections.defaultdict(tuple)
	pvals, adjust_pvals = {}, {}
	with open(args.input, 'r') as fIN:
		for line in fIN:
			tmp_line = line.strip().split("\t")
			data[int(tmp_line[1])] = ("_".join(tmp_line[0:2]), int(tmp_line[0][0]), int(tmp_line[1]))
			pvals[int(tmp_line[1])] = float(tmp_line[-3])
			adjust_pvals[int(tmp_line[1])] = float(tmp_line[-2])

	# get SNPs of interests
	snp_of_interests = []
	with open(args.interests_snps, 'r') as fINTEREST:
		for line in fINTEREST:
			snp_of_interests.append(line.strip())

	# getting p-value cutoff given FDR level
	unadjust_p_cutoff = calculate_pval_cutoff([pvals[k] for k in sorted(pvals.iterkeys())], args.adj_cutoff)
	if args.pdf:
		out_qqplot = args.outp + ".snps.fisher.qqplot.pdf"
		out_manhattan = args.outp + ".snps.fisher.manhattan.pdf"
	elif args.png:
		out_qqplot = args.outp + ".snps.fisher.qqplot.png"
		out_manhattan = args.outp + ".snps.fisher.manhattan.png"
	grdevices = rpackages.importr('grDevices')
	# making Q-Q plot
	raw_pvals_vector = robjects.FloatVector([pvals[k] for k in sorted(pvals.iterkeys())])
	adjust_pvals_vector = robjects.FloatVector([adjust_pvals[k] for k in sorted(adjust_pvals.iterkeys())])
	make_qqplots(grdevices, raw_pvals_vector, out_qqplot)

	# maing Manhattan plot
	make_manhattan(grdevices, data, raw_pvals_vector, snp_of_interests, out_manhattan)

def make_qqplots(grdevices, raw_pvals_vector, out_qqplot):
	''' making qqplot '''
	qqman = rpackages.importr('qqman')
	grdevices.pdf(out_qqplot)
#	robjects.r['par'](mfrow=robjects.IntVector([2,1]))
	qqman.qq(raw_pvals_vector, main="Q-Q plot for raw p-values using Fisher's Exact test")
	grdevices.dev_off()

def make_manhattan(grdevices, data, raw_pvals_vector, snp_of_interests, out_manhattan):
	snp_names = []
	snp_pos = []
	chr_names = []
	for pos in sorted(data.iterkeys()):
		snp_pos.append(pos)
		chr_names.append(int(data[pos][1]))
		snp_names.append(data[pos][0])
	od_raw = rlc.OrdDict([("SNP", robjects.StrVector(snp_names)),
					  ("CHR", robjects.IntVector(chr_names)),
					  ("BP", robjects.IntVector(snp_pos)),
					  ("P", robjects.FloatVector(raw_pvals_vector))])

	color_vector = robjects.StrVector(["blue4", "orange3"])
	sig_snps = robjects.StrVector(snp_of_interests)
	qqman = rpackages.importr('qqman')
	grdevices.pdf(out_manhattan)
#	robjects.r['par'](mfrow=robjects.IntVector([2,1]))
	qqman.manhattan(robjects.DataFrame(od_raw), highlight=sig_snps, col = color_vector, suggestiveline=3.716482, genomewideline=False, xlim=robjects.IntVector([20, 43]), ylim=robjects.IntVector([0,10]))
	grdevices.dev_off()

def assign_tables(tables, task_q, nproc):
	''' assigning each process a number of tables '''
	ntables = len(tables)
	ntables_per_proc = ntables/nproc + 1
	i = 0
	nth_job = 1
	while i + ntables_per_proc <= ntables:
		task_q.put((dict(sorted(tables.items())[i:i+ntables_per_proc+1]), nth_job))
		i += ntables_per_proc + 1
		nth_job += 1
	task_q.put((dict(sorted(tables.items())[i:]), nth_job))

def run_cmh_worker(task_q, result_q, ntables_per_snp):
	while True:
		try:
			table_part, nth_job = task_q.get()
			pvals, odds_ratios = {}, {}
			sys.stdout.write("[pool_gwas]: %s running Cochran-Mantel-Haenszel test on %d tables ...\n" %(mp.current_process().name, len(table_part)))
#			j = 0
			for pos in sorted(table_part.iterkeys()):
#				if j == 0:
#					print mp.current_process().name, pos, table_part[pos]
#					j += 1
				array = []
				i = 0
				while i <= len(table_part[pos])-4:
					if(i > 2 and sum(map(int, table_part[pos][i:i+4])) >= 8 and
					   int(table_part[pos][i])+int(table_part[pos][i+1]) >= 4 and
					   int(table_part[pos][i])+int(table_part[pos][i+2]) >= 4 and
					   int(table_part[pos][i+2])+int(table_part[pos][i+3]) >= 4 and
					   int(table_part[pos][i+1])+int(table_part[pos][i+3]) >= 4):
						array += map(int, table_part[pos][i:i+4])
						i += 4
					else:
						i += 1
				if len(array) == ntables_per_snp*4:
					dim_vector = robjects.IntVector([2, 2, ntables_per_snp])
					data = robjects.r['array'](robjects.IntVector(array), dim=dim_vector)
					rcmh = robjects.r['mantelhaen.test'](data, alternative='t', exact=True)
					pvals[pos] = float(rcmh[1][0])
					odds_ratios[pos] = float(rcmh[3][0])
			sys.stdout.write("[pool_gwas]: %s ran %d tests\n" %(mp.current_process().name, len(pvals)))
			result_q.put((pvals, odds_ratios))
		finally:
			task_q.task_done()

def create_cmh_procs(nproc, task_q, result_q, ntables_per_snp):
	''' initialize processes '''
	sys.stdout.write("[pool_gwas] Initializing processes ...")
	for _ in range(nproc):
		p = mp.Process(target=run_cmh_worker, args=(task_q, result_q, ntables_per_snp))
		p.daemon = True
		p.start()
	sys.stdout.write(" [done]\n")

def run_cmh(args):
	''' run Cochran-Mantel-Hasenzle test '''
	make_dirs_if_necessary(os.path.dirname(args.outp))
	allele_counts = {}
	pvals = {}
	tables = collections.defaultdict(list)
	ntests = 0
	tables, ntables_per_snp = count2table(args.table_file)
	sys.stdout.write("[pool_gwas] %d tables prepared\n" %(len(tables)))

	task_q = mp.JoinableQueue()
	result_q = mp.Queue()
	create_cmh_procs(args.nproc,task_q, result_q, ntables_per_snp)
	assign_tables(tables, task_q, args.nproc)

	# waiting for all tasks to be finished
	try:
		task_q.join()
	except KeyboardInterrupt:
		sys.stderr.write("[pool_gwas]: Terminated unexpectedly by keyboard\n")
		sys.exit()
	else:
		# merge results
		pvals, odds_ratios = {}, {}
		while args.nproc:
			pvals_split, odds_ratios_split = result_q.get()
			pvals.update(pvals_split)
			odds_ratios.update(odds_ratios_split)
			args.nproc -= 1

		# correcting raw p-values
		raw_pvals = [pvals[k] for k in sorted(pvals.iterkeys())]
		raw_pvals_vector = robjects.FloatVector(raw_pvals)
		padjust = robjects.r['p.adjust'](raw_pvals_vector, method=args.adj_method)

		# output p-values
		out_all = args.outp + ".snps.cmh.all"
		out_fdr = args.outp + ".snps.cmh.fdr%d" %(args.adj_cutoff*100)
		out_expect = args.outp + ".snps.cmh.fdr%d.expect" %(args.adj_cutoff*100)
		with open(out_all, 'w') as fALL, open(out_fdr, 'w') as fFDR, open(out_expect, 'w') as fEXPECT:
			for i, pos in enumerate(sorted(pvals.iterkeys())):
				if padjust[i] <= args.adj_cutoff:
					results_outputter(fFDR, pos, tables[pos][0], "\t".join(tables[pos][1:3]), tables[pos][3:], pvals[pos], padjust[i], odds_ratios[pos])
					if ((args.oddsr_direction == "greater" and odds_ratios[pos] > 1) or
						(args.oddsr_direction == "less" and odds_ratios[pos] < 1)):
						results_outputter(fEXPECT, pos, tables[pos][0], "\t".join(tables[pos][1:3]), tables[pos][3:], pvals[pos], padjust[i], odds_ratios[pos])
				results_outputter(fALL, pos, tables[pos][0], "\t".join(tables[pos][1:3]), tables[pos][3:], pvals[pos], padjust[i], odds_ratios[pos])

def results_outputter(fOUT, pos, chr, bases, tables, pval, corr_pval, odds_ratio):
	fOUT.write("%s\t%d\t%s" %(chr, pos, bases))
	i = 0
	while i <= len(tables) - 4:
		fOUT.write("\t%s" %(":".join(tables[i:i+4])))
		i += 4
	fOUT.write("\t%.8f\t%.8f\t%.8f\n" %(pval, corr_pval, odds_ratio))
	fOUT.flush()

def calculate_pval_cutoff(pvals, fdr):
	ntests = len(pvals)
	pre_pval = 0.0
	for i, pval in enumerate(sorted(pvals)):
		if pval > (float(i+1)/ntests)*fdr :
			return pre_pval
		pre_pval = pval

def multi_testing_correction(args):
	make_dirs_if_necessary(os.path.dirname(args.outp))
	info = {}
	pvals = []
	fFISHER = open(args.power_file, 'r')
	for line in fFISHER:
		tmp_line = re.split('\s+', line.strip())
		pvals.append(float(tmp_line[5]))

	p_cutoff = calculate_pval_cutoff(pvals, args.adj_cutoff)
	print p_cutoff

	out_sig = args.outp + ".fdr%d" %(args.adj_cutoff*100)
	out_expect = args.outp + ".fdr%d.expect" %(args.adj_cutoff*100)
	with open(out_sig, 'w') as fSIG, open(out_expect, 'w') as fEXPECT:
		table, header = "", ""
		fFISHER.seek(0)
		for line in fFISHER:
			tmp_line = re.split('\s+', line.strip())
			pvals.append(float(tmp_line[5]))
			pval = float(tmp_line[5])
			oddsr = float(tmp_line[6])
			if pval <= p_cutoff:
				fSIG.write("%s" %(line))
				if oddsr > 1:
					fEXPECT.write("%s" %(line))

def intersect(args):
	''' getting SNPs identified from both pools '''
	snp_a = collections.defaultdict(list)
	with open(args.file_a, 'r') as fA:
		for line in fA:
			tmp_line = line.strip().split("\t")
			snp_a[int(tmp_line[1])] = tmp_line
	sys.stdout.write("[pool_gwas]: %d SNPs parsed from %s\n" %(len(snp_a), os.path.basename(args.file_a)))

	num_intersection = 0
	with open(args.out, 'w') as fOUT:
		with open(args.file_b, 'r') as fB:
			for line in fB:
				tmp_line = line.strip().split("\t")
				if int(tmp_line[1]) in snp_a:
					num_intersection += 1
					fOUT.write("%s\t%s\n" %("\t".join(snp_a[int(tmp_line[1])]), "\t".join(tmp_line[-4:])))
	sys.stdout.write("[pool_gwas]: %d SNPs identified from both pools\n" %(num_intersection))

def diff(table1, table2, table_info, out):
	table1_base = os.path.basename(table1)
	table2_base = os.path.basename(table2)
	allele_counts = []
	with open(out, 'w') as fOUT:
		with open(table1, 'r') as fTABLE:
			for line in fTABLE:
				if line.startswith('>'):
					if allele_counts:
						if pos in table_info:
							tmp_header = header.split('|')
							tmp_table_info = table_info[pos][0].split('|')
#							if float(tmp_table_info[4]) >= 0.0012:
							fOUT.write(">%s:%s\t\t%s:%s\n" %(table1_base, header.lstrip('>'), table2_base, table_info[pos][0].lstrip('>')))
							fOUT.write("\t\tHigh\tLow\t\t\t\t\t\t\t\t\tHigh\tLow\n")
							fOUT.write("\t%s\t%s\t\t\t\t\t\t\t\t%s\t%s\n" %(tmp_header[2], "\t".join(allele_counts[:2]), tmp_header[2], "\t".join(table_info[pos][1][:2])))
							fOUT.write("\t%s\t%s\t\t\t\t\t\t\t\t%s\t%s\n" %(tmp_header[3], "\t".join(allele_counts[2:]), tmp_header[3], "\t".join(table_info[pos][1][2:])))
						allele_counts = []
					tmp_line = line.strip().split('|')
					header = line.strip()
					pos = tmp_line[1]
				else:
					tmp_line = line.strip().split("\t")
					if len(tmp_line) >= 3:
						allele_counts += [tmp_line[1], tmp_line[2]]
		if allele_counts:
			if pos in table_info:
				tmp_header = header.split('|')
				tmp_table_info = table_info[pos][0].split('|')
#				if int(tmp_table_info[4]) >= 0.0012:
				fOUT.write(">%s:%s\t\t%s:%s\n" %(table1_base, header.lstrip('>'), table2_base, table_info[pos][0].lstrip('>')))
				fOUT.write("\t\tHigh\tLow\t\t\t\t\t\t\t\t\tHigh\tLow\n")
				fOUT.write("\t%s\t%s\t\t\t\t\t\t\t\t%s\t%s\n" %(tmp_header[2], "\t".join(allele_counts[:2]), tmp_header[2], "\t".join(table_info[pos][1][:2])))
				fOUT.write("\t%s\t%s\t\t\t\t\t\t\t\t%s\t%s\n" %(tmp_header[3], "\t".join(allele_counts[2:]), tmp_header[3], "\t".join(table_info[pos][1][2:])))

def call_diff(args):
	table_info = collections.defaultdict(tuple)
	table_info = read_SNPs_from_pool(args.sigrep1_file, table_info)
	print len(table_info)
	diff(args.rep2_file, args.sigrep1_file, table_info, args.out)

def main():
	args = parse_cmd()
	args.func(args)

if __name__ == "__main__":
	main()
