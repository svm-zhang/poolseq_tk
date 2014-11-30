import sys
import os
import multiprocessing as mp
import argparse
import sz_utils

import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
import rpy2.rlike.container as rlc

def cmd_interface():
	usage = "Run Fisher's Exact Test at each SNP"
	fisher_parser = argparse.ArgumentParser(description=usage)
	fisher_parser.add_argument("-ac",
								metavar="FILE",
								dest="ac_file",
								help="allele counts for one pool")
	fisher_parser.add_argument("-outp",
								metavar="PREFIX",
								dest="outp",
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
								default="BH",
								choices=["holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"],
								help="specify the adjustment methods")
#	fisher_parser.add_argument("-direction", metavar="STR", dest="oddsr_direction", choices=["greater", "less"], help="specify whether odds ration greater, or less, than 1")

	return fisher_parser.parse_args()

def run_fisher():
	args = cmd_interface()
	''' run Fisher's Exact test '''
	szmake_dirs_if_necessary(os.path.dirname(args.outp))
	tables = count2table(args.ac_file)[0]

	task_q = mp.JoinableQueue()
	result_q = mp.Queue()
	_create_procs(args.nproc,task_q, result_q)
	sz_utils._assign_tables(tables, task_q, args.nproc)

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
					sz_utils.results_outputter(fFDR, pos, tables[pos][0], "\t".join(tables[pos][1:3]), tables[pos][3:], pvals[pos], padjust[i], odds_ratios[pos])
					if ((args.oddsr_direction == "less" and odds_ratios[pos] < 1) or
						 args.oddsr_direction == "greater" and odds_ratios[pos] > 1):
						sz_utils._results_outputter(fEXPECT, pos, tables[pos][0], "\t".join(tables[pos][1:3]), tables[pos][3:], pvals[pos], padjust[i], odds_ratios[pos])
				sz_utils._results_outputter(fALL, pos, tables[pos][0], "\t".join(tables[pos][1:3]), tables[pos][3:], pvals[pos], padjust[i], odds_ratios[pos])
		sys.stdout.write(" [done]\n")
		sys.stdout.write("[pool_gwas]: Program finishes successfully\n")

def _create_procs(nproc, task_q, result_q):
	''' initialize processes '''
	sys.stdout.write("[pool_gwas]: Initializing processes ...")
	for _ in range(nproc):
		p = mp.Process(target=_fisher_worker, args=(task_q, result_q))
		p.daemon = True
		p.start()
	sys.stdout.write(" [done]\n")

def _fisher_worker(task_q, result_q):
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
