'''
	poolseq_tk.py fisher
	Description: run Fisher's Exact test at each SNP
	Author: Simo V. Zhang

	Input: filtered or raw allele counts file

	Output: outprefix.all
			(1) chr
			(2) pos
			(3) ref base
			(4) alt base
			(5) allele counts
			(6) raw p-values
			(7) adjusted p-values
			(8) odds ratios

			outprefix.fdr: a set of SNPs survived FDR
'''

import sys
import os
import multiprocessing as mp
import argparse
import collections

import sz_utils
from colortext import ColorText

import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
import rpy2.rlike.container as rlc

def run_fisher(args):
	''' run Fisher's Exact test '''
	sz_utils.make_dirs_if_necessary(args.outp)
	sz_utils.check_if_files_exist(args.ac_file)
	tables = sz_utils._count2table(args.ac_file)[0]

	task_q = mp.JoinableQueue()
	result_q = mp.Queue()
	_create_procs(args.nproc,task_q, result_q)
	sz_utils._assign_tables(tables, task_q, args.nproc)

	try:
		task_q.join()
	except KeyboardInterrupt:
		ColorText().info("[poolseq_tk]: Terminated unexpectedly by keyboard\n",
						 "stderr")
		sys.exit()
	else:
		pvals, odds_ratios = {}, {}
		while args.nproc:
			pvals_split, odds_ratios_split = result_q.get()
			pvals.update(pvals_split)
			odds_ratios.update(odds_ratios_split)
			args.nproc -= 1
		ColorText().info("[poolseq_tk]: Running Fisher's Exact tests successfully\n", "stderr")

		# correcting raw p-values and make QQ plots
		ColorText().info("[poolseq_tk]: multi-testing correction using %s method at %d%% level ..."
						 %(args.adj_method, args.adj_cutoff*100), "stderr")
		raw_pvals = [pvals[k] for k in sorted(pvals.iterkeys())]
		raw_pvals_vector = robjects.FloatVector(raw_pvals)
		padjust = robjects.r['p.adjust'](raw_pvals_vector, method=args.adj_method)
		ColorText().info(" [done]\n", "stderr")
		ColorText().info("[poolseq_tk]: p-value cutoff using Benjamini.Hochberg procedure %.5e"
						 %(__bh_pvals_cutoff(pvals, args.adj_cutoff)), "stderr")
		ColorText().info(" [done]\n", "stderr")

		# output p-values
		ColorText().info("[poolseq_tk]: output to files ...", "stderr")
		out_all = args.outp + ".fisher.all"
		out_fdr = args.outp + ".fisher.fdr%d" %(args.adj_cutoff*100)
#		out_expect = args.outp + ".fisher.fdr%d.expect" %(args.adj_cutoff*100)
		with open(out_all, 'w') as fALL, open(out_fdr, 'w') as fFDR:
			for i, pos in enumerate(sorted(pvals.iterkeys())):
				if padjust[i] <= args.adj_cutoff:
					sz_utils._results_outputter(fFDR, pos, tables[pos][0], "\t".join(tables[pos][1:3]), tables[pos][3:], pvals[pos], padjust[i], odds_ratios[pos])
#					if ((args.oddsr_direction == "less" and odds_ratios[pos] < 1) or
#						 args.oddsr_direction == "greater" and odds_ratios[pos] > 1):
#						sz_utils._results_outputter(fEXPECT, pos, tables[pos][0], "\t".join(tables[pos][1:3]), tables[pos][3:], pvals[pos], padjust[i], odds_ratios[pos])
				sz_utils._results_outputter(fALL, pos, tables[pos][0], "\t".join(tables[pos][1:3]), tables[pos][3:], pvals[pos], padjust[i], odds_ratios[pos])
		ColorText().info(" [done]\n", "stderr")
		ColorText().info("[poolseq_tk]: Program finishes successfully\n", "stderr")

def __bh_pvals_cutoff(dPvals, fdr):
	'''
		Using BH procedure to calculate pvalue
		cutoff at a FDR level
	'''
	lPvals = [dPvals[k] for k in dPvals.iterkeys()]
	ntests = len(lPvals)
	sort_lPvals = sorted(lPvals)
	for i in range(len(sort_lPvals)):
		if sort_lPvals[i] > (float(i+1)/ntests)*fdr :
			return sort_lPvals[i - 1]
	if i == len(sort_lPvals):
		ColorText().error("[poolseq_tk] Fail to calculate pvalue cutoff\n")
		sys.exit()

def _create_procs(nproc, task_q, result_q):
	''' initialize processes '''
	ColorText().info("[poolseq_tk]: Initializing processes ...\n", "stderr")
	for _ in range(nproc):
		p = mp.Process(target=_fisher_worker, args=(task_q, result_q))
		p.daemon = True
		p.start()

def _fisher_worker(task_q, result_q):
	while True:
		try:
			tables, nth_job = task_q.get()
			ColorText().info("[poolseq_tk]: %s running Fisher's Exact test on %d tables ...\n"
							 %(mp.current_process().name, len(tables)), "stderr")
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
			ColorText().info("[poolseq_tk]: %s ran %d tests\n"
							 %(mp.current_process().name, len(pvals_split)),
							 "stderr")
			result_q.put((pvals_split, odds_ratios_split))
		finally:
			task_q.task_done()
