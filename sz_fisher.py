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
import math

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
	create_procs(args.nproc, task_q, result_q, args.outp)
	sz_utils._assign_tables(tables, task_q, args.nproc)

	try:
		task_q.join()
	except KeyboardInterrupt:
		ColorText().info("[poolseq_tk]: Terminated unexpectedly by keyboard\n",
						 "stderr")
		sys.exit()
	else:
		pvals, odds_ratios, log10_pvals = {}, {}, {}
		while args.nproc:
			file = result_q.get()
			with open(file, 'r') as fIN:
				for line in fIN:
					tmp_line = line.strip().split("\t")
					chr = tmp_line[0]
					pos = int(tmp_line[1])
					pval = float(tmp_line[2])
					odds_ratio = float(tmp_line[3])
					log10_pval = tmp_line[4]
					if (chr, pos) not in pvals:
						pvals[chr, pos] = pval
					if (chr, pos) not in odds_ratios:
						odds_ratios[chr, pos] = odds_ratio
					if (chr, pos) not in log10_pvals:
						log10_pvals[chr, pos] = log10_pval
			os.remove(file)
#			pvals_split, odds_ratios_split = result_q.get()
#			pvals.update(pvals_split)
#			odds_ratios.update(odds_ratios_split)
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
						 %(sz_utils.getFDR_BH(pvals, args.adj_cutoff)), "stderr")
		ColorText().info(" [done]\n", "stderr")

		# output p-values
		ColorText().info("[poolseq_tk]: output to files ...", "stderr")
		out_all = args.outp + ".fisher.all"
		out_fdr = args.outp + ".fisher.fdr%d" %(args.adj_cutoff*100)
		out_expect = args.outp + ".fisher.fdr%d.expect" %(args.adj_cutoff*100)
		with open(out_all, 'w') as fALL, \
			 open(out_fdr, 'w') as fFDR, \
			 open(out_expect, 'w') as fEXPECT:
			for i, k in enumerate(sorted(pvals.iterkeys())):
				chr = k[0]
				pos = k[1]
				raw_pval = pvals[k]
				log_pval = log10_pvals[k]
				odds_ratio = odds_ratios[k]
				if padjust[i] <= args.adj_cutoff:
					sz_utils._results_outputter(fFDR, pos, chr, "\t".join(tables[k][1:3]), tables[k][3:], raw_pval, log_pval, padjust[i], odds_ratio)
					if ((args.oddsr_direction == "greater" and odds_ratios[k] > 1) or
						(args.oddsr_direction == "less" and odds_ratios[k] < 1)):
						sz_utils._results_outputter(fEXPECT, pos, chr, "\t".join(tables[k][1:3]), tables[k][3:], raw_pval, log_pval, padjust[i], odds_ratio)
				sz_utils._results_outputter(fALL, pos, chr, "\t".join(tables[k][1:3]), tables[k][3:], raw_pval, log_pval, padjust[i], odds_ratio)
		ColorText().info(" [done]\n", "stderr")
		ColorText().info("[poolseq_tk]: Program finishes successfully\n", "stderr")

def create_procs(nproc, task_q, result_q, outp):
	''' initialize processes '''
	ColorText().info("[poolseq_tk]: Initializing processes ...\n", "stderr")
	for _ in range(nproc):
		p = mp.Process(target=fisher_worker, args=(task_q, result_q, outp))
		p.daemon = True
		p.start()

def fisher_worker(task_q, result_q, outp):
	while True:
		try:
			tables, nth_job = task_q.get()
			ColorText().info("[poolseq_tk]: %s running Fisher's Exact test on %d tables ...\n"
							 %(mp.current_process().name, len(tables)), "stderr")
			tmpFile = outp + "." + mp.current_process().name + ".fisher"
			fOUT = open(tmpFile, 'w')
			pvals_split, odds_ratios_split = {}, {}
			nTests = 0
			for k in sorted(tables.iterkeys()):
				oddsr = 0.0
				chr = k[0]
				pos = k[1]
				alt_base = tables[k][2]
				ref_base = tables[k][1]
				ref_ac1 = int(tables[k][3])
				alt_ac1 = int(tables[k][4])
				ref_ac2 = int(tables[k][5])
				alt_ac2 = int(tables[k][6])
				if (sum(map(int, tables[k][3:7])) >= 10 and
					alt_ac1 + ref_ac1 >= 5 and			# row subtotals
					alt_ac2 + ref_ac2 >= 5 and
					alt_ac1 + alt_ac2 >= 5 and			# column subtotals
					ref_ac1 + ref_ac2 >= 5):
					nTests += 1
					if (ref_ac1 == 0 or ref_ac2 == 0 or		# add pseudo counts in case
						alt_ac1 == 0 or alt_ac2 == 0):		# odds ratio goes to Inf
						ref_ac1 += 1
						ref_ac2 += 1
						alt_ac1 += 1
						alt_ac2 += 1
					data_vector = robjects.IntVector([ref_ac1, alt_ac1, ref_ac2, alt_ac2])
					table = robjects.r['matrix'](data_vector, ncol=2)
					rfisher = robjects.r['fisher.test'](table, alternative='t')
#					pvals_split[pos] = float(rfisher[0][0])
#					if (ref_ac1 == 0 or ref_ac2 == 0 or
#						alt_ac1 == 0 or alt_ac2 == 0):
#						oddsr = (float(ref_ac1+1)/(alt_ac1+1))/(float(ref_ac2+1)/(alt_ac2+1))
#					else:
					pvalue = float(rfisher[0][0])
					oddsr = rfisher[2][0]
#					odds_ratios_split[pos] = oddsr
					if pvalue == 0.0:
						fOUT.write("%s\t%d\t%.4g\t%.8f\tInf\n" %(chr, pos, pvalue, oddsr))
					elif pvalue == 1.0:
						fOUT.write("%s\t%d\t%.4g\t%.8f\t0.00000000\n" %(chr, pos, pvalue, oddsr))
					else:
						fOUT.write("%s\t%d\t%.8f\t%.8f\t%.8f\n" %(chr, pos, pvalue, oddsr, -1*math.log10(pvalue)))
			fOUT.close()
			ColorText().info("[poolseq_tk]: %s ran %d tests\n"
							 %(mp.current_process().name, nTests),
							 "stderr")
			result_q.put(tmpFile)
#			result_q.put((pvals_split, odds_ratios_split))
		finally:
			task_q.task_done()
