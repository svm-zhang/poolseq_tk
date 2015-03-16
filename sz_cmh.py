'''
	python poolseq_tk.py cmh
	Description: Running Cochran-Mantel-Hasenzle test at each SNPs
	Author: Simo V. Zhang

	Input: allele counts file

	Output: outprefix.all
			(1) chr
			(2) pos
			(3) ref base
			(4) alt base
			(5) allele counts (more than two columns)
			(6) raw p-values
			(7) corrected p-values
			(8) odds ratio

			outprefix.fdr5: a subset of above SNPs survived FDR
'''

import os
import sys
import multiprocessing as mp
import argparse
import collections

import sz_utils
from colortext import ColorText

import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
import rpy2.rlike.container as rlc

def run_cmh(args):
	''' run Cochran-Mantel-Hasenzle test '''

	sz_utils.make_dirs_if_necessary(args.outp)
	allele_counts = {}
	pvals = {}
	tables = collections.defaultdict(list)
	ntests = 0
	tables, ntables_per_snp = sz_utils._count2table(args.table_file)
	ColorText().info("[poolseq_tk] %d tables prepared\n" %(len(tables)), "stderr")

	task_q = mp.JoinableQueue()
	result_q = mp.Queue()
	_create_procs(args.nproc,task_q, result_q, ntables_per_snp, args.outp)
	sz_utils._assign_tables(tables, task_q, args.nproc)

	# waiting for all tasks to be finished
	try:
		task_q.join()
	except KeyboardInterrupt:
		ColorText().info("[poolseq_tk]: Terminated unexpectedly by keyboard\n", "stderr")
		sys.exit()
	else:
		# merge results
		pvals, odds_ratios = {}, {}
		while args.nproc:
			file = result_q.get()
			with open(file, 'r') as fIN:
				for line in fIN:
					tmp_line = line.strip().split("\t")
					chr = tmp_line[0]
					pos = int(tmp_line[1])
					pval = float(tmp_line[2])
					odds_ratio = float(tmp_line[3])
					if (chr, pos) not in pvals:
						pvals[chr, pos] = pval
					if (chr, pos) not in odds_ratios:
						odds_ratios[chr, pos] = odds_ratio
			os.remove(file)
#			pvals_split, odds_ratios_split = result_q.get()
#			pvals.update(pvals_split)
#			odds_ratios.update(odds_ratios_split)
			args.nproc -= 1
		ColorText().info("[poolseq_tk]: Running CMH tests successfully\n", "stderr")

		# correcting raw p-values
		ColorText().info("[poolseq_tk]: multi-testing correction using %s method at %d%% level ..."
						 %(args.adj_method, args.adj_cutoff*100), "stderr")
		raw_pvals = [pvals[chr, pos] for chr, pos in sorted(pvals.iterkeys())]
		raw_pvals_vector = robjects.FloatVector(raw_pvals)
		padjust = robjects.r['p.adjust'](raw_pvals_vector, method=args.adj_method)
		ColorText().info(" [done]\n", "stderr")
		pcutoff = sz_utils.getFDR_BH(pvals, args.adj_cutoff)
		ColorText().info("[poolseq_tk]: p-value cutoff using Benjamini.Hochberg procedure %.5e"
						 %(pcutoff), "stderr")
		ColorText().info(" [done]\n", "stderr")

		# output p-values
		ColorText().info("[poolseq_tk]: output to files ...", "stderr")
		out_all = args.outp + ".cmh.all"
		out_fdr = args.outp + ".cmh.fdr%d" %(args.adj_cutoff*100)
		out_expect = args.outp + ".cmh.fdr%d.expect" %(args.adj_cutoff*100)
		sz_utils.make_dirs_if_necessary(out_all, out_fdr)
		with open(out_all, 'w') as fALL, \
			 open(out_fdr, 'w') as fFDR, \
			 open(out_expect, 'w') as fEXPECT:
			for i, k in enumerate(sorted(pvals.iterkeys())):
				chr = k[0]
				pos = k[1]
				raw_pval = pvals[chr, pos]
				odds_ratio = odds_ratios[chr, pos]
				if padjust[i] <= args.adj_cutoff:
					sz_utils._results_outputter(fFDR, pos, chr, "\t".join(tables[chr, pos][1:3]), tables[chr, pos][3:], pval, padjust[i], odds_ratio)
					if ((args.oddsr_direction == "greater" and odds_ratios[pos] > 1) or
						(args.oddsr_direction == "less" and odds_ratios[pos] < 1)):
						sz_utils._results_outputter(fEXPECT, pos, chr, "\t".join(tables[chr, pos][1:3]), tables[chr, pos][3:], pval, padjust[i], odds_ratio)
				sz_utils._results_outputter(fALL, pos, chr, "\t".join(tables[chr, pos][1:3]), tables[chr, pos][3:], pval, padjust[i], odds_ratio)
		ColorText().info(" [done]\n", "stderr")
		ColorText().info("[poolseq_tk]: Program finishes successfully\n", "stderr")

def _cmh_worker(task_q, result_q, ntables_per_snp, outp):
	while True:
		try:
			table_part, nth_job = task_q.get()
			pvals, odds_ratios = {}, {}
			ColorText().info("[poolseq_tk]: %s running Cochran-Mantel-Haenszel test on %d tables ...\n"
							 %(mp.current_process().name, len(table_part)),
							 "stderr")
			tmpFile = outp + "." + mp.current_process().name + ".cmh"
			fOUT = open(tmpFile, 'w')
			nTests = 0
			for chr, pos in sorted(table_part.iterkeys()):
				array = []
				i = 0
				while i <= len(table_part[chr, pos])-4:
					if(i > 2 and sum(map(int, table_part[chr, pos][i:i+4])) >= 10 and
					   int(table_part[chr, pos][i])+int(table_part[chr, pos][i+1]) >= 5 and
					   int(table_part[chr, pos][i])+int(table_part[chr, pos][i+2]) >= 5 and
					   int(table_part[chr, pos][i+2])+int(table_part[chr, pos][i+3]) >= 5 and
					   int(table_part[chr, pos][i+1])+int(table_part[chr, pos][i+3]) >= 5):
						array += map(int, table_part[chr, pos][i:i+4])
						i += 4
					else:
						i += 1
				if len(array) == ntables_per_snp*4:
					dim_vector = robjects.IntVector([2, 2, ntables_per_snp])
					data = robjects.r['array'](robjects.IntVector(array), dim=dim_vector)
					rcmh = robjects.r['mantelhaen.test'](data, alternative='t')
					nTests += 1
					fOUT.write("%s\t%d\t%.8f\t%.8f\n" %(chr, pos, rcmh[1][0], rcmh[3][0]))
#					pvals[pos] = float(rcmh[1][0])
#					odds_ratios[pos] = float(rcmh[3][0])
			fOUT.close()
			ColorText().info("[poolseq_tk]: %s ran %d tests\n"
							 %(mp.current_process().name, nTests),
							 "stderr")
			result_q.put(tmpFile)
#			result_q.put((pvals, odds_ratios))
		finally:
			task_q.task_done()

def _create_procs(nproc, task_q, result_q, ntables_per_snp, outp):
	''' initialize processes '''
	ColorText().info("[poolseq_tk] Initializing processes ...\n", "stderr")
	for _ in range(nproc):
		p = mp.Process(target=_cmh_worker, args=(task_q, result_q, ntables_per_snp, outp))
		p.daemon = True
		p.start()
