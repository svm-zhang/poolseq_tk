import os
import sys
import multiprocessing as mp
import argparse
import sz_utils

import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
import rpy2.rlike.container as rlc

def _getopts():
	usage="run Cochran-Mantel-Haenszel test with multi-testing adjustment"
	cmh_parser = argparse.ArgumentParser(description=usage)
	cmh_parser.add_argument("-table",
							metavar="FILE",
							dest="table_file",
							required=True,
							help="output file with the table that CMH test run on")
	cmh_parser.add_argument("-outp",
							metavar="PREFIX",
							dest="outp",
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
#	cmh_parser.add_argument("-direction", metavar="STR", dest="oddsr_direction", choices=["greater", "less"], help="specify whether odds ration greater, or less, than 1")
	return cmh_parser.parse_args()

def run_cmh():
	''' run Cochran-Mantel-Hasenzle test '''

	args = _getopts()
	sz_utils.make_dirs_if_necessary(os.path.dirname(args.outp))
	allele_counts = {}
	pvals = {}
	tables = collections.defaultdict(list)
	ntests = 0
	tables, ntables_per_snp = count2table(args.table_file)
	sys.stdout.write("[pool_gwas] %d tables prepared\n" %(len(tables)))

	task_q = mp.JoinableQueue()
	result_q = mp.Queue()
	_create_procs(args.nproc,task_q, result_q, ntables_per_snp)
	sz_utils._assign_tables(tables, task_q, args.nproc)

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
					sz_utils._results_outputter(fFDR, pos, tables[pos][0], "\t".join(tables[pos][1:3]), tables[pos][3:], pvals[pos], padjust[i], odds_ratios[pos])
					if ((args.oddsr_direction == "greater" and odds_ratios[pos] > 1) or
						(args.oddsr_direction == "less" and odds_ratios[pos] < 1)):
						sz_utils._results_outputter(fEXPECT, pos, tables[pos][0], "\t".join(tables[pos][1:3]), tables[pos][3:], pvals[pos], padjust[i], odds_ratios[pos])
				sz_utils._results_outputter(fALL, pos, tables[pos][0], "\t".join(tables[pos][1:3]), tables[pos][3:], pvals[pos], padjust[i], odds_ratios[pos])

def _cmh_worker(task_q, result_q, ntables_per_snp):
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

def _create_procs(nproc, task_q, result_q, ntables_per_snp):
	''' initialize processes '''
	sys.stdout.write("[pool_gwas] Initializing processes ...")
	for _ in range(nproc):
		p = mp.Process(target=cmh_worker, args=(task_q, result_q, ntables_per_snp))
		p.daemon = True
		p.start()
	sys.stdout.write(" [done]\n")
