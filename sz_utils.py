import os
import sys
import shutil

from colortext import ColorText

def make_dirs_if_necessary(*dirs):
	for dir in dirs:
		if not os.path.exists(dir):
			os.makedirs(dir)

def check_if_files_exist(*files):
	for file in files:
		if not os.path.exists(file):
			ColorText().error("\n[poolseq_tk] ERROR: cannot find file %s\n"
							 %(os.path.realpath(file)))
			sys.exit(1)

def calculate_pval_cutoff(pvals, fdr):
	ntests = len(pvals)
	pre_pval = 0.0
	for i, pval in enumerate(sorted(pvals)):
		if pval > (float(i+1)/ntests)*fdr :
			return pre_pval
		pre_pval = pval

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

def _assign_tables(tables, task_q, nproc):
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

def _results_outputter(fHANDLE, pos, chr, bases, tables, pval, corr_pval, odds_ratio):
	fHANDLE.write("%s\t%d\t%s" %(chr, pos, bases))
	i = 0
	while i <= len(tables) - 4:
		fHANDLE.write("\t%s" %(":".join(tables[i:i+4])))
		i += 4
	fHANDLE.write("\t%.8f\t%.8f\t%.8f\n" %(pval, corr_pval, odds_ratio))
	fHANDLE.flush()
