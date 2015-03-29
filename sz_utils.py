import os
import sys
import shutil
import collections
import re

from colortext import ColorText

def make_dirs_if_necessary(*files):
	for file in files:
		dir = getdirs(file)
		if not os.path.exists(dir):
			os.makedirs(dir)

def getdirs(file):
	return os.path.dirname(os.path.realpath(file))

def check_if_files_exist(*files):
	for file in files:
		if not os.path.exists(file):
			ColorText().error("\n[poolseq_tk] ERROR: cannot find file %s\n"
							 %(os.path.realpath(file)))
			sys.exit(1)

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

def getFDR_BH(dPvals, fdr_level):
	'''
		Using BH procedure to calculate pvalue
		cutoff at a FDR level
	'''
	lPvals = [dPvals[k] for k in dPvals.iterkeys()]
	ntests = len(lPvals)
	sort_lPvals = sorted(lPvals)
	for i in xrange(len(sort_lPvals)):
		if sort_lPvals[i] > (float(i+1)/ntests)*fdr_level :
#			print i, ntests, (float(i+1)/ntests)*fdr_level, sort_lPvals[i], sort_lPvals[i-1]
			if i == 0:
				return 0.00000000
			else:
				return sort_lPvals[i - 1]
	if i == len(sort_lPvals):
		ColorText().error("[poolseq_tk] Fail to calculate pvalue cutoff\n")
		sys.exit()

def cat_split_files(file_list, out_file):
	''' concatenate split files '''
	with open(out_file, 'w') as fOUT:
		for file in sorted(file_list):
			shutil.copyfileobj(open(file, 'r'), fOUT)

def _count2table(ac_file, max_cov = 100):
	ColorText().info("[poolseq_tk]: reading counts and preparing 2*2 tables ...",
					 "stderr")
	tables = collections.defaultdict(list)
	ntables_per_snp = 0
	with open(ac_file, 'r') as fAC:
		for line in fAC:
			tmp_line = line.strip().split("\t")
			if ntables_per_snp == 0:
				ntables_per_snp = len(tmp_line[4:])
			chr = tmp_line[0]
			pos = int(tmp_line[1])
			base1 = tmp_line[2]
			base2 = tmp_line[3]
			tables[chr, pos] = [tmp_line[0], base1, base2]		# chr, allele1, allele2
			for counts in tmp_line[4:]:
				tmp_counts = counts.split(':')
				if sum(map(int, tmp_counts[0:2])) <= max_cov:
					tables[chr, pos] += tmp_counts[:2]
				if sum(map(int, tmp_counts[2:4])) <= max_cov:
					tables[chr, pos] += tmp_counts[2:]
#					tables[chr, pos] += counts.split(':')			# counts
			if len(tables[chr, pos]) < len(tmp_line[4:])*4+3:
				del tables[chr, pos]
	ColorText().info(" [done]\n", "stderr")
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

def _results_outputter(fHANDLE, pos, chr, bases, tables, pval, log10_pval, corr_pval, odds_ratio):
	fHANDLE.write("%s\t%d\t%s" %(chr, pos, bases))
	i = 0
	while i <= len(tables) - 4:
		fHANDLE.write("\t%s" %(":".join(tables[i:i+4])))
		i += 4
	fHANDLE.write("\t%.8f\t%s\t%.8f\t%.8f\n" %(pval, str(log10_pval), corr_pval, odds_ratio))
	fHANDLE.flush()

def parseReadsBases(ref_base, alt_base, reads_bases):
	''' collapse reads bases of the two pileup files '''
	reads_bases_parsed = ""
	if reads_bases != "":
		i = 0
		while i <= len(reads_bases)-1:
			if reads_bases[i] == '.':
				reads_bases_parsed += ref_base
				i += 1
			elif reads_bases[i] == ',':
				reads_bases_parsed += ref_base.lower()
				i += 1
			elif reads_bases[i] == alt_base:
				reads_bases_parsed += alt_base
				i += 1
			elif reads_bases[i] ==  alt_base.lower():
				reads_bases_parsed += alt_base.lower()
				i += 1
			elif reads_bases[i] in ['+', '-']:
				len_indel = int(re.search(r'\d+', reads_bases[i+1:i+3]).group())
				i += len_indel + len(str(len_indel)) + 1
			elif reads_bases[i] in ['N', 'n', '$', '*']:
				i += 1
			elif reads_bases[i] == '^':
				i += 2
			else:
				reads_bases_parsed += reads_bases[i]
				i += 1
	return reads_bases_parsed

def old_parseReadsBases(reads_bases, refBase, altBase):
	'''
		parsing the 5th column in a mpileup file
		5th column represents the bases on the reads
		covering the given site
	'''
	i = 0
	dIndels = {}
	dMultiBases = {}
	nRefBases, nAltBases = 0, 0
	nReadsBases = 0
	reads_bases_transformed = ""
	while i < len(reads_bases):
		# ref forward
		if reads_bases[i] == '.':
			reads_bases_transformed += refBase
			nRefBases += 1
			i += 1
		# ref reverse
		elif reads_bases[i] == ',':
			reads_bases_transformed += refBase.lower()
			nRefBases += 1
			i += 1
		# alt forward
		elif reads_bases[i] == altBase:
			reads_bases_transformed += altBase
			nAltBases += 1
			i += 1
		# alt reverse
		elif reads_bases[i] ==  altBase.lower():
			reads_bases_transformed += altBase.lower()
			nAltBases += 1
			i += 1
		# indels
		elif reads_bases[i] in ['+', '-', '*']:
			indel = ""
			if reads_bases[i] == '*':
				i += 1
				indel = '*'
			else:
				len_indel = int(re.search(r'\d+', reads_bases[i+1:i+3]).group())
				i += len_indel + len(str(len_indel)) + 1
				indel = reads_bases[i+2:i+2+len_indel].upper()
			if indel not in dIndels:
				dIndels[indel] = 1
			else:
				dIndels[indel] += 1
		# at begining
		elif reads_bases[i] == '^':
			i += 2
		# at ends, or gap
		elif reads_bases[i] in ['N', 'n', '$']:
			i += 1
		# other alleles
		else:
			if reads_bases[i] not in dMultiBases:
				dMultiBases[reads_bases[i]] = 1
			else:
				dMultiBases[reads_bases[i]] += 1
			reads_bases_transformed += reads_bases[i]
			i += 1
		nReadsBases += 1
	return reads_bases_transformed, nReadsBases, nRefBases, dMultiBases, dIndels
