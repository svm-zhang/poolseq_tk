import os
import sys
import collections

import sz_utils
from colortext import ColorText

def outVCFHeaders(samples, out):
	out.write("##fileformat=VCF4.1\n")
	out.write("##INFO=<ID=pval,Number=1,Type=Float,Description=\"raw p-value for each test\">\n")
	out.write("##INFO=<ID=corrPval,Number=1,Type=Float,Description=\"corrected p-value for each test\">\n")
	out.write("##INFO=<ID=ratio,Number=1,Type=Float,Description=\"odss ratio for each test\">\n")
	out.write("#CHROM\tPOS\tID\tREF\tALT\t"
			  "QUAL\tFILTER\tINFO\tFORMAT\t"
			  "%s\n" %(samples.replace(',', "\t")))

def getFilters(lfilters):
	dfilters = collections.defaultdict(tuple)
	for i in range(len(lfilters)):
		for j in range(len(lfilters[i].strip())):
			# for now only support two types of operators
			# I need a much clever way to parse expressions
			if lfilters[i][j] in ['>', '<']:
				operator = lfilters[i][j]
				literal = lfilters[i][:j].strip()
				value = float(lfilters[i][j+1:].strip())
				if literal not in dfilters:
					dfilters[literal] = (operator, value)
					break
	return dfilters

def getFst(ifst, dfst):
	with open(ifst, 'r') as fFST:
		for line in fFST:
			tmp_line = line.strip().split("\t")
			chr = tmp_line[0]
			pos = int(tmp_line[1])
			fst = tmp_line[2]
			if chr not in dfst:
				dfst[chr] = [(pos, fst)]
			else:
				dfst[chr].append((pos, fst))

	for k, v in dfst.items():
		v.sort()
	return dfst

def run_prepVCF(args):
	sz_utils.check_if_files_exist(args.infile)

	dfst = collections.defaultdict(list)
	if args.ifst:
		dfst = getFst(args.ifst, dfst)

	dfilters = getFilters(args.filters)
	outVCFHeaders(args.samples, args.out)
	with open(args.infile, "r") as fIN:
		for line in fIN:
			tmp_line = line.strip().split("\t")
			chr = "2L"
			pos = int(tmp_line[1])
			refBase = tmp_line[2]
			altBase = tmp_line[3]
			pval = float(tmp_line[-3])
			corrPval = float(tmp_line[-2])
			ratio = float(tmp_line[-1])
			fst = -1.0
			if chr in dfst:
				for j in range(len(dfst[chr])):
					if pos == dfst[chr][j][0]:
						fst = float(dfst[chr][j][1])
						dfst[chr].pop(j)
						break
			if "ratio" in dfilters:
				if ((dfilters["ratio"][0] == '<' and ratio >= dfilters["ratio"][1]) or
					 dfilters["ratio"][0] == '>' and ratio <= dfilters["ratio"][1]):
					continue
			if "pval" in dfilters:		# fix later
				pass
			if "corrPval" in dfilters:		# fix later
				pass
			args.out.write("2L\t%s\t.\t%s\t%s\t.\t.\t"
							%(pos, refBase, altBase))
			if fst == -1.0:
				args.out.write("pval=%.5g;corrPval=%.5f;ratio=%.5f;fst=N/A\t"
								%(pval, corrPval, ratio))
			else:
				args.out.write("pval=%.5g;corrPval=%.5f;ratio=%.5f;fst=%.5f\t"
								%(pval, corrPval, ratio, fst))
			args.out.write("GT:Table")
			poolIndex = 1
			for i in range(len(tmp_line[4:-3])):
				table = tmp_line[4:-3][i].replace(':', '-')
				args.out.write("\t./.:%s" %(table))
			args.out.write("\n")
	args.out.close()
