import os
import sys
import argparse

import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
import rpy2.rlike.container as rlc

def making_plot(args):
	''' making Q-Q plot and Manhattan plot '''
	# install qqman package if not installed
	if not rpackages.isinstalled("qqman"):
		rutils = rpackages.importr('utils')
		rutils.chooseCRANmirror(ind=84)
		rutils.install_packages("qqman")

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
	snps_to_highlight = []
	with open(args.highlight_snps, 'r') as fHIGHLIGHT:
		for line in fHIGHLIGHT:
			snps_to_highlight.append(line.strip())

	# getting p-value cutoff given FDR level
	unadjust_p_cutoff = calculate_pval_cutoff([pvals[k] for k in sorted(pvals.iterkeys())], args.adj_cutoff)
	if args.pdf:
		out_qqplot = args.outp + ".fisher.qqplot.pdf"
		out_manhattan = args.outp + ".fisher.manhattan.pdf"
	elif args.png:
		out_qqplot = args.outp + ".fisher.qqplot.png"
		out_manhattan = args.outp + ".fisher.manhattan.png"
	grdevices = rpackages.importr('grDevices')
	# making Q-Q plot
	raw_pvals_vector = robjects.FloatVector([pvals[k] for k in sorted(pvals.iterkeys())])
	adjust_pvals_vector = robjects.FloatVector([adjust_pvals[k] for k in sorted(adjust_pvals.iterkeys())])
	make_qqplots(grdevices, raw_pvals_vector, out_qqplot)

	# maing Manhattan plot
	make_manhattan(grdevices, data, raw_pvals_vector, snps_to_highlight, out_manhattan)

def make_qqplots(grdevices, raw_pvals_vector, out_qqplot):
	''' making qqplot '''
	qqman = rpackages.importr('qqman')
	grdevices.pdf(out_qqplot)
#	robjects.r['par'](mfrow=robjects.IntVector([2,1]))
	qqman.qq(raw_pvals_vector, main="Q-Q plot for raw p-values using Fisher's Exact test")
	grdevices.dev_off()

def make_manhattan(grdevices, data, raw_pvals_vector, snps_to_highlight, out_manhattan):
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
	sig_snps = robjects.StrVector(snps_to_highlight)
	qqman = rpackages.importr('qqman')
	grdevices.pdf(out_manhattan)
#	robjects.r['par'](mfrow=robjects.IntVector([2,1]))
	qqman.manhattan(robjects.DataFrame(od_raw), highlight=sig_snps, col = color_vector, suggestiveline=3.716482, genomewideline=False, xlim=robjects.IntVector([20, 43]), ylim=robjects.IntVector([0,10]))
	grdevices.dev_off()


