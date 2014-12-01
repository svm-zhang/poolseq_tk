import os
import sys
import argparse

import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
import rpy2.rlike.container as rlc

def getopts():
	usage = "Making Manhattan plot and QQ plot"
	plot_parser = argparse.ArgumentParser(description=usage)
	plot_parser.add_argument("-i", metavar="FILE", dest="input", required=True, help="input file of test results with all SNPs, e.g. *.fisher.all, *.cmh.all")
	plot_parser.add_argument("-interests", metavar="FILE", dest="interests_snps", help="file of a list of SNPs of interests. These SNPs will be highlighted in the Manhattan plot")
	plot_parser.add_argument("-outp", metavar="PREFIX", dest="outp", help="prefix of output file")
	plot_parser.add_argument("-adj_cutoff", metavar="FLOAT", dest="adj_cutoff", type=float, default=0.05, help="specify the cutoff below which adjusted p-values will be considered as significant")
	plot_mutual_group = plot_parser.add_mutually_exclusive_group(required=True)
	plot_mutual_group.add_argument("-pdf", dest="pdf", action="store_true", help="output qqplot in pdf format")
	plot_mutual_group.add_argument("-png", dest="png", action="store_true", help="output qqplot in pdf format")
	return plot_parser.parse_args()

def making_plot():
	''' making Q-Q plot and Manhattan plot '''

	args = getopts()

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


