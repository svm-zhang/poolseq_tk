'''
	poolseq_tk.py plot
	Description: Taking the test resutls and make Q-Q and Manhattan plots
	Author: Simo V. Zhang

	Input: tests results obtained from running:
		   python poolseq_tk.py fisher or
		   python poolseq_tk.py cmh

	Output: Q-Q and Manhattan plots in either PDF or PNG
'''

import os
import sys
import argparse
import collections
import math

from colortext import ColorText
import sz_utils

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
	ColorText().info("[poolseq_tk]: Extracting P-Values ... ", "stderr")
	data = collections.defaultdict(tuple)
	pvals, adjust_pvals = {}, {}
	with open(args.input, 'r') as fIN:
		for line in fIN:
			tmp_line = line.strip().split("\t")
			data[int(tmp_line[1])] = ("_".join(tmp_line[0:2]), int(tmp_line[0][0]), int(tmp_line[1]))
			pvals[int(tmp_line[1])] = float(tmp_line[-3])
	ColorText().info(" [done]\n", "stderr")

	# get FDR cutoff using BH if not provided through command line
	pcutoff = 0.0
	if not args.pcutoff:
		ColorText().info("[poolseq_tk]: Getting p-value cutoff at FDR %d%%: " %(args.fdrlevel*100), "stderr")
		pcutoff = sz_utils.getFDR_BH(pvals, args.fdrlevel)
		ColorText().info("%.5e\n" %(pcutoff), "stderr")
	else:
		pcutoff = args.pcutoff
		ColorText().info("[poolseq_tk]: p-value cutoff provided: %.5e\n" %(pcutoff), "stderr")

	# get SNPs to highlight
	snps_to_highlight = []
	if args.highlight_snps:
		ColorText().info("[poolseq_tk]: Getting SNPs to be highlighed in Manhattan plot ... ", "stderr")
		with open(args.highlight_snps, 'r') as fHIGHLIGHT:
			for line in fHIGHLIGHT:
				tmp_line = line.strip().split("\t")
				snps_to_highlight.append('_'.join(tmp_line[:2]))
	ColorText().info(" [done]\n", "stderr")

	if args.pdf:
		out_qqplot = args.outp + ".qqplot.pdf"
		out_manhattan = args.outp + ".manhattan.pdf"
	elif args.png:			# save to PNG probably wont work
		out_qqplot = args.outp + ".qqplot.png"
		out_manhattan = args.outp + ".manhattan.png"
	sz_utils.make_dirs_if_necessary(out_qqplot, out_manhattan)
	grdevices = rpackages.importr('grDevices')
	raw_pvals_vector = robjects.FloatVector([pvals[k] for k in sorted(pvals.iterkeys())])

	ColorText().info("[poolseq_tk]: Making Q-Q plot ...", "stderr")
	make_qqplots(grdevices, raw_pvals_vector, out_qqplot, args.qqtitle)
	ColorText().info(" [done]\n", "stderr")

	ColorText().info("[poolseq_tk]: Making Manhattan plot ...", "stderr")
	make_manhattan(grdevices, data, raw_pvals_vector,
				   snps_to_highlight, pcutoff,
				   out_manhattan, args.mantitle, args.manx, args.manxlim)
	ColorText().info(" [done]\n", "stderr")

def make_qqplots(grdevices, raw_pvals_vector, out_qqplot, title=""):
	''' making qqplot '''
	qqman = rpackages.importr('qqman')
	grdevices.pdf(out_qqplot)
	qqman.qq(raw_pvals_vector, main=title)
	grdevices.dev_off()

def make_manhattan(grdevices, data, raw_pvals_vector,
				   snps_to_highlight, padj_cutoff,
				   out_manhattan, title="", xlable="", xlim="-"):
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
	if xlim != "-":
		xmin = int(xlim.split(",")[0])
		xmax = int(xlim.split(",")[1])
		qqman.manhattan(robjects.DataFrame(od_raw), highlight=sig_snps,
						col = color_vector, suggestiveline=False,
						genomewideline=-1*math.log10(padj_cutoff),
						xlim=robjects.IntVector([xmin, xmax]),
						xlab=xlable, main=title)
	else:
		qqman.manhattan(robjects.DataFrame(od_raw), highlight=sig_snps,
						col = color_vector, suggestiveline=False,
						genomewideline=-1*math.log10(padj_cutoff),
						xlab=xlable, main=title)
	grdevices.dev_off()
