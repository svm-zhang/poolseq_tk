import os
import sys
import argparse
import collections

import sz_utils
from colortext import ColorText

#def run_collapse(args):
def run_collapse(isnp, m1, m2, out):
	'''
		Given two pileup files of the same region, like 2l+ and 2la,
		collapse the pileups at each corresponding SNP
		Some SNPs are not reported in one or the other pileup file.
		A full list of SNP positions are required
	'''
	m1_base = os.path.basename(m1)
	m2_base = os.path.basename(m2)

	# first, getting the full list of SNPs
	dSNPs = sz_utils.getSNPs(isnp)

	offset1 = 0
	offset2 = 20524057

	# second, reading each of the pileup files
	chr1, m1_info = read_mpileup(m1, offset1)
	chr2, m2_info = read_mpileup(m2, offset2)

	ColorText().info("[poolseq_tk] %s: %d SNPs parsed\n" %(m1_base, len(m1_info)), "stderr")
	ColorText().info("[poolseq_tk] %s: %d SNPs parsed\n" %(m2_base, len(m2_info)), "stderr")

#	fOUT = None
#	if args.out != sys.stdout:
#		outdir = os.path.dirname(os.path.realpath(args.out))
#		sz_utils.make_dirs_if_necessary(outdir)
#		fOUT = open(args.out, 'w')
#	else:
#		fOUT = args.out
	fOUT = open(out, 'w')
	ColorText().info("[poolseq_tk]: collapsing mpileups %s and %s ..."
					 %(m1_base, m2_base), "stderr")
	for (chr,pos) in sorted(dSNPs.iterkeys()):
		k = (chr, pos)
		reads_bases_collapsed = ""
		refBase1, refBase2 = "", ""
		if pos in m1_info:
			if m1_info[pos][0] == dSNPs[k][0]:
				refBase1 = m1_info[pos][0]
			else:
				ColorText().error("SNP position: %d %s\t\tMpileup position: %d %s\n"
								  %(pos, dSNPs[k][0], pos, m1_info[pos][0]))
				sys.exit()
		else:
			refBase1 = ""
		if pos in m2_info:
			if m2_info[pos][0] == dSNPs[k][1]:
				refBase2 = m2_info[pos][0]
			else:
				ColorText().error("SNP position: %d %s\t\tMpileup position: %d %s\n"
								  %(pos, dSNPs[k][1], pos, m2_info[pos][0]))
				sys.exit()
		else:
			refBase2 = ""

		if refBase1 != "" and refBase2 != "":
			reads_bases_collapsed1, nReadsBases1, nRefBases1, dMultiBases1, dIndels1 = sz_utils.parseReadsBases(m1_info[pos][1],
																												refBase1,
																												refBase2)
			reads_bases_collapsed2, nReadsBases2, nRefBases2, dMultiBases2, dIndels2 = sz_utils.parseReadsBases(m2_info[pos][1],
																												refBase2,
																												refBase1)
			reads_bases_collapsed = reads_bases_collapsed1 + reads_bases_collapsed2
			nReadsBases = nReadsBases1 + nReadsBases2
			nRefBases = nRefBases1 + (nReadsBases2-nRefBases2)
			dMultiBases = dict(dMultiBases1.items() + dMultiBases2.items())
			dIndels = dict(dIndels1.items() + dIndels2.items())
			nMultiBases = sum(dMultiBases.values()) + sum(dIndels.values())
			if (nReadsBases == nRefBases or
				nMultiBases <= 1):
				fOUT.write("%s\t%d\t%s\t%s\t%s\n" %(chr, pos, refBase1, refBase2,
													reads_bases_collapsed))
		elif refBase1 == "" and refBase2 != "":
			reads_bases_collapsed2, nReadsBases2, nRefBases2, dMultiBases2, dIndels2 = sz_utils.parseReadsBases(m2_info[pos][1],
																												refBase2,
																												refBase1)
			nMultiBases = sum(dMultiBases2.values()) + sum(dIndels2.values())
			if (nReadsBases2 == nRefBases2 or
				nMultiBases <= 1):
				fOUT.write("%s\t%d\t%s\t%s\t%s\n" %(chr, pos, dSNPs[k][0], refBase2,
													reads_bases_collapsed2))
		elif refBase1 != "" and refBase2 == "":
			reads_bases_collapsed1, nReadsBases1, nRefBases1, dMultiBases1, dIndels1 = sz_utils.parseReadsBases(m1_info[pos][1],
																												refBase1,
																												refBase2)
			nMultiBases = sum(dMultiBases1.values()) + sum(dIndels1.values())
			if (nReadsBases1 == nRefBases1 or
				nMultiBases <= 1):
				fOUT.write("%s\t%d\t%s\t%s\t%s\n" %(chr, pos, refBase1, dSNPs[k][1],
													reads_bases_collapsed1))
	ColorText().info(" [done]\n", "stderr")
	fOUT.close()

def read_mpileup(mpileup_file, offset):
	''' read certain columns in a pileup file into a dictionary of tuple '''
	ColorText().info("[poolseq_tk]: reading %s ..." %(mpileup_file), "stderr")
	mpileup_info = collections.defaultdict(tuple)
	chr = ""
	sz_utils.check_if_files_exist(mpileup_file)
	with open(mpileup_file, 'r') as fMPILEUP:
		for line in fMPILEUP:
			tmp_line = line.strip().split("\t")
			chr = tmp_line[0]
			'''
				key: SNP position in integer
				value: a tuple with two elements
					   1) ref base at the position
					   2) reads bases covering that position
						  if the coverage at that position > 0
						  else N/A
			'''
#			if int(tmp_line[3]) == 0:		# the forth column: coverage at a position
#				mpileup_info[int(tmp_line[1])+offset] = (tmp_line[2].upper(), "N/A")
			if int(tmp_line[3]) > 0:
				mpileup_info[int(tmp_line[1])+offset] = (tmp_line[2].upper(), tmp_line[4])
	ColorText().info(" [done]\n", "stderr")
	return chr, mpileup_info

def main():
	isnp = sys.argv[1]
	m1 = sys.argv[2]
	m2 = sys.argv[3]
	out = sys.argv[4]

	sz_utils.make_dirs_if_necessary(out)
	run_collapse(isnp, m1, m2, out)

if __name__ == "__main__":
	main()
