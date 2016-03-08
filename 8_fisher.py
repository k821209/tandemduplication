#!/usr/bin/python3
import sys
import numpy as np
from scipy import stats
np.seterr(divide='ignore', invalid='ignore')
file_gff 		= sys.argv[1]#'Vvinifera_145_gene.gff3'
file_csv		= sys.argv[2]#'Vvinifera_145.fa.ver2_tabout.csv.parsed.csv'
file_tandem		= sys.argv[3]#'Vvinifera_145.fa.pep.fa.bp.ev1e5.out7.tandemNnontandem.ver2.out.Ksadd'
spcs			= sys.argv[4]# 'Vvi'

class fisherinput:
	def __init__(self,strANN):
		self.repeatname = strANN
		self.tandem_w_repeat 	= 0
		self.tandem_wo_repeat 	= 0
		self.nontandem_w_repeat 	= 0
		self.nontandem_wo_repeat 	= 0
	
dicTandemlist = {}

for line in open(file_tandem):
	cell = line.strip().split('\t')
	gene_list = cell[5].split(',')
	for gene in gene_list:
		dicTandemlist[gene] = 1

dicGN2loc = {}
for line in open(file_gff):
	if line[0] == '#':
		continue
	cell = line.strip().split('\t')
	strChr = cell[0]
	intL1 = int(cell[3])
	intL2 = int(cell[4])
	strT = cell[2] 
	strInfo_list = cell[8].split(';')
	dicInfo = dict(zip([x.split('=')[0] for x in strInfo_list],[x.split('=')[1] for x in strInfo_list]))
	if strT == 'mRNA':
		dicGN2loc[dicInfo['Name']] = (strChr,intL1,intL2)
	
span   = 2000
window = 100000
dicChr2LTR = {}
for line in open(file_csv):
	cell = line.strip().split('\t')
	strChr = cell[0].split('_')[0]
	intL1 = int(cell[1])
	intL2 = int(cell[2])
	try:
		dicChr2LTR[(strChr,int(intL1/window))].append(cell)
	except KeyError:
		dicChr2LTR[(strChr,int(intL1/window))] = [cell]

def dicif(dic,key):
	try:
		if dic[key]:
			return(1)
	except KeyError:
		return(0)
dicRP2Class = {}

for strGN in dicGN2loc:
	strChr,intL1,intL2 = dicGN2loc[strGN]
	intL1 = intL1 - span 
	intL2 = intL2 + span
	bLTR = 0
	try: 
		if dicChr2LTR[(strChr,int(intL1/window))]:
			for ecell in dicChr2LTR[(strChr,int(intL1/window))]:
				eL1 = int(ecell[1])
				eL2 = int(ecell[2])
				eSig = ecell[3]
				eAnn = ecell[4]
				if eL1 < intL1 < eL2 or eL1 < intL2 < eL2:
					bLTR = 1
					if dicif(dicTandemlist,strGN):
						try:
							dicRP2Class[eAnn].tandem_w_repeat += 1
						except KeyError:
							cRP = fisherinput(eAnn)
							dicRP2Class[eAnn] = cRP
							dicRP2Class[eAnn].tandem_w_repeat += 1
					else:
						try:
							dicRP2Class[eAnn].nontandem_w_repeat += 1
						except KeyError:
							cRP = fisherinput(eAnn)
							dicRP2Class[eAnn] = cRP
							dicRP2Class[eAnn].nontandem_w_repeat += 1
	except KeyError:
		pass

total_tandem = len(list(dicTandemlist.keys()))
total_genes  = len(list(dicGN2loc.keys()))
total_nontandem = total_genes - total_tandem
#print("repeatname","tandem_w_repeat","nontandem_w_repeat","tandem_wo_repeat","cRP.nontandem_wo_repeat",spcs,sep='\t')
for strRP in dicRP2Class:
	cRP = dicRP2Class[strRP]
	tandem_wo_repeat = total_tandem - cRP.tandem_w_repeat
	nontandem_wo_repeat = total_nontandem - cRP.nontandem_w_repeat
	oddsratio, pvalue = stats.fisher_exact([[cRP.tandem_w_repeat,cRP.nontandem_w_repeat], [tandem_wo_repeat, nontandem_wo_repeat]])
	print(cRP.repeatname,cRP.tandem_w_repeat,cRP.nontandem_w_repeat,tandem_wo_repeat,nontandem_wo_repeat,spcs,oddsratio, pvalue,sep='\t')
					
			
