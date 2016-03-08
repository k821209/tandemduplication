#!/usr/bin/python3 

import kang,sys
#uniprot_sprot.fasta can be downloaded from uniprot webpage
file_fa = 'uniprot_sprot.fasta'
file_bo	= sys.argv[1]#'uniprot_sprot.fasta.bx.na1.out7'
Outfile	= open(file_bo+'.annot','w')

dicHD2seq 	= kang.Fasta2dic_all(file_fa)
dicID2des	= {}
for strHD in dicHD2seq:
	strID 	= strHD.split()[0]
	strDES	= ','.join(strHD.split()[1:])
	dicID2des[strID] = strDES

for line in open(file_bo):
	if line[0] == '#':
		continue
	cell = line.strip().split('\t')
	strQ = cell[0]	
	strT = cell[1]
	print(strQ,strT,dicID2des[strT],','.join(cell[2:]),sep='\t',file=Outfile)
Outfile.close()
