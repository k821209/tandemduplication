#!/usr/bin/python3

import sys

Infile       = 'Athaliana_167.fa.pep.fa.bp.ev1e5.out7'
Outfile      = open(Infile+'.tandemNnontandem.ver2.out','w')
Outfile_test = open('test','w')
Infile_annot = '/data/ref/Athaliana_167_annotation_info.txt'
Infile_gff   = '/data/ref/Athaliana_167_gene.gff3'
intSimCut    = 80
intGlCut     = 200
file_fa      = 'Athaliana_167.fa.pep.fa'
########################################################################
def Fasta2dic(file_fasta):
	dic     = {}
	ofile	= open(file_fasta)
	bulk    = ofile.read()
	ofile.close()
	bulk_list = bulk.split('>')
	for each in bulk_list:
		if each == '':
			continue
		strHD  = each.split('\n')[0].split()[0].split('|')[0]
		strSeq = ''.join(each.split('\n')[1:])
		dic[strHD] = strSeq
	return(dic)
#########################################################################
dicHD2seq = Fasta2dic(file_fa)
#########################################################################
dicGN2LEN   = {}
for line in open(Infile_gff):
	if line[0] == '#':
		continue
	cell               = line.strip().split('\t')
	strChr             = cell[0]
	strTP              = cell[2]
	intLOC1            = int(cell[3])
	intLOC2            = int(cell[4])
	if strTP == 'gene':
		strGN              = cell[-1].split(';')[1].replace('Name=','') 
		intLN              = intLOC2 - intLOC1
		dicGN2LEN[strGN]   = [strChr,intLOC1,intLN]
dicGN2LEN_list = list(dicGN2LEN.keys())
dicGN2LEN_list.sort(key=lambda x:dicGN2LEN[x][1])
dicGN2LEN_list.sort(key=lambda x:dicGN2LEN[x][0])
i = 0
for dicGN2LEN_key in dicGN2LEN_list:
	dicGN2LEN[dicGN2LEN_key].append(i)
	i += 1
#########################################################################
dicGN2ANNOT = {}
for line in open(Infile_annot):
	if line[0] == '#':
		continue
	cell               = line.strip().split('\t')
	strGN              = cell[0]
	print(strGN)
	#strNT              = cell[2].split('.')[1]
	strAnnot           = cell[-1]
	try:
		if dicGN2ANNOT[strGN]:
			pass
	except KeyError:
		dicGN2ANNOT[strGN] = strAnnot
#########################################################################
pre_HMs_list = []
HMs_list     = [] # Genename 2 homologs
i = 0
for line in open(Infile):
	if i % 500 == 0 :
		print(i)
	i += 1
	if line[0] == '#':
		continue
	cell  = line.strip().split('\t')
	strGN = cell[0].split('|')[0]
	strHM = cell[1].split('|')[0]
	fSim	= float(cell[2])
	intAL = int(cell[3])
	intGL = len(dicHD2seq[strGN])
	#########################################################################
	#########################################################################
	if strGN == strHM:
		continue
	#if intAL < intAL_limit:
	#	continue
	if intGL < intGlCut:
		continue
	if fSim < intSimCut:
		continue
	if intAL/intGL < 0.8: ## intGL : Gene length
		continue
	pre_HMs_list.append([strGN,strHM,0])
i = 0
for HM in pre_HMs_list:
	if pre_HMs_list[i][2] == 1:
		i += 1
		continue
	pre_HMs_list[i][2] = 1
	i += 1
	preexist = 0
	strGN      = HM[0]
	strHM      = HM[1]
	pre_len = 0
	HMs_pooling = [strGN,strHM]
	while pre_len != len(HMs_pooling):
		j = 0
		print(pre_len,len(HMs_pooling))
		pre_len = len(HMs_pooling)
		for each in pre_HMs_list:
			if each[0] in HMs_pooling and each[1] not in HMs_pooling:
				HMs_pooling.append(each[1])
				if pre_HMs_list[j][2] == 1:
					print('error')
					exit()
				pre_HMs_list[j][2] = 1
			elif each[0] not in HMs_pooling and each[1] in HMs_pooling:
				HMs_pooling.append(each[0])
				if pre_HMs_list[j][2] == 1:
					print('error')
					exit()
				pre_HMs_list[j][2] = 1
			elif each[0] in HMs_pooling and each[1] in HMs_pooling:
				pre_HMs_list[j][2] = 1
				pass
			j += 1
	HMs_list.append(HMs_pooling)		
			
			
			
#########################################################################
	
def homolog_parse(HM_list):
	tandem_list = []
	nontandem_list = []
	dicCHR2GNs = {}
	for strTGN in HM_list:	
		strGN = strTGN.split('.')[0]
		strChr = dicGN2LEN[strGN][0]
		intORD = dicGN2LEN[strGN][3]
		try:
			dicCHR2GNs[strChr].append([strTGN,intORD,0])
		except KeyError:
			dicCHR2GNs[strChr] = [[strTGN,intORD,0]]
	for dicCHR2GNs_key in dicCHR2GNs:
		dicCHR2GNs_value = dicCHR2GNs[dicCHR2GNs_key]
		i = 0
		for GN_infos in dicCHR2GNs_value:
			print(dicCHR2GNs_key,dicCHR2GNs_value,file=Outfile_test)
			strGN  = GN_infos[0]
			intORD = GN_infos[1]
			bCheck = GN_infos[2]
			if bCheck == 1:
				continue
			temp_list = [strGN]
			j = 0 + i
			for target_GN_infos in dicCHR2GNs_value[i:]:
				strTGN  = target_GN_infos[0]
				intTORD = target_GN_infos[1]
				bTCheck = target_GN_infos[2]
				if bTCheck == 1:
					j += 1
					continue
				for temp_each in temp_list:
					strGN = temp_each.split('.')[0]
					temp_each_ord = dicGN2LEN[strGN][3]
					if 0 < abs(temp_each_ord-intTORD) <= 10:
						temp_list.append(strTGN)
						dicCHR2GNs_value[j][2] = 1 # Add bCheck
						break
				j += 1
			if len(temp_list) == 1:
				pass
			else : 
				dicCHR2GNs_value[i][2] = 1 # Add bCheck
				tandem_list.append(temp_list)
			i += 1
		dicCHR2GNs[dicCHR2GNs_key] = dicCHR2GNs_value
	for dicCHR2GNs_key in dicCHR2GNs:
		dicCHR2GNs_value = dicCHR2GNs[dicCHR2GNs_key]
		for each in dicCHR2GNs_value:
			if each[2] == 0:
				nontandem_list.append(each[0])
	return(tandem_list,nontandem_list)
		


def main():
	i = 0
	print('#Cluster','Rep. genename','Annotation','Number of homologs','List of homologs',sep='\t',file=Outfile)
	HMs_list_parsed = []
	HMs_list.sort(key=lambda x:len(x),reverse=True)
	t = 0
	nt = 0
	for each in HMs_list:
		tandem_list, nontandem_list = homolog_parse(each)
		print(tandem_list,file=Outfile_test)
		print('Homology cluster %d'%i,'whole homologs',each[0],dicGN2ANNOT[each[0]],len(each),','.join(each),file=Outfile,sep='\t')

		for tandem in tandem_list:
			tandem.sort(key=lambda x:len(dicHD2seq[x]),reverse=True)
			nontandem_list.append(tandem[0])
			if len(tandem) <2:
				continue
			print('Homology cluster %d'%i,'Tandem_cluster %d'%t,tandem[0],dicGN2ANNOT[tandem[0]],len(tandem),','.join(tandem),file=Outfile,sep='\t')
			t += 1
		nontandem_list.sort(key=lambda x:len(dicHD2seq[x]),reverse=True)
		try:
			if len(nontandem_list) < 2:
				i += 1
				continue
			print('Homology cluster %d'%i,'Nontandem_cluster %d'%nt,nontandem_list[0],dicGN2ANNOT[nontandem_list[0]],len(nontandem_list),','.join(nontandem_list),file=Outfile,sep='\t')
			nt += 1
			
		except:
			pass
		i += 1 
main()

	
