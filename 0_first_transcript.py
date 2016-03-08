#!/usr/bin/python3 

import kang,sys

file_gff 	= sys.argv[1] #'Pepper.v.1.5.total.gff3'
file_fa		= sys.argv[2] #'Pepper.v.1.5.scaffold.fa'

dicHD2seq = kang.Fasta2dic(file_fa)

class transcript():
	def _init_(self):
		self.strChr	= ''
		self.strID	= ''	
		self.strTN	= '' # transcript name
		self.strCDS_list	= [] # coding sequence location [[l1,l2],[l1,l2],..]
		self.strSTRD = ''

dicT2cls = {} # trancript to class	
for line in open(file_gff):
	if line.strip() == '':
		continue
	if line[0] == '#':
		continue
	cell = line.strip().split('\t')
	strChr	= cell[0]
	strT	= cell[2]
	strL1	= cell[3]
	strL2	= cell[4]
	strST	= cell[6] # strand
	strINFO	= cell[8]
	dicINFO = {}
	for each in strINFO.split(';'):
		try:
			key 	= each.split('=')[0]
			val	= each.split('=')[1]
			dicINFO[key] = val
		except IndexError:
			pass
	if strT == 'mRNA':
		if dicINFO['longest'] == '1':
			pass
		else: continue
		strID = dicINFO['ID']
		strTN = dicINFO['Name']
		try:
			cT = dicT2cls[strID]
			print('!!!!')
			exit()
		except KeyError:
			cT 		= transcript()
			cT.strTN 	= strTN
			cT.strSTRD 	= strST
			cT.strChr	= strChr
			cT.strCDS_list = []
			dicT2cls[strID] = cT
	elif strT == 'CDS':
		strID = dicINFO['Parent']
		try:
			cT = dicT2cls[strID]
			cT.strCDS_list.append([int(strL1),int(strL2)])
		except KeyError:
			#print('!!')
			#exit()
			continue
	

Outfile = open(file_fa+'.cds.fa','w')
for ID in dicT2cls:
	cT 		= dicT2cls[ID]
	cds_locs 	= cT.strCDS_list
	cds_locs.sort(key = lambda x : x[0])
	cdsseq_pre = ''

	start_loc 	= cds_locs[0]
	intL1, intL2 	= start_loc
	try:
		if dicHD2seq[cT.strChr]:
			pass
	except KeyError:
		continue
	#cdsseq_pre 	= dicHD2seq[cT.strChr][intL1-1-1000:intL1-1] 
	for cds_loc in cds_locs:
		intL1, intL2 = cds_loc
		cdsseq_pre += dicHD2seq[cT.strChr][intL1-1:intL2]
	end_loc		= cds_locs[-1]
	intL1, intL2	= end_loc
	#cdsseq_pre += dicHD2seq[cT.strChr][intL2:intL2+1000]

	if cT.strSTRD == '+':
		print('>'+cT.strTN,file=Outfile)
		print(cdsseq_pre,file=Outfile)
	elif cT.strSTRD == '-':
		print('>'+cT.strTN,file=Outfile)
		print(kang.rev_comp(cdsseq_pre),file=Outfile)

