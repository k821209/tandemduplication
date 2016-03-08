#!/usr/bin/python3 

import glob,subprocess, sys, kang, numpy, os

#Cluster        Rep. genename   Annotation      Number of homologs      List of homologs
#Homology cluster 0      Tandem_cluster 0        AT4G15330       cytochrome P450, putative, expressed    4       AT4G15330,AT4G15360,AT4G15380,AT4G15350
#Homology cluster 0      Tandem_cluster 1        AT3G20080       cytochrome P450, putative, expressed    6       AT3G20080,AT3G20120,AT3G20090,AT3G20110,AT3G20130,AT3G20100
#Homology cluster 0      Tandem_cluster 2        AT3G20950       cytochrome P450, putative, expressed    3       AT3G20950,AT3G20940,AT3G20960
#Homology cluster 0      Tandem_cluster 3        AT2G27000       cytochrome P450, putative, expressed    2       AT2G27000,AT2G27010

outdir = './output'
os.system('mkdir ./output')
def fasta2dic(file_cds):
	bulk = open(file_cds).read()
	bulk_list	= bulk.split('>')
	dic	= {}
	for each in bulk_list:
		strHD	= each.split('\n')[0].split('|')[0]
		strSeq	= ''.join(each.split('\n')[1:])
		dic[strHD] = strSeq	
	return(dic)


def convert_axt(filename):
	file_list = glob.glob('%s/%s'%(outdir,filename)) # plank result in fas, second result
	Outfile = open('%s/%s.result_of_2.axt'%(outdir,filename),'w')
	for each_file in file_list:
		result = each_file + '\n'
		bulk = open(each_file).read()
		bulk_list = bulk.split('>')
		dicH2Seq = {} # header2seq dictionary
		for bulk_each in bulk_list:
			if bulk_each.strip() == '':
				continue
			header = bulk_each.split('\n')[0]
			seq = ''.join(bulk_each.split('\n')[1:])
			dicH2Seq[header] = seq
		for dicH2Seq_key in dicH2Seq:
			dicH2Seq_value = dicH2Seq[dicH2Seq_key]
			result += dicH2Seq_value + '\n'
		print(result.strip()+'\n',file=Outfile)
	Outfile.close()

	
def get_Ks_list(seq1name,seq2name,seq1,seq2):
	file_temp = open('%s/temp.%s.%s.fa'%(outdir,seq1name,seq2name),'w')
		
	print('>seq1',file=file_temp)
	print(seq1 ,file=file_temp)
	print('>seq2',file=file_temp)	
	print(seq2 ,file=file_temp)
	file_temp.close()
		

	subprocess.call('/data/program/prank/bin/prank -quiet -DNA -translate -d=%s/%s -o=%s/%s.plankout '\
			%(outdir,   	'temp.%s.%s.fa'%(seq1name,seq2name),   	outdir,  	'temp.fa'),shell=True)		

	convert_axt('temp.fa.plankout.best.nuc.fas')

	subprocess.call('/data/program/KaKs_Calculator1.2/src/KaKs_Calculator -i %s/temp.fa.plankout.best.nuc.fas.result_of_2.axt -o %s/%s.%s.result_of_2.axt.LWL -m LWL '\
			%(outdir,	outdir,seq1name,seq2name),	shell=True)
	
	result = open('%s/%s.%s.result_of_2.axt.LWL'%(outdir,seq1name,seq2name)).read()
	strKs  = result.split('\n')[1].split('\t')[3]
	return(strKs)

#Cluster        Rep. genename   Annotation      Number of homologs      List of homologs
#Homology cluster 0      Tandem_cluster 0        AT4G15330       cytochrome P450, putative, expressed    4       AT4G15330,AT4G15360,AT4G15380,AT4G15350
#Homology cluster 0      Tandem_cluster 1        AT3G20080       cytochrome P450, putative, expressed    6       AT3G20080,AT3G20120,AT3G20090,AT3G20110,AT3G20130,AT3G20100
#Homology cluster 0      Tandem_cluster 2        AT3G20950       cytochrome P450, putative, expressed    3       AT3G20950,AT3G20940,AT3G20960
#Homology cluster 0      Tandem_cluster 3        AT2G27000       cytochrome P450, putative, expressed    2       AT2G27000,AT2G27010

file_cds	= sys.argv[1] 	# cds file
dicGN2seq 	= fasta2dic(file_cds)

def main():
	file_in	= sys.argv[2]	# example table  
	Outfile	= open(file_in+'.Ksadd','w')
	for line in open(file_in):
		if line[0] == '#':
			continue	
		cell = line.strip().split('\t')
		strHC	= cell[0] # Homology cluster
		strSC	= cell[1] # Sub cluster
		if 'Tandem' in strSC:
			pass
		else: continue
		strGN	= cell[2]
		strDes	= cell[3] # Description
		strNH	= cell[4] # Number of homologs
		strGN_list	= cell[5].split(',')
		fKs_list	= []
		i = 1
		for gene1 in strGN_list:
			for gene2 in strGN_list[i:]:
				j = 0 
				while 1:
					try:
						seq1	= dicGN2seq[gene1] # modify accoridng to the genename format
						break
					except KeyError:
						j += 1
				j = 0 
				while 1:
					try:
						seq2	= dicGN2seq[gene2]
						break
					except KeyError:
							j += 1
				strKs = get_Ks_list(gene1,gene2,seq1,seq2)
				if strKs != 'NA':
					fKs_list.append(float(strKs))
			i += 1
		fKs_median	= numpy.median(fKs_list)
		print(line.strip(),fKs_median,','.join(map(str,fKs_list)),sep='\t',file=Outfile)
main()
