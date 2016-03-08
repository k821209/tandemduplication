#!/usr/bin/python3

import pickle,sys
try:
	dicSig2annot = pickle.load( open( "../Sig2Annot.p", "rb" ) )
except IOError:
	dicSig2annot = {}
Outfile = open(sys.argv[1]+'.parsed.csv','w')
for line in open(sys.argv[1]): # Gmax_275_v2.0.fa.ver2_tabout.csv
	if line[0:7] == 'element':
		continue
	cell 	= line.strip().split('\t')
	strSL 	= cell[0] 	# start location
	strEL	= cell[1] 	# end location
	strChr	= cell[3]
	sig_list = cell[-1].split('/')
	if len(sig_list) > 1:
		pass
	else:
		print(strChr,strSL,strEL,'no_sig','no_annot',sep='\t',file=Outfile) 
		continue
	sig_list = [x.split('_')[0] for x in sig_list]
	sig = ';'.join(sig_list)
	try:
		print(strChr,strSL,strEL,sig,dicSig2annot[sig],sep='\t',file=Outfile)
	except KeyError:
		if False in [(x in sig_list) for x in 'AP;INT;RT;RNaseH;GAG'.split(';')]:
			value = 'Other'
		else:
			print(sig)
			print(set(list(dicSig2annot.values())))
			value = input("suggestion?")
		dicSig2annot[sig] = value
		print(strChr,strSL,strEL,sig,value,sep='\t',file=Outfile)
pickle.dump(dicSig2annot, open( "../Sig2Annot.p", "wb" ) )	
	
