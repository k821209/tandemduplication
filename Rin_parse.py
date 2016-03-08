#!/usr/bin/python3 

import sys

Outfile = open(sys.argv[1]+'.Rin','w')
for line in open(sys.argv[1]): # Alyrata_107_cds_primaryTranscriptOnly.fa.selfmegablast.eval1e100.out7.tandemNnontandem.Ksadd
	cell = line.strip().split('\t')
	try:
		for each in cell[7].split(','):
			print(sys.argv[2], each,sep='\t',file=Outfile) # Species name
	except IndexError:
		pass
	
	
