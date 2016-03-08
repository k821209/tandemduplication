#!/usr/bin/python3 
import sys
file_interpro = sys.argv[1]
dicGN2GO = {}
Outfile = open(file_interpro.split('/')[-1]+'.bingo','w')
for line in open(file_interpro):
	cell = line.strip().split('\t')
	strGN = cell[1]
	try:
		strGO = cell[9].split(',')
		if strGO == []:
			continue
		try:
			dicGN2GO[strGN] += strGO
		except KeyError:
			dicGN2GO[strGN] = strGO
	except IndexError:
		pass
for strGN in dicGN2GO:
	for eachGO in set(dicGN2GO[strGN]):
		if eachGO == '':
			continue
		print(strGN,'=',eachGO.replace('GO:',''),file=Outfile)
	
