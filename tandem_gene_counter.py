#!/usr/bin/python3
#Homology cluster 0      whole homologs  AT4G20540.1     TKL_IRAK_DUF26-lc.12 - DUF26 kinases have homology to DUF26 containing loci, expressed  15      AT4G20540.1,AT4G20630.1,AT4G20530.1,AT4G20650.1,AT4G20640.1,AT4G20620.1,AT4G20600.1,AT4G20580.1,AT4G20610.1,AT4G20590.1,AT4G20570.1,AT4G20550.1,AT4G20560.1,AT4G20670.1,AT3G21990.1

#Homology cluster 0      Tandem_cluster 0        AT4G20650.1     TKL_IRAK_DUF26-lc.12 - DUF26 kinases have homology to DUF26 containing loci, expressed  14      AT4G20650.1,AT4G20640.1,AT4G20620.1,AT4G20600.1,AT4G20580.1,AT4G20610.1,AT4G20590.1,AT4G20570.1,AT4G20550.1,AT4G20560.1,AT4G20540.1,AT4G20630.1,AT4G20530.1,AT4G20670.1
import sys
i = 0 
for line in open(sys.argv[1]):
	if line[0] == '#':
		continue
	cell = line.strip().split('\t')
	if 'Tandem' in cell[1]:
		i += len(cell[-1].split(','))
print(i)	
