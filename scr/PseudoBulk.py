import numpy as np


f = open('./data/Sample_scRNA.txt')
data = f.readlines();f.close()
del data[0]
gene = [];
for i in range(len(data)):
	data[i] = data[i].split('\t')
	gene.append(data[i][0]);del data[i][0]
	for j in range(len(data[i])):
		data[i][j] = float(data[i][j])
	data[i] = np.sum(data[i])
SumC = np.sum(data)
g = open('./Sample/Sample.txt','w')
for i in range(len(gene)):
	CPM = data[i]/SumC*1000000
	g.write(gene[i]+'\t'+str(CPM)+'\n')
g.close()

f = open('./data/Sample_scATAC.txt')
data = f.readlines();f.close()
del data[0]
peak = []
for i in range(len(data)):
	data[i] = data[i].split('\t')
	peak.append(data[i][0]);del data[i][0]
	for j in range(len(data[i])):
		data[i][j] = float(data[i][j])
	data[i] = np.sum(data[i])
SumC = np.sum(data)
g = open('./Sample/Sample_ATAC.txt','w')
for i in range(len(peak)):
	CPM = data[i]/SumC*1000000
	g.write(peak[i]+'\t'+str(CPM)+'\n')
g.close()
