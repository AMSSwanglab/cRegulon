f = open('./data/Sample_scRNA.txt')
data = f.readlines();f.close()
g = open('Sample_GE.txt','w')
g.write(data[0]);del data[0]
gene = []
for i in range(len(data)):
	gene.append(data[i].split('\t')[0])

f = open('./scr/Mm10_GoodGene.bed')
a = f.readlines();f.close()
for i in range(len(a)):
	a[i] = a[i].strip('\n').split('\t')
	if a[i][4] in gene:
		indel = gene.index(a[i][4])
		g.write(data[indel])
g.close()

import numpy as np

f = open('Sample_GE.txt')
data = f.readlines();f.close()
g = open('Sample_GE.txt','w')
g.write(data[0]);del data[0]
gene = []
for i in range(len(data)):
	data[i] = data[i].split('\t')
	gene.append(data[i][0]);del data[i][0]
	for j in range(len(data[i])):
		data[i][j] = float(data[i][j])
data = np.array(data).T
for i in range(data.shape[0]):
	data[i] = data[i]/np.sum(data[i])*1000000
data = data.T
for i in range(len(gene)):
	g.write(gene[i])
	for j in range(data.shape[1]):
		g.write('\t'+str(data[i][j]))
	g.write('\n')
g.close()
