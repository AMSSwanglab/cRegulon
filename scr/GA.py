import numpy as np

f = open('./data/Sample_scATAC.txt')
data = f.readlines();f.close()
g = open('Sample_GA.txt','w')
g.write('Cells\t'+data[0]);del data[0]
Peak = []
for i in range(len(data)):
	data[i] = data[i].split('\t')
	Peak.append(data[i][0])
	del data[i][0]
	for j in range(len(data[i])):
		data[i][j] = float(data[i][j])
data = np.array(data)

f = open('Sample_Gene_NearPeak.txt')
a = f.readlines();f.close()
f = open('Sample_Gene_NearPeakDis.txt')
b = f.readlines();f.close()
NearPeak = [[],[],[]]
for i in range(len(a)):
	a[i] = a[i].strip('\n').split('\t')
	a[i][0] = a[i][0].split('--')
	NearPeak[0].append(a[i][0][1])
	a[i][1] = a[i][1].split(',')
	b[i] = b[i].strip('\n').split('\t')[1].split(',')
	pk = [];ds = []
	for j in range(len(a[i][1])):
		pk.append(Peak.index(a[i][1][j]))
		ds.append(float(b[i][j]))
	NearPeak[1].append(pk);NearPeak[2].append(ds)

f = open('./scr/Mm10_GoodGene.bed')
a = f.readlines();f.close()
GA = [];gene = []
for i in range(len(a)):
	a[i] = a[i].strip('\n').split('\t')
	if a[i][4] in NearPeak[0]:
		gene.append(a[i][4])
		ga = []
		indel = NearPeak[0].index(a[i][4])
		ds = data[NearPeak[1][indel]];ds1 = (np.e)**(-1/(1+ds/2))
		ga.append(np.dot(NearPeak[2][indel],ds))
		GA.append(ga[0])
	else:
		gene.append(a[i][4])
		ga = [0.0 for j in range(data.shape[1])]
		GA.append(ga)
GA = np.array(GA)

for i in range(GA.shape[1]):
	GA[:,i] = GA[:,i] / np.sum(GA[:,i]) * 1000000
for i in range(GA.shape[0]):
	g.write(gene[i])
	for j in range(GA.shape[1]):
		g.write('\t'+str(GA[i][j]))
	g.write('\n')
g2.close()
