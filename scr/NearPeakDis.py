import numpy as np

def Dis(x,y):
        x = x.split('_');x[1] = int(x[1]);x[2] = int(x[2])
        y = y.split('_');y[1] = int(y[1]);y[2] = int(y[2])
        gx = list(range(int(x[1]),int(x[2])+1));gy = list(range(int(y[1]),int(y[2])+1));
        if len(set(gx)&set(gy)) > 0:
                return 0
        else:
                return max(x[1],y[1])-min(x[2],y[2])

scale = 20000/500

f = open('Sample_Gene_NearPeak.txt')
a = f.readlines();f.close()
g = open('Sample_Gene_NearPeakDis.txt','w')
for i in range(len(a)):
	a[i] = a[i].strip('\n').split('\t')
	a[i][0] = a[i][0].split('--')
	a[i][1] = a[i][1].split(',')
	ds = []
	for j in range(len(a[i][1])):
		ds.append(str(np.e**(-Dis(a[i][0][0],a[i][1][j])/scale)))
	g.write(a[i][0][1]+'\t'+','.join(ds)+'\n')
g.close()
