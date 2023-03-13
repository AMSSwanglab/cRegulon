import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

f = open('./Input/Sample_TFName.txt')
TF = f.readlines();f.close()
TF = [TF[i].strip('\n') for i in range(len(TF))]

X = np.loadtxt('./Results/Sample/X.txt')
factor1 = [1.5,0.4,1.5,1.5]
factor2 = [1,1,1,1]
for i in range(4):
    X1 = X[:,i].reshape(X.shape[0],1)#/np.max(X[:,i])
    XXT = (np.dot(X1,X1.T))**(1/factor1[i])*factor2[i]
    sample = []
    for j in range(X.shape[0]):
        for k in range(j+1,X.shape[0]):
            #if XXT[j][k] > 0:
                sample.append(XXT[j][k])
    par = stats.gamma.fit(sample)
    g1 = open('./Results/Sample/Sample_TFPairs.txt','w')
    g2 = open('./Results/Sample/Sample_TFs.txt','w')
    N = [[],[]]
    for j in range(X.shape[0]):
        for k in range(j+1,X.shape[0]):
            pval = stats.gamma.sf(XXT[j][k],par[0],par[1],par[2])
            if pval <= 0.01:
                g1.write(TF[j]+'\t'+TF[k]+'\t'+str(XXT[j][k])+'\t'+str(pval)+'\n')
                N[0].append(TF[j]+'\t'+str(X[j,i]));N[0].append(TF[k]+'\t'+str(X[k,i]))
    N[0] = list(set(N[0]))
    for j in range(len(N[0])):
        N[1].append(float(N[0][j].split('\t')[1]))
    for j in range(len(N[1])):
        indel = N[1].index(max(N[1]))
        g2.write(N[0][indel]+'\n')
        N[1][indel] = -100
    g1.close();g2.close()
