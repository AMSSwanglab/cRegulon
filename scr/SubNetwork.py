import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

f = open('./Input/Sample_TGName.txt')
TG = f.readlines();f.close()
TG = [TG[i].strip('\n') for i in range(len(TG))]

f = open('./Sample/Sample_network.txt')
net = f.readlines();f.close()
del net[0]
for i in range(len(net)):
    net[i] = net[i].strip('\n').split('\t')
    net[i][4] = net[i][4].split(';')

W1 = np.loadtxt("./Results/Sample/Sample_W1.txt").T
W2 = np.loadtxt("./Results/Sample/Sample_W2.txt").T
W0 = 0.9*W1 + 0.1*W2

for t in range(W0.shape[0]):
    W = W0[t]
    me = np.mean(W);sd = np.std(W)
    cTG = []
    for i in range(len(TG)):
        x = (W[i]-me)/sd
        p = 1-stats.norm.cdf(x)
        if p <= 0.05 or W[i]>=0.95:
            cTG.append(TG[i])
    f = open('./Results/Sample/Sample_TFs.txt')
    cTF = f.readlines();f.close()
    cTF = [cTF[i].split('\t')[0] for i in range(len(cTF))]
    g = open('./Results/Sample/Sample_SubNetwork.txt','w')
    cRE = []
    for i in range(len(net)):
        if net[i][0] in cTF and net[i][1] in cTG:
            g.write(net[i][0]+'\t'+net[i][1]+'\t'+net[i][2]+'\t')
            RE1 = []
            for j in range(len(net[i][4])):
                if net[i][4][j] != "":
                    RE1.append(net[i][4][j])
            cRE += RE1
            g.write(';'.join(RE1)+'\n')
    g.close()
