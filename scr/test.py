#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 13 15:34:48 2023

@author: Zhanying
"""

import numpy as np
from scipy import stats
from sklearn.decomposition import NMF
from numpy import linalg as LA

def Filter1(T):
	TFSTD = [];TGSTD = []
	for i in range(len(T)):
		TFSTD.append(np.std(T[i],1))
		TGSTD.append(np.std(T[i],0))
	TFSTD = np.array(TFSTD);TFSTD = np.max(TFSTD,0);I1 = np.where(TFSTD > 0)[0]
	TGSTD = np.array(TGSTD);TGSTD = np.max(TGSTD,0);I2 = np.where(TGSTD > 0)[0]
	return I1, I2
def connection_specificity_index(R):
	R0 = np.log2(R+1)
	R1 = stats.zscore(R0,1);R1[np.isnan(R1)] = 0.0
	R2 = stats.zscore(R0,0);R2[np.isnan(R2)] = 0.0
	R3 = (R1+R2)/2;R3[R3<0] = 0
	P = np.corrcoef(R3);P[np.isnan(P)] = 0.0
	csi = np.zeros((R3.shape[0],R3.shape[0]))
	for i in range(csi.shape[0]):
		for j in range(i,csi.shape[0]):
			n = 0;cutoff = P[i][j] - 0.05
			csi[i][j] = np.sum((P[i]<cutoff)*(P[j]<cutoff))/csi.shape[0]
			csi[j][i] = csi[i][j]
	return csi

def Filter2(C, I1, I2):
	CSTD = []
	for i in range(len(C)):
		CSTD.append(np.std(C[i],1))
	CSTD = np.array(CSTD);CSTD = np.max(CSTD,0)
	I11 = I1[np.where(CSTD > 0)];I22 = I2.copy()
	return I11, I22

def CSI(Name):
	f = open('./PseudoBulk/'+Name+'_CellType.txt')
	ct = f.readlines();f.close()
	ct = [ct[c].strip('\n') for c in range(len(ct))]
	### Load TF and TG Names
	f = open('./Networks/'+ct[0]+'_TFName.txt')
	TF = f.readlines();f.close()
	TF = [TF[i].strip('\n') for i in range(len(TF))]
	f = open('./Networks/'+ct[0]+'_TGName.txt')
	TG = f.readlines();f.close()
	TG = [TG[i].strip('\n') for i in range(len(TG))]
	### Load TRS matrixes
	AllTRS = [];TFSTD = [];TGSTD = []
	for c in range(len(ct)):
		TRS = np.loadtxt('./Networks/'+ct[c]+'_TRS.txt')
		AllTRS.append(TRS)
		TFSTD.append(np.std(TRS,1))
		TGSTD.append(np.std(TRS,0))
	### First round filtering TFs and TGs
	IndelTF1, IndelTG1 = Filter1(AllTRS) ### First round filtering TFs and TGs
	### Second round filtering TFs and TGs
	AllCSI = [connection_specificity_index(AllTRS[c][IndelTF1,][:,IndelTG1]) for c in range(len(ct))]
	IndelTF2, IndelTG2 = Filter2(AllCSI,IndelTF1,IndelTG1)
	TFF = [TF[IndelTF2[i]] for i in range(IndelTF2.shape[0])]
	TGF = [TG[IndelTG2[i]] for i in range(IndelTG2.shape[0])]
	### Final CSI
	CSIF = []
	for c in range(len(ct)):
		CF = connection_specificity_index(AllTRS[c][IndelTF2,][:,IndelTG2])
		CSIF.append(CF)
	return CSIF, TFF, TGF

def TFC(x,y):
	eps = 1e-4
	if y > 0:
		fc = x/y-0.5
	else:
		fc = x/eps-0.5
	if fc > 1:
		fc = 1
	if fc < 0:
		fc = 0
	return fc

def TFAS(Name,TF,RNA,meta):
	f = open('./PseudoBulk/'+Name+'_CellType.txt')
	ct = f.readlines();f.close()
	ct = [ct[c].strip('\n') for c in range(len(ct))]
	f = open(rna_meta)
	meta = f.readlines();f.close()
	meta = [Name+"_"+meta[i].strip('\n').split('\t')[1] for i in range(len(meta))]
	Indel_ct = [[] for c in range(len(ct))]
	for i in range(len(meta)):
		Indel_ct[ct.index(meta[i])].append(i)
	f = open(rna)
	data = f.readlines();f.close()
	del data[0];gene = []
	for i in range(len(data)):
		data[i] = data[i].split('\t')
		gene.append(data[i][0]);del data[i][0]
	data = np.array(data).astype('float')
	AllTFAS = []
	for i in range(len(TF)):
		AllTFAS.append([])
		indel = gene.index(TF[i])
		for c in range(len(ct)):
			e1 = np.mean(data[indel][Indel_ct[c]])
			e2 = np.mean(np.delete(data[indel],Indel_ct[c]))
			AllTFAS[i].append(TFC(e1,e2))
	return np.array(AllTFAS).T

def CSIS(CS,TS):
	CSTS = []
	for c in range(len(CS)):
		TFES = np.dot(TS[c].reshape((TS[c].shape[0],1)),TS[c].reshape((1,TS[c].shape[0])))
		CSTS.append((CS[c]*TFES)**1/3)
	return CSTS

eps = 0.001
def Norm(X):
	X = X/np.sqrt(np.sum(X**2,0).reshape(1,X.shape[1]))
	X = X/(np.sum(X,1).reshape(X.shape[0],1)+eps)
	return X
def NormL(L):
	L = L/np.sum(L)
	return L
def Loss(C,x0,l0,mu1):
	t1 = 0;t2 = 0
	for c in range(len(C)):
		t1 += pow(LA.norm(C[c] - x0.dot(l0[c]).dot(x0.T), ord='fro'), 2)
		t2 += mu1*np.sum(C[c] * x0.dot(l0[c]).dot(x0.T))
	return t1 - t2

def RunModule(C,K):
	I = len(C)
	C0 = np.zeros(C[0].shape)
	for c in range(I):
		C0 += C[c]
	C0 = C0/I
	err1 = np.zeros(50)
	for i in range(50):
		model = NMF(n_components=K, init='random', random_state=i, solver='cd', max_iter=50)
		X0 = model.fit_transform(C0);X1 = model.components_
		err1[i] = LA.norm(C0 - np.dot(X0, X1), ord='fro')
	model = NMF(n_components=K, init='random', random_state=np.argmin(err1), solver='cd', max_iter=1000)
	X0 = model.fit_transform(C0);X1 = model.components_
	L0 = []
	for c in range(I):
		Li = np.diag([1/K for i in range(K)])
		L0.append(Li)
	L1 = 0;L2 = 0
	for c in range(I):
		L1 += pow(LA.norm(C[c] - X0.dot(L0[c]).dot(X0.T), ord='fro'), 2)
		L2 += np.sum(C[c] * (X0.dot(L0[c]).dot(X0.T)))
	mu1 = L2/L1
	err = 1000;
	X = X0.copy();L = L0.copy();loss1 = Loss(C,X,L,mu1)
	while err > 1e-4:
		XTX = X.T.dot(X)
		XX1 = np.zeros(X.shape);XX2 = np.zeros(X.shape)
		for c in range(I):
			XX1 += C[c].dot(X).dot(L[c])
			XX2 += X.dot(L[c]).dot(XTX).dot(L[c])
		XNext = X * ((1+0.5*mu1)*XX1)/(eps+XX2)
		XNext = Norm(XNext)
		LNext = []
		for c in range(I):
			LNext.append(L[c]*((1+mu1)*X.T.dot(C[c]).dot(X))/(eps+XTX.dot(L[c]).dot(XTX)))
			LNext[c] = NormL(LNext[c])
		loss2 = Loss(C,XNext,LNext,mu1)
		err = np.abs(loss2-loss1)
		loss1 = loss2
		X = XNext.copy();L = LNext.copy()
	return loss1, X, L

def TFModule(C,MinK,MaxK):
	if MinK != MaxK:
		### Select best number of TF modules
		FinalLoss = []
		for kk in range(MinK,MaxK+1):
			Err, X, L = RunModule(C,kk)
			FinalLoss.append(Err)
		K = np.argmin(FinalLoss)+MinK
		Err, X, L = RunModule(C,K)
	else:
		K = MinK
		Err, X, L = RunModule(C,K)
	return X, L

def WriteXL(X,L,TF,Name):
	np.savetxt('./cRegulon/'+Name+'_X.txt',X,delimiter='\t',fmt='%1.8f')
	for i in range(X.shape[1]):
		D1 = []
		for j in range(X.shape[0]):
			if X[j,i] > 0.05:
				D1.append(X[j,i])
		me = np.mean(D1);sd = np.std(D1)
		for j in range(X.shape[0]):
			X[j,i] = (X[j,i]-me)/sd
	g = open('./cRegulon/'+Name+'_XN.txt','w')
	g.write('TF\t'+'\t'.join(["M"+str(c+1) for c in range(X.shape[1])])+'\n')
	for i in range(len(TF)):
		g.write(TF[i])
		for j in range(X.shape[1]):
			g.write('\t'+str(X[i][j]))
		g.write('\n')
	g.close()
	LL = []
	for c in range(len(L)):
		LL.append(np.diag(L[c]))
	LL = np.array(LL)
	f = open('./PseudoBulk/'+Name+'_CellType.txt')
	ct = f.readlines();f.close()
	ct = [ct[c].strip('\n') for c in range(len(ct))]
	g = open('./cRegulon/'+Name+'_A.txt','w')
	g.write('CellType\t'+'\t'.join(["M"+str(c+1) for c in range(LL.shape[1])])+'\n')
	for c in range(LL.shape[0]):
		g.write(ct[c])
		for j in range(LL.shape[1]):
			g.write('\t'+str(LL[c][j]))
		g.write('\n')
	g.close()

if __name__ == '__main__':
	name = "CL";rna = "./Input/CL_scRNA.txt";rna_meta = "./Input/CL_scRNA_Cluster.txt";minK=7;maxK=7
	AllCSI,AllTF,AllTG = CSI(name)
	AllTFAS = TFAS(name,AllTF,rna,rna_meta)
	AllCSIS = CSIS(AllCSI,AllTFAS)
	FinalX, FinalL = TFModule(AllCSIS,minK,maxK)
	WriteXL(FinalX,FinalL,AllTF,name)
