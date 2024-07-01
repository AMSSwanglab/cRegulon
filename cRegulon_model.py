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
import warnings
from sklearn.exceptions import ConvergenceWarning
warnings.filterwarnings("ignore", category=ConvergenceWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
import os
import argparse
import sys

parser = argparse.ArgumentParser(description='cRegulon inference from all cell cluster')
parser.add_argument('--name','-n',type=str, default = "Run",required=True,help="Task name")
parser.add_argument('--module_number','-mn',type=int,default = -1,required=False,help='The number of TF modules')
parser.add_argument('--module_max','-mmax',type=int,default = -1,required=False,help='The maximum number of TF modules')
parser.add_argument('--module_min','-mmin',type=int,default = -1,required=False,help='The minimum number of TF modules')

args = parser.parse_args()

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
	return R3, csi

def find_intersection(lists):
	intersection_set = set(lists[0])
	for lst in lists[1:]:
		intersection_set &= set(lst)
	return list(intersection_set)

def CSI(name):
	f = open('./PseudoBulk/'+name+'_CellType.txt')
	ct = f.readlines();f.close()
	AllTF = [];AllTG = [];AllTRS = []
	for c in range(len(ct)):
		ct[c] = ct[c].strip('\n')
		AllTRS.append(np.loadtxt('./Networks/'+ct[c]+'/TFTG_regulationScore.txt'))
		f = open("./Networks/"+ct[c]+"/TFName.txt")
		a = f.readlines();f.close();a = [a[i].strip('\n') for i in range(len(a))]
		AllTF.append(a)
		f = open("./Networks/"+ct[c]+"/TGName.txt")
		a = f.readlines();f.close();a = [a[i].strip('\n') for i in range(len(a))]
		AllTG.append(a)

	TF0 = find_intersection(AllTF);TG0 = find_intersection(AllTG);TF0.sort();TG0.sort()
	TFI = [[AllTF[c].index(TF0[i]) for i in range(len(TF0))] for c in range(len(ct))]
	TGI = [[AllTG[c].index(TG0[i]) for i in range(len(TG0))] for c in range(len(ct))]
	AllTRS = [AllTRS[c][TFI[c],:][:,TGI[c]] for c in range(len(ct))]

	TFSTD = np.array([np.std(AllTRS[c],1) for c in range(len(ct))]).max(0)
	TGSTD = np.array([np.std(AllTRS[c],0) for c in range(len(ct))]).max(0)
	TF1 = [];TFI = []
	for i in range(len(TF0)):
		if TFSTD[i] > 0:
			TF1.append(TF0[i])
			TFI.append(i)
	TG1 = [];TGI = []
	g = open("./CSI/TGName.txt","w")
	for i in range(len(TG0)):
		if TGSTD[i] > 0:
			TG1.append(TG0[i]);g.write(TG0[i]+"\n")
			TGI.append(i)
	g.close()
	AllTRS = [AllTRS[c][TFI,:][:,TGI] for c in range(len(ct))]

	STD = np.array([np.std(connection_specificity_index(AllTRS[c])[1],1) for c in range(len(ct))]).max(0)
	TF2 = [];TFI = []
	g = open("./CSI/TFName.txt",'w')
	for i in range(len(TF1)):
		if STD[i] > 0:
			TF2.append(TF1[i]);g.write(TF1[i]+'\n')
			TFI.append(i)
	g.close()
	AllTRS = [AllTRS[c][TFI,:] for c in range(len(ct))]
    
	AllCSI = []
	for c in range(len(ct)):
		TRS, CSI = connection_specificity_index(AllTRS[c])
		np.savetxt('./CSI/'+ct[c]+'_CSI.txt',CSI,delimiter='\t',fmt='%1.8f')
		np.savetxt('./CSI/'+ct[c]+'_TRS.txt',TRS,delimiter='\t',fmt='%1.8f')
		AllCSI.append(CSI)
	return AllCSI, TF2, TG1

def TFAS(Name,TF):
	f = open('./PseudoBulk/'+Name+'_CellType.txt')
	ct = f.readlines();f.close()
	ct = [ct[c].strip('\n') for c in range(len(ct))]
	
	f = open('./PseudoBulk/'+Name+'_TFES.txt')
	TFAS0 = f.readlines();f.close()
	del TFAS0[0]
	TF1 = []
	for i in range(len(TFAS0)):
		TFAS0[i] = TFAS0[i].split('\t')
		TF1.append(TFAS0[i][0]);del TFAS0[i][0]
	TFAS0 = np.array(TFAS0)

	AllTFAS = []
	for i in range(len(TF)):
		indel = TF1.index(TF[i])
		AllTFAS.append(TFAS0[indel])
	return np.array(AllTFAS).T.astype('float')

def CSIS(CS,TS,Name):
	f = open('./PseudoBulk/'+Name+'_CellType.txt')
	ct = f.readlines();f.close()
	ct = [ct[c].strip('\n') for c in range(len(ct))]
	CSTS = []
	for c in range(len(CS)):
		TFES = np.dot(TS[c].reshape((TS[c].shape[0],1)),TS[c].reshape((1,TS[c].shape[0])))
		csist = (CS[c]*TFES)**1/3
		CSTS.append(csist)
		np.savetxt('./CSI/'+ct[c]+'_CSIS.txt',csist,delimiter='\t',fmt='%1.8f')
	return CSTS


from sklearn.decomposition import NMF
from numpy import linalg as LA
import warnings
from sklearn.exceptions import ConvergenceWarning
warnings.filterwarnings("ignore", category=ConvergenceWarning)
import os

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

def RunModuleK(C,K):
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
		loss1 = loss2;
		X = XNext.copy();L = LNext.copy()
	return loss1, X, L

def TFPairs(X, TF, Name):
	X = X.T
	cTF = []
	for kk in range(X.shape[0]):
		XK = X[kk]
		XT = [];TFT = []
		for i in range(80):
			indel = np.argmax(XK)
			XT.append(XK[indel]);TFT.append(TF[indel])
			XK[indel] = -1000
		XT = np.array(XT);XT = XT.reshape(XT.shape[0],1)
		XXT = np.dot(XT,XT.T)
		sample = []
		for j in range(XXT.shape[0]):
			for k in range(j+1,XXT.shape[0]):
				sample.append(XXT[j][k])
		par = stats.gamma.fit(sample)
		g = open("./Results/"+Name+'/cRegulon'+str(kk+1)+'_TFModule.txt','w')
		cTF.append([])
		for j in range(XXT.shape[0]):
			for k in range(j+1,XXT.shape[0]):
				pval = stats.gamma.sf(XXT[j][k],par[0],par[1],par[2])
				if pval <= 0.05:
					g.write(TFT[j]+'\t'+TFT[k]+'\t'+str(XXT[j][k])+'\t'+str(pval)+'\n')
					cTF[kk].append(TFT[j]);cTF[kk].append(TFT[k])
		cTF[kk] = list(set(cTF[kk]))
		g.close();
	return cTF
        
def WriteXL(X, L, TF, TG, Name, cutoff=0.1):
	f = open('./PseudoBulk/'+Name+'_CellType.txt')
	ct = f.readlines();f.close()
	ct = [ct[c].strip('\n') for c in range(len(ct))]
	TRS = [np.loadtxt('./CSI/'+ct[c]+'_TRS.txt') for c in range(len(ct))]
    
	os.makedirs('./Results/'+Name, exist_ok=True)
	g = open('./Results/'+Name+'/X.txt','w')
	g.write('TF\t'+'\t'.join(["M"+str(c+1) for c in range(X.shape[1])])+'\n')
	for i in range(len(TF)):
		g.write(TF[i])
		for j in range(X.shape[1]):
			g.write('\t'+str(round(X[i][j],8)))
		g.write('\n')
	g.close()
	cTFs = TFPairs(X, TF, Name)

	g = open('./Results/'+Name+'/A.txt','w')
	g.write('CellType\t'+'\t'.join(["M"+str(c+1) for c in range(X.shape[1])])+'\n')
	for c in range(len(L)):
		L[c] = np.diag(L[c])
		g.write(ct[c])
		for j in range(X.shape[1]):
			g.write('\t'+str(round(L[c][j],8)))
		g.write('\n')
	g.close()

	os.makedirs('./Results/'+Name+'/Annotation/', exist_ok=True)
	for c in range(len(ct)):
		f = open('./Networks/'+ct[c]+'/'+ct[c]+'_network.txt')
		net = f.readlines();f.close();del net[0]
		for i in range(len(net)):
			net[i] = net[i].strip('\n').split('\t')
			net[i][4] = net[i][4].split(';')
		TGX = np.dot((X.dot(np.diag(L[c]))).T,TRS[c])
		TGX = (TGX-np.min(TGX))/(np.max(TGX)-np.min(TGX))
		for j in range(X.shape[1]):
			if L[c][j] >= cutoff:
				W = TGX[j];me = np.mean(W);sd = np.std(W)
				cTG  = []
				for k in range(len(TG)):
					x = (W[k]-me)/sd
					p = 1-stats.norm.cdf(x)
					if p <= 0.05 or W[k]>=0.95:
						cTG.append(TG[k])    
				g = open('./Results/'+Name+'/Annotation/'+ct[c]+'_M'+str(j+1)+"_subnetwork.txt",'w')
				for k in range(len(net)):
					if net[k][0] in cTFs[j] and net[k][1] in cTG:
						RE = []
						for l in range(len(net[k][4])):
							if net[k][4][l] != "":
								RE.append(net[k][4][l])
						if len(RE) > 0:
							RE = ';'.join(RE)
							g.write(net[k][0]+'\t'+net[k][1]+'\t'+net[k][2]+'\t'+RE+'\n')
				g.close()
	for i in range(X.shape[1]):
		D1 = []
		for j in range(X.shape[0]):
			if X[j,i] > 0.05:
				D1.append(X[j,i])
		me = np.mean(D1);sd = np.std(D1)
		for j in range(X.shape[0]):
			X[j,i] = (X[j,i]-me)/sd
	g = open('./Results/'+Name+'/XN.txt','w')
	g.write('TF\t'+'\t'.join(["M"+str(c+1) for c in range(X.shape[1])])+'\n')
	for i in range(len(TF)):
		g.write(TF[i])
		for j in range(X.shape[1]):
			g.write('\t'+str(X[i][j]))
		g.write('\n')
	g.close()

def fit_line_through_points(points):
    x_coords = np.array([point[0] for point in points])
    y_coords = np.array([point[1] for point in points])
    A = np.vstack([x_coords, np.ones(len(x_coords))]).T
    m, b = np.linalg.lstsq(A, y_coords, rcond=None)[0]
    return m

def SelectK(C,MinK,MaxK):
	if MinK != MaxK:
		### Select best number of TF modules
		FinalLoss = [];FinalK = []
		for kk in range(MinK,MaxK+1):
			K = kk;print("Running",K,"TF modules...")
			Err, X, L = RunModuleK(C,K)
			FinalLoss.append(Err);FinalK.append(K)
		slope = []
		for i in range(MaxK-MinK):
			if i < MaxK-MinK-1:
				points = [[FinalK[i],FinalLoss[i]],[FinalK[i+1],FinalLoss[i+1]],[FinalK[i+2],FinalLoss[i+2]]]
				slope.append(fit_line_through_points(points))
			else:
				points = [[FinalK[i],FinalLoss[i]],[FinalK[i+1],FinalLoss[i+1]]]
				slope.append(fit_line_through_points(points))
		diff = []
		for i in range(MaxK-MinK):
			if i==0:
				diff.append(np.abs(slope[i]))
			else:
				diff.append((slope[i]-slope[i-1])*(float(slope[i]<0)*2-1)/100)
		KF = 0
		for i in range(MaxK-MinK-1):
			if diff[i+1]<1:
				KF = FinalK[i];break;
		return KF
	else:
		KF = MinK
		return MinK

if __name__ == '__main__':
	CSIF, TFF, TGF = CSI(args.name)
	AllTFAS = TFAS(args.name,TFF)
	CSISF = CSIS(CSIF,AllTFAS,args.name)
	if args.module_number != -1:
		lossF, XF, LF = RunModuleK(CSISF,args.module_number)
		WriteXL(XF, LF, TFF, TGF, args.name)
	elif args.module_max != -1 and args.module_min != -1:
		SK = SelectK(CSISF,args.module_min,args.module_max)
		lossF, XF, LF = RunModuleK(CSISF,SK)
		WriteXL(XF, LF, TFF, TGF, args.name)
	else:
		print("Please provide the TF module number or its range!")