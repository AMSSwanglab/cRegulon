import numpy as np
from sklearn.decomposition import NMF
from numpy import linalg as LA
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
from sklearn.exceptions import ConvergenceWarning
warnings.filterwarnings("ignore", category=ConvergenceWarning)
import os
import argparse

parser = argparse.ArgumentParser(description='cRegulon annotation for only scRNA-seq data')
parser.add_argument('--name','-n',type=str, default = "Run",required=True,help="Task name")
parser.add_argument('--path_rna','-rna',type=str,default = '',required=True,help='The number of TF modules')
parser.add_argument('--module_number','-mn',type=int,default = 25,required=True,help='The maximum number of TF modules')
args = parser.parse_args()

def Annot(name,path_to_rna,K):
    f = open(path_to_rna)
    E = f.readlines();f.close()
    cell = E[0].strip("\n").split("\t");del cell[0];del E[0]
    Gene = []
    for i in range(len(E)):
        E[i] = E[i].strip('\n').split('\t')
        Gene.append(E[i][0]);del E[i][0]
    E = np.array(E).astype('float')

    f = open("./Data/HumanAtlasX.txt")
    C0 = f.readlines();f.close()
    MName = C0[0].strip('\n').split('\t');del MName[0];del C0[0]
    TF = []
    for i in range(len(C0)):
        C0[i] = C0[i].split('\t')
        TF.append(C0[i][0]);del C0[i][0]
    C0 = np.array(C0).astype('float')

    TFI = [[],[],[]]
    for i in range(len(TF)):
        TF[i] = TF[i].strip('\n')
        if TF[i] in Gene:
            TFI[0].append(TF[i])
            TFI[1].append(i)
            TFI[2].append(Gene.index(TF[i]))

    C0 = C0[TFI[1],]

    PS = E.sum(1)/np.sum(E)*1000000
    TG = [];TGI = []
    for i in range(len(Gene)):
        if PS[i]>=1:
            TGI.append(i);TG.append(Gene[i])
    E0 = E[TGI,];E1 = E[TFI[2],]

    E0 = np.log(1 + E0);E0 = E0 / np.max(E0)
    model = NMF(n_components=K, init='random', random_state=6, solver='cd', max_iter=100)
    W0 = model.fit_transform(E0)
    H0 = model.components_
    R0 = np.corrcoef(E1,E0)[0:E1.shape[0],E1.shape[0]:];R0[np.isnan(R0)] = 0.0;R0[R0<0] = 0
    X0 = C0[:,np.random.choice(range(C0.shape[1]),K)]
    A0 = np.ones((X0.shape[1],C0.shape[1]))
    A0 = A0/C0.shape[1]

    mu = pow(LA.norm(C0 - np.dot(X0, A0), ord='fro'), 2) / pow(LA.norm(E0 - R0.T.dot(X0).dot(H0), ord='fro'), 2)
    
    err = 1000;eps = 0.001;
    def Loss(x,r,h,a,e):
        t1 = 1/2 * pow(LA.norm(e - r.T.dot(x).dot(h),ord='fro'), 2)
        t2 = mu/2 * pow(LA.norm(C0 - x.dot(a),ord='fro'), 2)
        return t1 - t2
    def NormX(xx):
        xx = xx/np.sqrt(np.sum(xx**2,0).reshape(1,xx.shape[1]))
        xx = xx/(np.sum(xx,1).reshape(xx.shape[0],1)+eps)
        return xx
    def NormA(aa):
        aa = aa/(np.sum(aa,1).reshape(aa.shape[0],1)+eps)
        return aa


    X = X0.copy();R = R0.copy();H = H0.copy();A = A0.copy()
    loss1 = Loss(X,R,H,A,E0)
    epoch = 0
    while err > 1e-4:
        XNext = X * (R.dot(E0).dot(H.T)+mu*C0.dot(A.T)) / (eps + R.dot(R.T).dot(X).dot(H).dot(H.T)+mu*X.dot(A).dot(A.T))
        XNext = NormX(XNext)
        RNext = R * (X.dot(H).dot(E0.T)) / (eps + X.dot(H).dot(H.T).dot(X.T).dot(R))
        HNext = H * (X.T.dot(R).dot(E0)) / (eps + X.T.dot(R).dot(R.T).dot(X).dot(H))
        ANext = A * (X.T.dot(C0)) / (eps + X.T.dot(X).dot(A))
        ANext = NormA(ANext)

        M = np.zeros((K, K))
        for z in range(K):
            M[z, z] = LA.norm(HNext[z, :])
        HNext = np.dot(LA.inv(M), HNext)

        loss2 = Loss(XNext,RNext,HNext,ANext,E0)
        err = np.abs(loss2-loss1)
        loss1 = loss2
        X = XNext.copy();R = RNext.copy();H = HNext.copy();A = ANext.copy()
    
    os.makedirs('./Results/'+name, exist_ok=True)
    np.savetxt('./Results/'+name+'/X.txt',X,delimiter='\t')
    np.savetxt('./Results/'+name+'/R.txt',R,delimiter='\t')
    np.savetxt('./Results/'+name+'/H.txt',H,delimiter='\t')
    np.savetxt('./Results/'+name+'/A.txt',A,delimiter='\t')

if __name__ == '__main__':
    Annot(args.name, args.path_rna, args.module_number)
