import os
import gzip
import numpy as np
from scipy.sparse import csr_matrix

import argparse
import sys
import os
import numpy as np

parser = argparse.ArgumentParser(description='Pseudo-bulk for each cell cluster')
parser.add_argument('--name','-n',type=str, default = "Run",required=False,help="Task name")
parser.add_argument('--rna','-r',type=str,default = "",required=True,help='Path to the scRNA-seq data file')
parser.add_argument('--rna_meta','-rm',type=str, default = "",required=True,help='Path to the scRNA-seq meta file')
parser.add_argument('--atac','-a',type=str,default = "",required=True,help='Path to the scATAC-seq data file')
parser.add_argument('--atac_meta','-am',type=str, default = "",required=True,help='Path to the scATAC-seq meta file')
parser.add_argument('--species','-g',type=str, default = "",required=True,help='The analyzing species')

args = parser.parse_args()
AllArgs = [args.rna,args.rna_meta,args.atac,args.atac_meta]
ArgsErr = ["ERROR! scRNA-seq data file missing","ERROR! scRNA-seq data cluster file missing","ERROR! scATAC-seq data file missing","ERROR! scATAC-seq data cluster file missing"]
ERRORS = 0
for i in range(4):
	if os.path.exists(AllArgs[i]) == False:
		print(ArgsErr[i])
		ERRORS += 1
if ERRORS > 0:
	sys.exit(1)

def Read10X_ATAC(path):
	try:
		with gzip.open(path+'/peaks.bed.gz', 'rt') as f:
			pk = f.readlines();f.close()
			pk = [pk[i].strip('\n').replace('\t','_') for i in range(len(pk))]
	except OSError:
		try:
			with open(path+'/peaks.bed', 'rt') as f:
				pk = f.readlines();f.close()
				pk = [pk[i].strip('\n').replace('\t','_') for i in range(len(pk))]
		except OSError:
			print("No correct peak file!")
	try:
		with gzip.open(path+'/barcodes.tsv.gz', 'rt') as f:
			bc = f.readlines();f.close()
			bc = [bc[i].strip('\n') for i in range(len(bc))]
	except OSError:
		try:
			with open(path+'/barcodes.tsv', 'rt') as f:
				bc = f.readlines();f.close()
				bc = [bc[i].strip('\n') for i in range(len(bc))]
		except OSError:
			print("No correct barcode file!")
	try:
		with gzip.open(path+'/matrix.mtx.gz', 'rt') as f:
			count = np.loadtxt(f, skiprows=2).astype('int')
	except OSError:
		try:
			with open(path+'/matrix.mtx', 'rt') as f:
				count = np.loadtxt(f, skiprows=2).astype('int')
		except OSError:
			print("No correct count file!")

	if count[0][1] != len(bc) or count[0][0] != len(pk):
		print("Error!!!")
	else:
		count = np.delete(count, 0, axis=0).T
		count = csr_matrix((count[2], (count[0]-1, count[1]-1)), shape=(len(pk), len(bc)))
		return pk, bc, count

#peak,cell,data = Read10X_ATAC("./scATAC/")

def Read10X_RNA(path):
	try:
		with gzip.open(path+'/features.tsv.gz', 'rt') as f:
			gs = f.readlines();f.close()
			gs = [gs[i].strip('\n').split('\t')[1] for i in range(len(gs))]
	except OSError:
		try:
			with open(path+'/features.tsv', 'rt') as f:
				gs = f.readlines();f.close()
				gs = [gs[i].strip('\n').split('\t')[1] for i in range(len(gs))]
		except OSError:
			print("No correct peak file!")
	try:
		with gzip.open(path+'/barcodes.tsv.gz', 'rt') as f:
			bc = f.readlines();f.close()
			bc = [bc[i].strip('\n') for i in range(len(bc))]
	except OSError:
		try:
			with open(path+'/barcodes.tsv', 'rt') as f:
				bc = f.readlines();f.close()
				bc = [bc[i].strip('\n') for i in range(len(bc))]
		except OSError:
			print("No correct barcode file!")
	try:
		with gzip.open(path+'/matrix.mtx.gz', 'rt') as f:
			count = np.loadtxt(f, skiprows=2).astype('int')
	except OSError:
		try:
			with open(path+'/matrix.mtx', 'rt') as f:
				count = np.loadtxt(f, skiprows=2).astype('int')
		except OSError:
			print("No correct count file!")

	if count[0][1] != len(bc) or count[0][0] != len(gs):
		print("Error!!!")
	else:
		count = np.delete(count, 0, axis=0).T
		count = csr_matrix((count[2], (count[0]-1, count[1]-1)), shape=(len(gs), len(bc)))
		return gs, bc, count

def PseudoBulk(name,rna,rna_meta,atac,atac_meta):
	folder = "./PseudoBulk/"
	f = open(rna_meta)
	meta = f.readlines();f.close()
	cmeta = [meta[i].strip('\n').split('\t')[0] for i in range(len(meta))]
	tmeta = [meta[i].strip('\n').split('\t')[1] for i in range(len(meta))]
	ct = list(set(tmeta));ct.sort()
	print("We are processing scRNA-seq data with "+str(len(ct))+" cell clusters...")
	gene,cell,data = Read10X_RNA(rna)
	gc = open(folder+name+"_"+'CellType.txt','w')
	for i in range(len(ct)):
		gc.write(name+"_"+ct[i]+'\n')
		Indel = [cell.index(cmeta[indel]) for indel,x in enumerate(tmeta) if x==ct[i]]
		datac = data[:,Indel]
		datac = np.sum(datac,1)/np.sum(datac)*1000000
		g = open(folder+name+"_"+ct[i]+'_PS_RNA.txt','w')
		for j in range(len(gene)):
			g.write(gene[j]+'\t'+str(datac[j,0])+'\n')
		g.close()
	gc.close()

	f = open(atac_meta)
	meta = f.readlines();f.close()
	cmeta = [meta[i].strip('\n').split('\t')[0] for i in range(len(meta))]
	tmeta = [meta[i].strip('\n').split('\t')[1] for i in range(len(meta))]
	ct = list(set(tmeta));ct.sort()
	print("We are processing scATAC-seq data with "+str(len(ct))+" cell clusters...")
	peak,cell,data = Read10X_ATAC(atac)
	for i in range(len(ct)):
		Indel = [cell.index(cmeta[indel]) for indel,x in enumerate(tmeta) if x==ct[i]]
		datac = data[:,Indel]
		datac = np.sum(datac,1)/np.sum(datac)*1000000
		g = open(folder+name+"_"+ct[i]+'_PS_ATAC.txt','w')
		for j in range(len(peak)):
			if datac[j,0] >= 2 and "GL" not in peak[j] and "JH" not in peak[j]:
				g.write(peak[j].replace(':','_').replace('-','_')+'\t'+str(datac[j,0])+'\n')

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

def TFAS(name,rna,rna_meta,genome):
	f = open('./PseudoBulk/'+name+'_CellType.txt')
	ct = f.readlines();f.close()
	ct = [ct[c].strip('\n') for c in range(len(ct))]
	f = open("./Data/AllTFs_"+genome+".txt")
	TF = f.readlines();f.close()
	TF = [TF[i].strip('\n') for i in range(len(TF))]

	f = open(rna_meta)
	meta = f.readlines();f.close()
	cmeta = [meta[i].strip('\n').split('\t')[0] for i in range(len(meta))]
	tmeta = [name+"_"+meta[i].strip('\n').split('\t')[1] for i in range(len(meta))]

	gene,cell,data = Read10X_RNA(rna)
	cindex = [cell.index(cmeta[i]) for i in range(len(cmeta))];data = data[:,cindex]
	Indel_ct = [[] for c in range(len(ct))]
	for i in range(len(cmeta)):
		Indel_ct[ct.index(tmeta[i])].append(i)

	AllTFAS = []
	for i in range(len(TF)):
		AllTFAS.append([])
		indel = gene.index(TF[i])
		for c in range(len(ct)):
			e1 = np.mean(data[indel][:,Indel_ct[c]])
			e2 = np.mean(np.delete(data[indel].toarray(),Indel_ct[c]))
			AllTFAS[i].append(TFC(e1,e2))
	g = open("./PseudoBulk/"+name+"_TFES.txt",'w')
	g.write(name+"\t"+"\t".join(ct)+'\n')
	for i in range(len(TF)):
		g.write(TF[i])
		for c in range(len(ct)):
			g.write("\t"+str(AllTFAS[i][c]))
		g.write("\n")
	g.close()

if __name__ == '__main__':
	PseudoBulk(args.name,args.rna,args.rna_meta,args.atac,args.atac_meta)
	TFAS(args.name,args.rna,args.rna_meta,args.species)
