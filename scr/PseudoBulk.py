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

def PseudoBulk(name,rna,rna_meta,atac,atac_meta):
	folder = "./PseudoBulk/"
	f = open(rna_meta)
	meta = f.readlines();f.close()
	meta = [meta[i].strip('\n').split('\t')[1] for i in range(len(meta))]
	ct = list(set(meta));ct.sort()
	print("We are processing scRNA-seq data with "+str(len(ct))+" cell clusters...")
	f = open(rna)
	data = f.readlines();f.close()
	del data[0]
	gene = [];
	for i in range(len(data)):
		data[i] = data[i].split('\t')
		gene.append(data[i][0]);del data[i][0]
	data = np.array(data).astype('float')
	gc = open(folder+name+"_"+'CellType.txt','w')
	for i in range(len(ct)):
		gc.write(name+"_"+ct[i]+'\n')
		Indel = [indel for indel,x in enumerate(meta) if x==ct[i]]
		datac = data[:,Indel]
		datac = np.sum(datac,1)/np.sum(datac)*1000000
		g = open(folder+name+"_"+ct[i]+'_PS_RNA.txt','w')
		for j in range(len(gene)):
			g.write(gene[j]+'\t'+str(datac[j])+'\n')
		g.close()
	gc.close()

	f = open(atac_meta)
	meta = f.readlines();f.close()
	meta = [meta[i].strip('\n').split('\t')[1] for i in range(len(meta))]
	ct = list(set(meta));ct.sort()
	print("We are processing scATAC-seq data with "+str(len(ct))+" cell clusters...")
	f = open(atac)
	data = f.readlines();f.close()
	del data[0]
	peak = [];
	for i in range(len(data)):
		data[i] = data[i].split('\t')
		peak.append(data[i][0]);del data[i][0]
	data = np.array(data).astype('float')
	AllPeaks = []
	for i in range(len(ct)):
		Indel = [indel for indel,x in enumerate(meta) if x==ct[i]]
		datac = data[:,Indel]
		datac = np.sum(datac,1)/np.sum(datac)*1000000
		g = open(folder+name+"_"+ct[i]+'_PS_ATAC.txt','w')
		for j in range(len(peak)):
			if datac[j] >= 2:
				AllPeaks.append(peak[j])
				g.write(peak[j].replace(':','_').replace('-','_')+'\t'+str(datac[j])+'\n')
	AllPeaks = list(set(AllPeaks))
	g = open(folder+name+"_"+'Peaks.bed','w')
	for i in range(len(AllPeaks)):
		AllPeaks[i] = AllPeaks[j].replace(':','_').replace('-','_')
		g.write(AllPeaks[i].replace('_','\t')+'\t'+AllPeaks[i]+'\n')
	g.close()
	
if __name__ == '__main__':
	PseudoBulk(args.name,args.rna,args.rna_meta,args.atac,args.atac_meta)
