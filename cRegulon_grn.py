import numpy as np
import pandas as pd
from ismember import ismember
import scipy.sparse as sparse
import os
import scipy.io as scio
from numpy.matlib import repmat
import numpy_groupies as npg
from collections import Counter
import mpmath
import subprocess
from multiprocessing import Process
import pybedtools
import shutil
import argparse
import sys

parser = argparse.ArgumentParser(description='GRN construction for each cell cluster')
parser.add_argument('--name','-n',type=str, default = "Run",required=True,help="Task name")
parser.add_argument('--celltype','-ct',type=str,default = "",required=True,help='Cell type name')
parser.add_argument('--genome','-g',type=str, default = "",required=True,help='Genome build')
parser.add_argument('--cores','-p',type=int,default = 4,required=True,help='Path to the scATAC-seq data file')

args = parser.parse_args()

def mfbs(MCFile, TFName, Element_name, motifName, motifWeight, Match2):
	MC = pd.read_csv(MCFile, sep='\t', header=None)
	Mf3 = MC.iloc[:,2]
	Md1, Mf1 = ismember(MC.iloc[:,0], Element_name)
	Md2, Mf2 = ismember(MC.iloc[:,1], pd.DataFrame(motifName))
	Mf1_2 = np.zeros(shape=(MC.shape[0],1))
	Mf2_2 = np.zeros(shape=(MC.shape[0],1))
	Mf1_2[np.where(Md1 == True)[0], 0] = Mf1
	Mf2_2[np.where(Md2 == True)[0], 0] = Mf2
	Mt1 = np.setdiff1d(np.arange(len(motifName)), np.unique(Mf2))
	Mf2 = np.concatenate((np.squeeze(Mf2_2[(Md1 & Md2)]), Mt1 ))
	Mf1 = np.concatenate((np.squeeze(Mf1_2[(Md1 & Md2)]), np.ones(len(Mt1))))
	Mf3 = np.concatenate((Mf3[np.where(Md1 & Md2)[0]], np.zeros(len(Mt1))))
	Mt1 = np.setdiff1d(np.arange(len(Element_name)), np.unique(Mf1))
	Mf1 = np.concatenate((Mf1, Mt1))
	Mf2 = np.concatenate((Mf2, np.ones(len(Mt1))))
	Mf3 = np.concatenate((Mf3, np.zeros(len(Mt1))))
	MMotif_binding = sparse.csr_matrix((Mf3, (Mf2, Mf1)), shape=(len(motifName), len(Element_name)))
	MMotif_binding = sparse.diags(np.squeeze(np.asarray(1/(motifWeight + 0.1)))) * MMotif_binding
	MMotif_binding = np.log(MMotif_binding.todense()+1)
	MTF_binding = np.zeros((len(TFName), len(Element_name)))
	Mf1_2 = np.zeros(shape=(Match2.shape[0], 1))-1
	Mf2_2 = np.zeros(shape=(Match2.shape[0], 1))-1
	Md1, Mf1 = ismember(pd.DataFrame(Match2[:, 0]), motifName)
	Md2, Mf2 = ismember(pd.DataFrame(Match2[:, 1]), TFName)
	Mf1_2[np.where(Md1 == True)[0], 0] = Mf1
	Mf2_2[np.where(Md2 == True)[0], 0] = Mf2
	Mf1 = np.squeeze(Mf1_2)
	Mf2 = np.squeeze(Mf2_2)
	for i in range(len(TFName)):
		Ma = np.where(Mf2 == i)[0]
		if len(Ma) > 1:
			MTF_binding[i, :] = MMotif_binding[Mf1[Ma].astype(np.int64), :].max(axis=0)
		elif len(Ma) == 1:
			MTF_binding[i, :] = MMotif_binding[Mf1[Ma].astype(np.int64), :]
		else:
			MTF_binding[i, :] = np.zeros((len(Element_name),))
	MTF_binding = sparse.csr_matrix(MTF_binding)
	return MTF_binding

def mfbs_c(N,TFName, Element_name, motifName, motifWeight, Match2):
	TFB = sparse.csr_matrix(([], ([], [])), shape=(len(TFName), len(Element_name)))
	for i in range(N):
		TFB1 = mfbs(".MotifTarget"+str(i+1)+".txt",TFName, Element_name, motifName, motifWeight, Match2)
		TFB = TFB + TFB1
	return TFB

def mf_collect(genome):
	if "hg" in genome:
		species = "human"
	elif "mm" in genome:
		species = "mouse"
	filename = 'knownResults.txt'
	fileID = open(filename)
	C = np.loadtxt(fileID, dtype=str, skiprows=1)
	fileID.close()

	for i in range(C.shape[0]):
		C[i, 6] = float(C[i, 6].strip('%'))
		C[i, 8] = float(C[i, 8].strip('%'))
	Score = np.zeros(shape=(C.shape[0], 3))
	Score[:, 0] = np.array([float(-mpmath.log10(mpmath.mpf(C[i, 2]))) for i in range(C.shape[0])])
	Score[:, 1] = (C[:, 6].astype(float) + 0.1) / (C[:, 8].astype(float) + 0.1)
	Score[Score[:, 0] > 100, 0] = 100
	Score[:, 2] = np.sqrt(Score[:, 0] * Score[:, 1])

	f = np.argsort(Score[:, 2])[::-1]
	Score = Score[f, :]
	Name = C[f, 0]

	filename = 'knownResults_rank.txt'
	fid = open(filename, 'w')
	fid.write('Motif\t-log10(p)\tFoldChange\tScore\n')
	for i in range(Name.shape[0]):
		fid.write('{}\t{}\t{}\t{}\n'.format(Name[i], Score[i, 0], Score[i, 1], Score[i, 2]))
	fid.close()

	Motifmatch = scio.loadmat('../../../Data/MotifMatch_{}_rmdup.mat'.format(species))

	Match2 = np.empty((Motifmatch['Match2'].shape[0], 2), dtype='U100')
	for i in range(Motifmatch['Match2'].shape[0]):
		Match2[i, 0] = Motifmatch['Match2'][i, 0].item()
		Match2[i, 1] = Motifmatch['Match2'][i, 1].item()

	d, f = ismember(Match2[:, 0], Name)
	Score1 = Score[f, :]
	Match2 = Match2[np.where(d == True), :]
	Match2 = np.squeeze(Match2)
	TF, ic = np.unique(Match2[:, 1], return_inverse=True)
	TFScore = npg.aggregate(ic, Score1[:, 2], func='max')
	f = np.argsort(TFScore)[::-1]
	Results = np.column_stack((TF[f], TFScore[f]))
	filename = 'knownResults_TFrank.txt'
	fid = open(filename, 'w')
	for i in range(Results.shape[0]):
		fid.write('{}\t{}\n'.format(Results[i, 0], Results[i, 1]))
	fid.close()
    

def Split(N):
	f = open("region.txt")
	pn = f.readlines();f.close()
	NN = len(pn);NM = int(NN/N)+1
	for i in range(N):
		g = open(".region"+str(i+1)+'.txt','w')
		for j in range(len(pn)):
			if int(j/NM) == i:
				g.write(pn[j])
		g.close()

def Homer_prior(peak,genome,output):
	command="mkdir "+output+";findMotifsGenome.pl "+peak+" "+genome+" ./"+output+" -size given -find ../../Data/all_motif_rmdup -preparsedDir ../../Data/Homer/ > ./"+output+"/"+output+".bed;cat ./"+output+"/"+output+".bed | awk 'NR>1'|cut -f 1,4,6 > "+output+".txt;rm -rf "+output
	logfile=open(output+".log",'w')
	ph=subprocess.Popen(command,shell=True,stderr=logfile)
	return_code=ph.wait()

def Homer(peak,genome,output):
	command="mkdir "+output+";findMotifsGenome.pl "+peak+" "+genome+" ./"+output+" -size given -find ../../Data/all_motif_rmdup > ./" + output + "/" + output + ".bed;cat ./" + output + "/" + output + ".bed | awk 'NR>1'|cut -f 1,4,6 > " + output + ".txt;rm -rf " + output
	logfile=open(output+".log",'w')
	ph = subprocess.Popen(command, shell=True,stderr=logfile)
	return_code = ph.wait()

def MotifFind(genome, num_processes, prior):
	"""
	Args:
		genome (str): 
		num_processes (int): 线程数
		prior (int): 是否有已经扫好的homer motif文件
	"""
	processes = []
	Split(num_processes)
	if prior == 1:
		for p in range(num_processes):
			process = Process(target=Homer_prior, args=(".region"+str(p+1)+'.txt',str(genome),".MotifTarget"+str(p+1)))
			processes.append(process)
	else:
		for p in range(num_processes):
			process = Process(target=Homer, args=(".region"+str(p+1)+'.txt',str(genome),".MotifTarget"+str(p+1)))
			processes.append(process)
	for process in processes:
		process.start()
	for process in processes:
		process.join()
        
def GRN(name, celltype, genome,num_processes):
	if "hg" in genome:
		species = "human"
	elif "mm" in genome:
		species = "mouse"
	C = pd.read_csv('openness2.bed', sep='\t', header=None)
	Element_name = C.iloc[:, 0]
	Opn = C.iloc[:, 1]
	Opn_median = C.iloc[:, 2]

	# 读取mat文件
	# Match2, motifName, motifWeight
	MotifMatch_mouse_rmdup = scio.loadmat('../../Data/MotifMatch_{}_rmdup.mat'.format(species))
	# Exp_median, List, R2, TFExp_median, TFName
	TFTG_corr_mouse = scio.loadmat('../../Data/TFTG_corr_{}.mat'.format(species))

	Match2 = np.empty((MotifMatch_mouse_rmdup['Match2'].shape[0], 2), dtype='U100')
	for i in range(MotifMatch_mouse_rmdup['Match2'].shape[0]):
		Match2[i, 0] = MotifMatch_mouse_rmdup['Match2'][i, 0].item()
		Match2[i, 1] = MotifMatch_mouse_rmdup['Match2'][i, 1].item()

	motifName = np.empty((MotifMatch_mouse_rmdup['motifName'].shape[0], 1), dtype='U100')
	for i in range(MotifMatch_mouse_rmdup['motifName'].shape[0]):
		motifName[i, 0] = MotifMatch_mouse_rmdup['motifName'][i, 0].item()

	motifWeight = MotifMatch_mouse_rmdup['motifWeight']

	TFName = np.empty((TFTG_corr_mouse['TFName'].shape[0], 1), dtype='U100')
	for i in range(TFTG_corr_mouse['TFName'].shape[0]):
		TFName[i, 0] = TFTG_corr_mouse['TFName'][i, 0].item()

	List = np.empty((TFTG_corr_mouse['List'].shape[0], 1), dtype='U100')
	for i in range(TFTG_corr_mouse['List'].shape[0]):
		List[i, 0] = TFTG_corr_mouse['List'][i, 0].item()

	R2 = TFTG_corr_mouse['R2']
	Exp_median = TFTG_corr_mouse['Exp_median']
	# ---------------------------

	N = num_processes
	TF_binding = mfbs_c(N,TFName, Element_name, motifName, motifWeight, Match2)

	# gene expr
	C = pd.read_csv("../../PseudoBulk/{}_{}_PS_RNA.txt".format(name,celltype), sep='\t', header=None)
	Symbol = C.iloc[:, 0]
	G = C.iloc[:, 1]

	alhfa = 0.5
	Opn_median = np.log2(1 + Opn_median)
	Opn1 = np.log2(1 + Opn)
	Opn = Opn1 * (Opn1 / (Opn_median + 0.5))

	geneName = np.intersect1d(List, Symbol)
	d, f = ismember(geneName, List)
	R2 = R2[:, f]
	Exp_median = Exp_median[f]

	d, f = ismember(geneName, Symbol)
	G = G[f]
	d1 = sorted(G)
	f1 = np.argsort(np.array(G))
	d2 = sorted(Exp_median)
	f2 = np.argsort(np.squeeze(Exp_median))
	G1 = np.empty(len(G))
	for i in range(len(G)):
		G1[f1[i]] = d2[i][0]
	G = np.multiply((np.power(G1, alhfa)), (G1 / (np.squeeze(Exp_median) + 0.5)))

	d, f = ismember(TFName, geneName)
	TFName = TFName[np.where(d == True)]
	TF_binding = TF_binding[np.where(d == True)[0], :]
	TFExp = G[f]
	R2 = R2[np.where(d == True)[0], :]

	C = pd.read_csv('./Enrichment/knownResults_TFrank.txt', header=None, sep='\t')
	d, f = ismember(TFName, C.iloc[:, 0])
	TF_motif = np.zeros(len(TFName))
	TF_motif[np.where(d == True)[0]] = C.iloc[:, 1][f]
	TFExp = TFExp * TF_motif

	C = pd.read_csv('peak_gene_100k_corr.bed', header=None, sep='\t')
	d, f = ismember(C.iloc[:, 0], Element_name)
	d1, f1 = ismember(C.iloc[:, 1], geneName)
	f_2 = np.zeros(shape=(C.shape[0], 1))
	f1_2 = np.zeros(shape=(C.shape[0], 1))
	f_2[np.where(d == True)[0], 0] = f
	f1_2[np.where(d1 == True)[0], 0] = f1

	f2, ia, ic = np.unique(np.hstack((f_2[(d & d1)], f1_2[(d & d1)])), axis=0, return_index=True, return_inverse=True)
	c3 = npg.aggregate(ic, C.iloc[d & d1, 2], func='min')
	c4 = npg.aggregate(ic, C.iloc[d & d1, 3], func='min')
	c4[np.where(c4 < 0.2)] = 0
	d0 = 500000
	c = np.exp(-1 * c3 / d0) * c4
	Opn[np.isnan(Opn)] = 0
	H1 = sparse.csr_matrix((c, (f2[:, 1], f2[:, 0])), shape=(len(geneName), len(Element_name)))
	TFO = np.multiply(TF_binding.todense(), np.tile(Opn.T, (np.size(TF_binding, 0), 1)))
	H1Tdense = H1.todense().T
	BOH = np.matmul(TFO, H1Tdense)
	Score = np.multiply(np.multiply((np.dot(TFExp.reshape(-1, 1), G.reshape(-1, 1).T)), (2 ** np.abs(R2))), BOH)
	Score[np.isnan(Score)] = 0
	np.savetxt('TFTG_regulationScore.txt', Score, delimiter='\t')
	np.savetxt('TFName.txt', TFName, fmt='%s', delimiter='\n')
	np.savetxt('TGName.txt', geneName, fmt='%s', delimiter='\n')

	TFTGControl = scio.loadmat('../../Data/TFTG_{}_nagetriveControl.mat'.format(species))
	Back_net = TFTGControl['Back_net']
	for i in range(Back_net.shape[0]):
		for j in range(Back_net.shape[1]):
			Back_net[i, j] = Back_net[i, j].item()
	d, f = ismember(Back_net[:, 0], TFName)
	d1, f1 = ismember(Back_net[:, 1], geneName)
	f_2 = np.zeros(shape=(Back_net.shape[0], 1))
	f1_2 = np.zeros(shape=(Back_net.shape[0], 1))
	f_2[np.where(d == True)[0], 0] = f
	f1_2[np.where(d1 == True)[0], 0] = f1
	f2 = np.hstack((f_2[(d & d1)], f1_2[(d & d1)]))
	Score_T_1col = Score.T.reshape(-1, 1)
	aa = (f2[:, 1]) * Score.shape[0] + f2[:, 0]
	Back_score = Score_T_1col[aa.astype(np.int64)].squeeze().T

	Cut = np.percentile(np.asarray(Back_score), 99)
	[b, a] = np.where(Score.T > Cut)
	c = np.where(Score_T_1col > Cut)[0]
	c1 = Score_T_1col[c]
	Net = np.column_stack((TFName[a], geneName[b]))
	a1 = np.sort(H1.T.todense(), axis=0)[::-1]
	a2 = np.argsort(H1.T.todense(), axis=0)[::-1]
	a1 = a1[:10, :]
	a2 = a2[:10, :]
	TFTG_RE = [';'.join(Element_name[np.asarray(a2[np.where((TFO[a[i], a2[:, b[i]]] > 0) & (a1[:, b[i]] > 0))[0], b[i]]).squeeze()]) for i in range(len(a))]
	for i in range(len(a)):
		kk = np.asarray(a2[np.where((TFO[a[i], a2[:, b[i]]] > 0) & (a1[:, b[i]] > 0))[0], b[i]]).squeeze()
		if kk.size == 1:
			TFTG_RE[i] = TFTG_RE[i].replace(';', '')

	d = np.sort(np.asarray(c1).squeeze())[::-1]
	f = np.argsort(np.asarray(c1).squeeze())[::-1]
	Net = np.column_stack((Net[f], np.asarray(d).squeeze(), np.asarray(TFTG_RE)[f]))
	filename = '{}_{}_network.txt'.format(name,celltype)
	with open(filename, 'wt') as fid:
		fid.write('\t'.join(['TF', 'TG', 'Score', 'FDR', 'REs']) + '\n')
		for i in range(Net.shape[0]):
			fid.write('\t'.join([Net[i, 0], Net[i, 1], str(Net[i, 2]), str((np.sum(Back_score > d[i]) + 1) / len(Back_score)),Net[i, 3]]) + '\n')

def network(name,celltype,self_genome, num_processes=20, prior=0):
	
	folder_name = "./Networks/"+name+"_"+celltype
	os.makedirs(folder_name, exist_ok=True)
	os.chdir('./{}'.format(folder_name))
	
	# openness, region 
	df = pd.read_csv("../../PseudoBulk/{}_{}_PS_ATAC.txt".format(name,celltype), sep = '\t', header=None, names=['col1', 'col2', 'col3'])
	df[['chr', 'start', 'end']] = df['col1'].str.split('_', expand=True)
	df = df[['chr', 'start', 'end', 'col2']]
	df.to_csv('openness1.bed', sep='\t', header=False, index=False)
	df = pd.read_csv('openness1.bed', sep='\t', header=None)
	sorted_df = df.sort_values(by=[0, 1])
	sorted_df.to_csv('openness1.bed', sep='\t', header=False, index=False)
	sorted_df[3] = sorted_df[0] + '_' + sorted_df[1].astype(str) + '_' + sorted_df[2].astype(str)
	sorted_df[[0, 1, 2, 3]].to_csv('region.txt', sep='\t', header=False, index=False)
	sorted_df[[0, 1, 2]].to_csv('region.bed', sep='\t', header=False, index=False)

	print("step1: motif binding...")
	#MotifFind(self_genome, num_processes, prior)
	
	print("step2: Opn files...")
	openness1 = pybedtools.BedTool('openness1.bed')
	opn_median = pybedtools.BedTool('../../Data/Opn_median_{}.bed'.format(self_genome))
	result = openness1.intersect(opn_median, wa=True, wb=True, sorted=True).cut([0, 1, 2, 3, 7]).to_dataframe()
	result.iloc[:, 0] = result.iloc[:, 0].astype(str) + '_' + result.iloc[:, 1].astype(str) + '_' + result.iloc[:,2].astype(str) + '_' + result.iloc[:, 3].astype(str)
	result = result.iloc[:, [0, 4]]
	max_values = result.groupby('chrom')['score'].max().reset_index()
	max_values2 = max_values.iloc[:, 0].str.rsplit('_', n=1, expand=True)
	max_values2['score'] = max_values.iloc[:, 1]
	max_values2.to_csv("openness2.bed", sep='\t', header=False, index=False)

	os.mkdir("Enrichment")
	data = pd.read_csv("openness2.bed", sep='\t', header=None)
	data['new_col'] = (data[1] + 0.5) / (data[2] + 0.5)
	sorted_data = data.sort_values('new_col', ascending=False)
	filtered_data = sorted_data[sorted_data.iloc[:, 3] - sorted_data.iloc[:, 2] < 2000].head(10000)
	filtered_data[['chr', 'start', 'end']] = filtered_data.iloc[:, 0].str.split('_', expand=True)
	filtered_data[['chr', 'start', 'end']].to_csv("./Enrichment/region.bed", sep='\t', header=False, index=False)

	print("step3: Motif collect...")
	os.chdir("./Enrichment/")
	if prior==1:
	    command = "findMotifsGenome.pl region.bed {} ./. -size given -mask -nomotif -mknown ../../../Data/all_motif_rmdup -preparsedDir ../../../Data/Homer/ -p 30".format(self_genome)
	    logfile=open("homer.log",'w')
	    subprocess.run(command, shell=True, stderr=logfile)
	else:
	    command = "findMotifsGenome.pl region.bed {} ./. -size given -mask -nomotif -mknown ../../../Data/all_motif_rmdup -p 30".format(self_genome)
	    logfile=open("homer.log",'w')
	    subprocess.run(command, shell=True, stderr=logfile)
	mf_collect(self_genome)

	print("step4: Prior...")
	os.chdir('../')
	a = pybedtools.BedTool('region.bed')
	b = pybedtools.BedTool('../../Data/RE_gene_corr_{}.bed'.format(self_genome))
	result = a.intersect(b, wa=True, wb=True, sorted=True).cut([0, 1, 2, 6, 7, 8]).to_dataframe()
	result.iloc[:, 0] = result.iloc[:, 0].astype(str) + '_' + result.iloc[:, 1].astype(str) + '_' + result.iloc[:,2].astype(str)
	result = result.iloc[:, [0, 3, 4, 5]]
	result.to_csv("peak_gene_100k_corr.bed", sep='\t', header=False, index=False)

	b = pybedtools.BedTool('../../Data/Enhancer_RE_gene_corr_{}.bed'.format(self_genome))
	result = a.intersect(b, wa=True, wb=True, sorted=True).cut([0, 1, 2, 6, 7, 8]).to_dataframe()
	result.iloc[:, 0] = result.iloc[:, 0].astype(str) + '_' + result.iloc[:, 1].astype(str) + '_' + result.iloc[:,2].astype(str)
	result = result.iloc[:, [0, 3, 4, 5]]
	result.to_csv("peak_gene_100k_corr.bed", mode='a', sep='\t', header=False, index=False)
	
	print("step5: Network...")
	GRN(name,celltype,self_genome,num_processes)
	
	# 清理中间文件
	
	for filename in os.listdir():
	    if filename.startswith('.'):
	        os.remove(filename)
	shutil.rmtree('Enrichment')
	os.remove("region.bed");os.remove("region.txt");os.remove("peak_gene_100k_corr.bed");os.remove("openness1.bed");os.remove("openness2.bed");
	print('GRN for'+celltype+' done!')

if __name__ == '__main__':
	network(name=args.name,celltype=args.celltype,self_genome=args.genome,num_processes=args.cores, prior=1)
