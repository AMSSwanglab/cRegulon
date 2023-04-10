#/bin/bash

# Pseudo-bulk PECA2 v1.0.0 updated Mar 11th 2023 by Zhanying Feng
# step 1: Merge all cells into pseudo-bulk sample
# step 2: motif binding
# step 3: calculate opn
# step 4: corr+dist
# step 5: score(Exp,binding,Opn,weight)

input=$1
genome=$2

echo step 1: Pseudo-bulk generating....
mkdir ${input}
sed s/Sample/${input}/g ./scr/PseudoBulk.py > PseudoBulk.py
python3 PseudoBulk.py;rm PseudoBulk.py
cd ${input}/
cat ${input}_ATAC.txt | grep -v GL | grep -v JH | awk '$2>2' | tr '_' '\t' | sortBed | awk -v OFS='\t' '{print $1"_"$2"_"$3,$4}' > openness.bed
cat ${input}_ATAC.txt | grep -v GL | grep -v JH | awk '$2>2' | tr '_' '\t' | sortBed | awk -v OFS='\t' '{print $1,$2,$3}' > region.bed
cat ${input}_ATAC.txt | grep -v GL | grep -v JH | awk '$2>2' | tr '_' '\t' | sortBed | awk -v OFS='\t' '{print $1,$2,$3,$1"_"$2"_"$3}' > region.txt

if [ `echo $genome|grep hg|wc -l` -gt 0 ]
then
	species=hs
	speciesFull=human
else
	species=mm
	speciesFull=mouse
fi

echo step 2: motif binding....
findMotifsGenome.pl region.txt ${genome} ./. -p 16 -size given -find ../PECA_scr/all_motif_rmdup -preparsedDir ../Homer/ > MotifTarget.bed
cat MotifTarget.bed|awk 'NR>1'|cut -f 1,4,6 > MotifTarget.txt
rm MotifTarget.bed
rm motifFindingParameters.txt

echo step 2: calculate opn...
cat openness.bed |tr '_' '\t' > openness1.bed
bedtools intersect -a openness1.bed -b ../PECA_scr/Opn_median_${genome}.bed -wa -wb -sorted|cut -f 1-4,8|sed 's/\t/_/1'|sed 's/\t/_/1'|sed 's/\t/_/1'|awk 'BEGIN{OFS="\t"}{ if ($2>a[$1] ) a[$1]=$2 }END{for (i in a) print i,a[i]}'|sed 's/_/\t/3' > openness2.bed
mkdir Enrichment
cat openness2.bed|awk 'BEGIN{OFS="\t"}{print $1,($2+0.5)/($3+0.5)}'|sort -k2nr|cut -f 1|tr '_' '\t'|awk 'BEGIN{OFS="\t"}{if ($3-$2 < 2000) print $0}'|head -10000 > ./Enrichment/region.bed
sed "s/species/${speciesFull}/g" ../PECA_scr/mf_collect.m > ./Enrichment/mf_collect.m 
cd ./Enrichment/
findMotifsGenome.pl region.bed ${genome} ./. -p 16 -size given -mask -nomotif -mknown ../../PECA_scr/all_motif_rmdup -preparsedDir ../../Homer/
module load matlab
matlab -nodisplay -nosplash -nodesktop -r "mf_collect; exit"
cd ../


echo step 3: Prior....
bedtools intersect -a region.bed -b ../PECA_scr/RE_gene_corr_${genome}.bed -wa -wb -sorted|cut -f 1-3,7-9|sed 's/\t/\_/1'|sed 's/\t/\_/1'>peak_gene_100k_corr.bed
bedtools intersect -a region.bed -b ../PECA_scr/Enhancer_RE_gene_corr_${genome}.bed -wa -wb -sorted|cut -f 1-3,7-9|sed 's/\t/\_/1'|sed 's/\t/\_/1'>>peak_gene_100k_corr.bed

echo step 4: Network....
cp ../PECA_scr/mfbs.m ./.
sed "s/toreplace/${input}/g" ../PECA_scr/PECA_network_${genome}.m > PECA_network.m
module load matlab
matlab -nodisplay -nosplash -nodesktop -r "PECA_network; exit"
echo region > CRbinding_region
cat region.txt|cut -f 4 > CRbinding_region1
cat CRbinding_region CRbinding_region1> CRbinding_region2
paste -d '\t' CRbinding_region2 CR_binding_pval.txt > CRB_pval.txt
rm CRbinding_region*
rm CR_binding_pval.txt
cat TGName.txt |tr '\n' '\t'> TG_head
echo >> TG_head
cat TG_head TFTG_regulationScore.txt > TG_score
echo TFName > TG_1
cat TG_1 TFName.txt >TG_2
paste -d '\t' TG_2 TG_score > TFTG_score.txt
rm TG_*
echo ${input} PS_PECA done
