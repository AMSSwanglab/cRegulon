#/bin/bash

Name=$1
genome=$2

if [ `echo $genome|grep hg|wc -l` -gt 0 ]
then
	species=hs
	speciesFull=human
else
	species=mm
	speciesFull=mouse
fi

mkdir $Name;
cd $Name
cat ../PseudoBulk/${Name}_PS_RNA.txt > ${Name}.txt
cat ../PseudoBulk/${Name}_PS_ATAC.txt | grep -v GL | grep -v JH | awk '$2>2' | tr '_' '\t' | sortBed | awk -v OFS='\t' '{print $1"_"$2"_"$3,$4}' > openness.bed
cat ../PseudoBulk/${Name}_PS_ATAC.txt | grep -v GL | grep -v JH | awk '$2>2' | tr '_' '\t' | sortBed | awk -v OFS='\t' '{print $1,$2,$3}' > region.bed
cat ../PseudoBulk/${Name}_PS_ATAC.txt | grep -v GL | grep -v JH | awk '$2>2' | tr '_' '\t' | sortBed | awk -v OFS='\t' '{print $1,$2,$3,$1"_"$2"_"$3}' > region.txt

cat openness.bed |tr '_' '\t' > openness1.bed
bedtools intersect -a openness1.bed -b ../Data/Opn_median_${genome}.bed -wa -wb -sorted|cut -f 1-4,8|sed 's/\t/_/1'|sed 's/\t/_/1'|sed 's/\t/_/1'|awk 'BEGIN{OFS="\t"}{ if ($2>a[$1] ) a[$1]=$2 }END{for (i in a) print i,a[i]}'|sed 's/_/\t/3' > openness2.bed
mkdir Enrichment
cat openness2.bed|awk 'BEGIN{OFS="\t"}{print $1,($2+0.5)/($3+0.5)}'|sort -k2nr|cut -f 1|tr '_' '\t'|awk 'BEGIN{OFS="\t"}{if ($3-$2 < 2000) print $0}'|head -10000 > ./Enrichment/region.bed
sed "s/species/${speciesFull}/g" ../scr/mf_collect.m > ./Enrichment/mf_collect.m 
cd ./Enrichment/
findMotifsGenome.pl region.bed ${genome} ./. -p 16 -size given -mask -nomotif -mknown ../../Data/all_motif_rmdup -preparsedDir ../../Data/Homer/
module load matlab
matlab -nodisplay -nosplash -nodesktop -r "mf_collect; exit"
cd ../

bedtools intersect -a region.bed -b ../Data/RE_gene_corr_${genome}.bed -wa -wb -sorted|cut -f 1-3,7-9|sed 's/\t/\_/1'|sed 's/\t/\_/1'>peak_gene_100k_corr.bed
bedtools intersect -a region.bed -b ../Data/Enhancer_RE_gene_corr_${genome}.bed -wa -wb -sorted|cut -f 1-3,7-9|sed 's/\t/\_/1'|sed 's/\t/\_/1'>>peak_gene_100k_corr.bed

cp ../scr/mfbs.m ./.
sed "s/toreplace/${Name}/g" ../scr/PECA_network_${speciesFull}.m > PECA_network.m
module load matlab
matlab -nodisplay -nosplash -nodesktop -r "PECA_network; exit"
cd ../
