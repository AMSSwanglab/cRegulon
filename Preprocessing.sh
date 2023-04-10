#!/bin/bash

echo "Preprocessing scATAC-seq data..."
Sample=$1;Genome=$2
cat ./data/${Sample}_scATAC.txt | awk '{print $1}' | sed '1d' | grep -v GL | grep -v JH | awk -F'_' -v OFS='\t' '{print $1,$2,$3,$0}' > ${Sample}_Peaks.bed
bedtools intersect -wa -wb -a ./scr/GoodGene_${Genome}.bed -b ${Sample}_Peaks.bed | awk -v OFS='\t' '{print $4"--"$5,$9}' > ${Sample}_Gene_NearPeak.txt
./scr/MergeLine ${Sample}_Gene_NearPeak.txt
sed s/Sample/${Sample}/g ./scr/NearPeakDis.py > NearPeakDis.py
python3 NearPeakDis.py;rm NearPeakDis.py
sed s/Sample/${Sample}/g ./scr/GA.py | sed s/Species/${Genome}/g > GA.py
python3 GA.py;rm GA.py ${1}_Gene_NearPeakDis.txt ${1}_Gene_NearPeak.txt ${1}_Peaks.bed

echo "Preprocessing scRNA-seq data..."
sed s/Sample/${Sample}/g ./scr/GE.py | sed s/Species/${Genome}/g > GE.py
python3 GE.py;rm GE.py
echo "Preprocessing end: ${Sample}_GE.txt and ${Sample}_GA.txt generated!"
