#!/bin/bash

echo "Preprocessing scRNA-seq and scATAC-seq data..."
Sample=$1
cat ./data/${Sample}_scATAC.txt | awk '{print $1}' | sed '1d' | grep -v GL | grep -v JH | awk -F'_' -v OFS='\t' '{print $1,$2,$3,$0}' > ${Sample}_Peaks.bed
bedtools intersect -wa -wb -a ./scr/Mm10_GoodGene.bed -b ${Sample}_Peaks.bed | awk -v OFS='\t' '{print $4"--"$5,$9}' > ${Sample}_Gene_NearPeak.txt
./scr/MergeLine ${Sample}_Gene_NearPeak.txt
sed s/Sample/${Sample}/g ./scr/NearPeakDis.py > NearPeakDis.py
python NearPeakDis.py;rm NearPeakDis.py
sed s/Sample/${Sample}/g ./scr/GA.py > GA.py
python GA.py;rm GA.py
sed s/Sample/${Sample}/g ./scr/GE.py > GE.py
python GE.py;rm GE.py
echo "Preprocessing end: ${Sample}_GE.txt and ${Sample}_GA.txt generated!"
