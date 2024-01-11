#/bin/bash

input=$1
genome=$2

echo step 1: Pseudo-bulk generating....
python3 ./scr/PseudoBulk.py -n CL -r ./Input/CL_scRNA.txt -rm ./Input/CL_scRNA_Cluster.txt -a ./Input/CL_scATAC.txt -am ./Input/CL_scATAC_Cluster.txt

cd PseudoBulk
findMotifsGenome.pl region.txt hg38 ./. -p 16 -size given -find ../Data/all_motif_rmdup -preparsedDir ../Data/Homer/ > MotifTarget.bed
cat MotifTarget.bed|awk 'NR>1'|cut -f 1,4,6 > MotifTarget.txt
rm MotifTarget.bed
rm motifFindingParameters.txt
cd ../

echo step 2: Generating GRN...
for ct in ls `cat ./PseudoBulk/CL_CellType.txt`
do
	source PECA.sh $ct $genome
	mv ./${ct}/${ct}_network.txt ./CL_Networks/
	mv ./${ct}/TFName.txt ./Networks/${ct}_TFName.txt
	mv ./${ct}/TGName.txt ./Networks/${ct}_TGName.txt
	mv ./${ct}/TFTG_regulationScore.txt ./Networks/${ct}_TRS.txt
	rm -rf $ct
done
