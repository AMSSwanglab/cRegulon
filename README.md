# cRegulon

## Introduction
This is cRegulon software: an optimization model to identify combinatorial regulon from single cell expression and chromatin accessibility data.
## Requirements:
1. Python >=3.0 with packages: numpy, sklearn, and scipy <br>
2. matlab >= 2021
3. Homer

### Installing cRegulon with the following command:
```bash
wget https://github.com/fengzhanying/cRegulon/archive/master.zip
unzip master.zip
cd cRegulon-master
wget https://www.dropbox.com/s/0h1wxlu7iqheajo/cRegulon.tar.gz
tar -xzvf cRegulon.tar.gz
```
## Training mode of cRegulon
### Input single cell data
The typic input file (CL_scRNA.txt) of scRNA-seq data is a gene by cell count matrix: <br>
<table>
  <tr>
    <td>scRNA</td>
    <td>RNACell1</td>
    <td>RNACell2</td>
    <td>RNACell3</td>
  </tr>
  <tr>
    <td>Gene1</td>
    <td>5</td>
    <td>0</td>
    <td>3</td>
  </tr>
  <tr>
    <td>Gene2</td>
    <td>0</td>
    <td>2</td>
    <td>0</td>
  </tr>
  <tr>
    <td>Gene3</td>
    <td>1</td>
    <td>0</td>
    <td>0</td>
  </tr>
</table>
The typic input file (CL_scATAC.txt) of scATAC-seq data is a peak by cell count matrix:
<table>
  <tr>
    <td>scATAC</td>
    <td>ATACell1</td>
    <td>ATACell2</td>
    <td>ATACell3</td>
    <td>ATACell4</td>
  </tr>
  <tr>
    <td>Peak1</td>
    <td>1</td>
    <td>0</td>
    <td>1</td>
    <td>0</td>
  </tr>
  <tr>
    <td>Peak2</td>
    <td>0</td>
    <td>1</td>
    <td>0</td>
    <td>1</td>
  </tr>
  <tr>
    <td>Peak3</td>
    <td>1</td>
    <td>0</td>
    <td>0</td>
    <td>0</td>
  </tr>
</table>
The peaks are in the format of "chr_start_end". <br>

### Input cell type meta data
The typic cell type meta file (CL_scRNA_Cluster.txt) of scRNA-seq data is as follows: <br>
<table>
  <tr>
    <td>RNACell1</td>
    <td>C1</td>
  </tr>
  <tr>
    <td>RNACell2</td>
    <td>C2</td>
  </tr>
  <tr>
    <td>RNACell3</td>
    <td>C3</td>
  </tr>
</table>
The typic cell type meta file (CL_scATAC_Cluster.txt) of scATAC-seq data is as follows: <br>
<table>
  <tr>
    <td>ATACell1</td>
    <td>C2</td>
  </tr>
  <tr>
    <td>ATACell2</td>
    <td>C3</td>
  </tr>
  <tr>
    <td>ATACell3</td>
    <td>C1</td>
  </tr>
  <tr>
    <td>ATACell4</td>
    <td>C2</td>
  </tr>
</table>

### Step 1: GRN construction
We run the following script to make the gene expression matrix and gene activity matrix (current we support hg38 and mm10):

```bash
source runNet.sh CL hg38
```
This process will produce GRN files (network.txt, TFName.txt, TGName.txt, TRS.txt) for each cell cluster.

### Step 2: Running cRegulon model
With the input of TF-TF combinatorial network (RAd4_CSI.txt), normalized TF-TG regulatory strength matrix (RAd4_TRS.txt), gene expression matrix (RAd4_GE.txt), and gene activity matrix (RAd4_GA.txt), we run the following cRegulon model:
```bash
python cRegulon.py CL hg38
```
This will output: <br>
1. TF combinatorial effects in each cRegulon: X.txt <br>
2. Association matrix between cell clusters and cRegulons: L.txt <br>
3. TF module of each cRegulon: TFs (*TF.txt) and TF pairs (*TFPair.txt).
5. Regulatory sub-network of each cRegulon: *SubNet.txt

## Annotation mode of cRegulon

## Citation:
If you use cRegulon software or cRegulon associated concepts, please cite:

Zhanying Feng, et al. Modeling combinatorial regulation from single-cell multi-omics provides units underpinning cell type landscape. 2024.
