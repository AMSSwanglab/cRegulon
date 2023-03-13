# cRegulon
cRegulon is an optimization model to identify combinatorial regulon from single cell expression and chromatin accessibility data.

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
## Step 1: single cell data preprocessing
The typic input file (RAd4_scRNA.txt) of scRNA-seq data is a gene by cell count matrix: <br>
<table>
  <tr>
    <td>scRNA</td>
    <td>Cell1</td>
    <td>Cell2</td>
    <td>Cell3</td>
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
The typic input file (RAd4_scATAC.txt) of scATAC-seq data is a peak by cell count matrix:
<table>
  <tr>
    <td>scATAC</td>
    <td>Cell1</td>
    <td>Cell2</td>
    <td>Cell3</td>
    <td>Cell4</td>
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
We run the following script to make the gene expression matrix and gene activity matrix:

```bash
source Preprocessing.sh RAd4
```
This process will produce gene expression file (RAd4_GE.txt) and gene activity file (RAd4_GA.txt)

## Step 2: Constructing TF-TG regulatory network by pseudo-bulk strategy
With the input files are (RAd4_scRNA.txt) and (RAd4_scATAC.txt), we run the following script:
```bash
source PS_PECA.sh RAd4 mm10
```
This process will produce the TF-REs-TG triplets files (RAd4_network.txt) and TF-TG regulatory strength file (RAd4_TRS.txt).
## Step 3: Constructing TF-TF combinatorial network
With the input TF-TG regulatory strength file (RAd4_TRS.txt), we run the following script:
```bash
source runCSI.sh RAd4
```
This will generate normalized TF-TG regulatory strength file (RAd4_TRS.txt) and TF-TF combinatorial network (RAd4_CSI.txt).
## Step 4: Running cRegulon model
With the input of TF-TF combinatorial network (C.txt), normalized TF-TG regulatory strength matrix (R.txt), gene expression matrix (GE.txt), and gene activity matrix (GA.txt), we run the following cRegulon model:
```bash
source cRegulon.sh RAd4
```
This will output: <br>
1. TF combinatorial effects in all cRegulons: X.txt <br>
2. cRegulon combination coefficients for scRNA-seq: H1.txt <br>
3. cRegulon combination coefficients for scATAC-seq: H2.txt <br>
4. TF modules of cRegulons: TFs (*TF.txt) and TF pairs (*TFPair.txt).
5. Regulatory sub-network of each cRegulon: *SubNet.txt


## Citation:
If you use cRegulon software or cRegulon associated concepts, please cite

Zhanying Feng, et al. Modeling combinatorial regulon from single cell gene expression and chromatin accessibility data. 2023.
