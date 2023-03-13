#/bin/bash

# TF-TF CSI combinatorial network v1.0.0 updated Mar 11th 2023 by Zhanying Feng
sed s/Sample/${1}/g ./scr/CSI.py > CSI.py
python CSI.py;rm CSI.py
rm -f Sample_GA.txt Sample_GE.txt
