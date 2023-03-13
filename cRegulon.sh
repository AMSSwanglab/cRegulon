#/bin/bash

# cRegulon software v1.0.0 updated Mar 11th 2023 by Zhanying Feng

mkdir ./Results/${1}

echo "Begining cRegulon optimization..."
sed s/Sample/${i}/g ./scr/cRegulon.py > cRegulon.py
python cRegulon.py;rm cRegulon.py
echo "cRegulon optimization DONE!"

echo "Begining extracting TF module..."
sed s/Sample/${i}/g ./scr/TFModule.py > TFModule.py
python TFModule.py;rm TFModule.py
echo "TF module DONE!"

echo "Begining extracting subnetwork..."
sed s/Sample/${i}/g ./scr/SubNetwork.py > SubNetwork.py
python SubNetwork.py;rm SubNetwork.py
echo "Subnetwork DONE!"
