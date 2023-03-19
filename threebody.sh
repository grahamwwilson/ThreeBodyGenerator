#!/bin/sh

# Compile the code and link with ROOT 
./clroot.sh

# Read each line of TChiWZ.dat and calculate corresponding normalization constants
# using ThreeBodyIntegral executable
# Output file for reweighting is TChiWZgridWeights.dat

while read -r a b; do
    echo "$a" "$b"
    ./ThreeBodyIntegral -a ${a} -b ${b} -n 1000000
done < TChiWZ.dat 

exit
