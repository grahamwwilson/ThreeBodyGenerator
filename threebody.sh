#!/bin/sh

./clroot.sh

while read -r a b; do
    echo "$a" "$b"
    ./ThreeBodyIntegral -a ${a} -b ${b} -n 1000000
done < TChiWZ.dat 

exit
