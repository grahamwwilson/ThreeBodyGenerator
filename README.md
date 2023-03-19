# ThreeBodyGenerator

threebody.sh compiles and runs ThreeBodyGenerator.cpp to create a 
file like TChiWZgridWeights-V1-S13579.dat that contains the normalization 
constants such that the calculated pdf (see ReWeight.cpp) is correctly 
normalized.

ReWeight.cpp is an indicative implementation of reweighting. 
The user still needs to implement reading of the reweighting file.
