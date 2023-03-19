# ThreeBodyGenerator

ReWeight.cpp is an indicative implementation of reweighting. 
The user still needs to implement reading of the supplied 
reweighting file. The supplied TChiWZgridWeights-V1-S13579.dat uses 
one million in-bounds (in allowed Dalitz region) MC events per grid point.

If needed for other grids or higher statistics, 
threebody.sh compiles and runs ThreeBodyGenerator.cpp and uses an input csv 
like file to create a file like TChiWZgridWeights-V1-S13579.dat that contains 
the normalization constants such that the calculated pdf (see ReWeight.cpp) is correctly 
normalized.

The sampling is done in (x,z) which is much more efficient than using (x,y). 
Typically 66% of candidate phase space points are valid (in-bounds).

Graham W. Wilson
