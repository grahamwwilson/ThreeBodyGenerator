
double threeBodyRelativeWeight(double x, double z, double mA, double mB, int signValue, double constant){

// For the process N_A -> N_B f fbar
// Arguments: Dalitz plot parameters (x, z). See Nojiri and Yamada paper (PRD 60 (1999) 015006))
//            Neutralino masses (mA, mB)
//            Mass matrix relative sign (+1 or -1)
//            constant: normalization constant. 
//            So that the mean weight of the modeled events coming from the allowed phase-space region is 1.0.

// Three choices for signValue = -1, 0, 1.
//  0  gives phase-space
//  1  gives "same-sign"
// -1  gives "opposite-sign"
 
    double mZ = 91.1876;
    double pdf;
    double rB = double(signValue)*mB/mA;
    double rBsq = rB*rB;
    double rZ = mZ/mA;
    double rZsq = rZ*rZ;
    double y = 1.0 + rBsq - x - z;
    
    if (signValue == -1 || signValue == 1){
    // Opposite sign (-1) or Same-Sign (+1)
        pdf = (1.0/constant) * ( (1.0-x)*(x-rBsq) + (1.0-y)*(y-rBsq) + 2.0*rB*z ) / std::pow(z - rZsq, 2.0);
    }
    else if (signValue == 0){
    // Phase-space only. Set to a constant. 
        pdf = 1.0;
    }
    else{
        std::cout << "shouldn't be here " << std::endl;
        pdf = 0.0;
    }
    return pdf;
}

int main(){

// Should read these values from the file (TChiWZgridWeights-V1-S13579.dat) 
// hard-coded here for illustration

    double mA = 300.0;
    double mB = 270.0;    
    double constPlus  = 2.82536;
    double constMinus = 0.90392;

// Not needed but for completeness
//    double fDalitz = 0.666147;
//    double percentPlus = 0.0383108;      // This is the fractional error in % on constPlus etc
//    double percentMinus = 0.0622303;
//    double percentDalitz = 0.05778;  
    
// Here we assume that x and z are valid values within the 3-body phase-space region (may need checks similar to those in ThreeBodyIntegral.cpp)   
    double OSwt = threeBodyRelativeWeight(x, z, mA, mB, -1, constMinus);
    double SSwt = threeBodyRelativeWeight(x, z, mA, mB,  1, constPlus);
    double PSwt = 1.0;
    
    double weightDesired = Newwt/Oldwt   // (eg. SSwt/PSwt).
 
// We need to take care that we don't end up over-weighting or under-weighting the events. 
// If there are N events generated each with weight 1, one needs to come up with a factor F such 
// that the sum of the new weights adds to N.  I had done that for the toy MC tests.
// This new implementation guarantees that the mean weight over all phase-space is 1.0.

// This new implementation contains the normalization constants.

// Without the normalization constant the pdf integral would be proportional to 
// the total width which will vary for the OS and SS models. 
// Since we only care that the particle has decayed somehow according to the 3 models, 
// rather than about what the modeled width is, we need the normalization constant to assure 
// that the number of decays remains constant.
   
}
