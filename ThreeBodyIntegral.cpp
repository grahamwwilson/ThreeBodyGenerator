#include "CLI11.hpp"
//#include <CLI11.hpp>   // CLI11 command line interface stuff (V2.3.1). See https://github.com/CLIUtils/CLI11
#include <cmath>
#include "TH2.h"
#include "TFile.h"
#include <iostream>
#include <fstream>
#include <random>
#include <set>
#include <map>
#include <string>

typedef std::mt19937 RandomNumberGenerator;

double kallen(double x, double y, double z){
// Kallen kinematical function. See Byckling and Kajantie II.6.3. 
    return std::pow(x - y - z, 2.0) - 4.0*y*z;
}

double angleFormula(double s, double mi, double mj, double sj, double sk, double si){
   //
   // Compute the opening angles in the 3-body rest frame using cyclic permutations of BK equation V.1.4.
   //                double costh12 = angleFormula(s, m1, m2, s23, s31, s12); 

       double costh;
       
       costh = ( (s + mi*mi - sj) * (s + mj*mj - sk) + 2.0*s*(mi*mi + mj*mj - si) ) /
               std::sqrt(kallen(s, mi*mi, sj) * kallen(s, mj*mj, sk) );

       return costh;
}

double decayAngleFormula(double s, double mi, double mj, double mk, double sj, double si){
   //
   // Compute the decay angle in the 2-particle rest frame using cyclic permutations of BK equation V.1.9.
   // Using the BK notation, Rij is the rest-frame for particles i and j, with pi + pj = 0.
   //
   // double costh23star = decayAngleFormula(s, m2, m3, m1, s23, s12);

       double costh;
       
       costh = ( (s - sj - mk*mk) * (sj + mi*mi - mj*mj) + 2.0*sj*(mk*mk + mi*mi - si) ) /
               std::sqrt(kallen(s, sj, mk*mk) * kallen(sj, mi*mi, mj*mj) );

       return costh;
}

std::pair<double, double> Calculations(std::string which, long int denom, double normConstant, double wtsum, double wtsumsq){

    std::cout << " " << std::endl;
    std::cout << "Calculations for " << which << " with " << denom << " in-bounds events " << std::endl; 
 //   std::cout << " " << std::endl;    
    double fbar = wtsum/double(denom);
    double ffbar = wtsumsq/double(denom);
    double sigf = std::sqrt(ffbar - fbar*fbar);
 //   std::cout << "MC integral " << volume*fbar << " +- " << volume*sigf/std::sqrt(double(denom)) << std::endl;
    std::cout << "MC integral normalization check " << fbar << " +- " << sigf/std::sqrt(double(denom)) << std::endl;
    std::cout << "Deviation from unity: " << fbar - 1.0 << std::endl;
    
    double c = normConstant*fbar;
    double dc = normConstant*sigf/std::sqrt(double(denom));    
    
    std::cout << "Current normconstant: " << normConstant << " Suggested revised value of normConstant : " 
              << c << " +- " << dc << std::endl;
    double dp = 100.0*dc/c;
    std::cout << "Percentage uncertainty " << dp << "% " << std::endl;
              
    std::pair<double, double> p = std::make_pair(c, dp);
    return p;      
}

void BiGenerator(int nevents, unsigned long int seed, double mA, double mB, 
                 double mf, int signValue, double WTMAX, double pdfave, 
                 bool passthrough, double normConstantPlus, double normConstantMinus, double normConstantFlat){

// Sample from the 3-body phase space for A -> l1 l2 B based on the presentation 
// in Nojiri and Yamada. This currently neglects fermion masses.
// Although at least it should be relatively trivial to incorporate 
// fermion masses on the allowed phase-space if not the 2-d matrix element-squared itself.
//
// x = [m(l1,N)/mA]**2
// y = [m(l2,N)/mA]**2
// z = [m(l1,l2)/mA]**2
// 
// We define r = +- (mB/mA) where both signs are feasible
// rZ = mZ/mA

    std::ofstream fout;   
    
    int imA = int(mA + 1.0e-8);
    int imB = int(mB + 1.0e-8);
    if(imA == 100 && imB == 50){
// Write header line to file (hard-coded that the first one is (100, 50)
        fout.open("TChiWZgridWeights.dat");
        fout << "m_N2 m_N1 constPlus percentPlus constMinus percentMinus fDalitz percentDalitz" << std::endl;    
    }
    else{
        fout.open("TChiWZgridWeights.dat", std::ios_base::app);  // append to existing file   
    }
    fout.precision(6);

    double mZ = 91.1876;
    double maxweight = -1.0;
    double rB = (mB/mA);
    double rBsq = rB*rB;
    double rZsq = mZ*mZ/(mA*mA);
    double zsqmax = std::pow( (mA-mB)/mA, 2.0);

    RandomNumberGenerator gx(seed);
    std::uniform_real_distribution<double> uniformx(rBsq, 1.0);
    RandomNumberGenerator gz(seed+1);
    std::uniform_real_distribution<double> uniformz(0.0, zsqmax); 
    RandomNumberGenerator gu(seed+2);
    std::uniform_real_distribution<double> uniform;    
        
    double u1,u2;
    
    long int ntrials = 0;
    long int ngenerated = 0;
    long int noutofbounds = 0;
    long int ninbounds = 0;
    
    double testQuantity;
    
    double Eneutmax = (mA*mA + mB*mB)/(2.0*mA);
    
    std::unique_ptr<TFile> myFile( TFile::Open("PhaseSpace.root","RECREATE") );
    TH1D *h_mll = new TH1D("h_mll","Dilepton Mass (GeV)",120,0.0,mA-mB);
    TH1D *h_mllp = new TH1D("h_mllp","Dilepton Mass (GeV)",120,0.0,mA-mB);
    TH1D *h_mllu = new TH1D("h_mllu","Dilepton Mass (GeV)",120,0.0,mA-mB);    
    TH1D *h_pneut = new TH1D("h_pneut","Neutralino Momentum (GeV)",114,0.0,std::sqrt(Eneutmax*Eneutmax - mB*mB));    
    TH1D *h_mllR = new TH1D("h_mllR","Reweighted Dilepton Mass (GeV)",120,0.0,mA-mB);    
    TH1D *h_costh12 = new TH1D("h_costh12","costh12",100,-1.0,1.0);
    TH1D *h_costh23 = new TH1D("h_costh23","costh23",100,-1.0,1.0);
    TH1D *h_costh31 = new TH1D("h_costh31","costh31",100,-1.0,1.0); 
    TH1D *h_pt = new TH1D("h_pt","pT (lepton) (GeV)",100,0.0,0.5*(mA-mB));
    TH1D *h_costh12star = new TH1D("h_costh12star","costh12star",100,-1.0,1.0);
    TH1D *h_costh12starR = new TH1D("h_costh12starR","Reweighted costh12star",100,-1.0,1.0);    
    TH1D *h_costh23star = new TH1D("h_costh23star","costh23star",100,-1.0,1.0); 
    TH1D *h_costh31star = new TH1D("h_costh31star","costh31star",100,-1.0,1.0); 
    TH2D *h_xy = new TH2D("h_xy","y vs x",100,rBsq,1.0,100,rBsq,1.0); 
    TH2D *h_xz = new TH2D("h_xz","x vs z",50,0.0,zsqmax,50,rBsq,1.0);
    TH2D *h_xzR = new TH2D("h_xzR","Reweighted x vs z",100,0.0,zsqmax,100,rBsq,1.0);
    
    TH1D *h_weight = new TH1D("h_weight","weight",600,0.0,6.0);
    TH1D *h_lowweight = new TH1D("h_lowweight","weight",500,0.0,0.1);    
    TH1D *h_accweight = new TH1D("h_accweight","accepted weight",600,0.0,6.0); 
    TH1D *h_wtweighted = new TH1D("h_wtweighted","weighted weight",600,0.0,6.0);
    
    TH1D *h_weight1 = new TH1D("h_weight1","weight1",600,0.0,6.0);
    TH1D *h_lowweight1 = new TH1D("h_lowweight1","weight1",500,0.0,0.1);     
    TH1D *h_wtweighted1 = new TH1D("h_wtweighted1","weighted weight1",600,0.0,6.0);  
    
    TH1D *h_weight2 = new TH1D("h_weight2","weight2",600,0.0,6.0);
    TH1D *h_lowweight2 = new TH1D("h_lowweight2","weight2",500,0.0,0.1);     
    TH1D *h_wtweighted2 = new TH1D("h_wtweighted2","weighted weight2",600,0.0,6.0);  
    
    TH1D *h_weight3 = new TH1D("h_weight3","weight3",600,0.0,6.0);
    TH1D *h_lowweight3 = new TH1D("h_lowweight3","weight3",500,0.0,0.1);     
    TH1D *h_wtweighted3 = new TH1D("h_wtweighted3","weighted weight3",600,0.0,6.0);                                  
    
    double pdfsum = 0.0; 
    double pdfsumsq = 0.0;
    double wt1sum = 0.0;
    double wt1sumsq = 0.0;
    double wt2sum = 0.0;
    double wt2sumsq = 0.0;
    double wt3sum = 0.0;
    double wt3sumsq = 0.0;        
    
    double volume = (1.0-rBsq)*zsqmax;   // Volume of the MC integral
    std::cout << "Cell volume  = " << volume << std::endl;
    
    while (ngenerated < nevents){
    
// Try (x,z on [rB**2, 1.0] x [0.0, (dM/M)**2])
        double x = uniformx(gx);
        double z = uniformz(gz);
        double y = 1.0 + rBsq - x - z;
        testQuantity = z*(x*y - rBsq);
        
        ntrials++;
        
        if (ntrials <= 5){
            std::cout << "trial " << ntrials << " " << x << " " << y << " " << z << " " << testQuantity << std::endl;
        }        
        
        if (testQuantity >= 0.0) {
  // We're in the allowed region of the Dalitz plot. 
            ninbounds++;
            double pdfPlus =  (1.0/normConstantPlus) * ( (1.0-x)*(x-rBsq) + (1.0-y)*(y-rBsq) + 2.0*std::abs(rB)*z ) / std::pow(z - rZsq, 2.0);
            double pdfMinus = (1.0/normConstantMinus) * ( (1.0-x)*(x-rBsq) + (1.0-y)*(y-rBsq) - 2.0*std::abs(rB)*z ) / std::pow(z - rZsq, 2.0);
            double pdfFlat =  (1.0/normConstantFlat) * 1.0;
                        
            h_weight1->Fill(pdfPlus);
            h_lowweight1->Fill(pdfPlus);
            h_wtweighted1->Fill(pdfPlus, pdfPlus);

            h_weight2->Fill(pdfMinus);
            h_lowweight2->Fill(pdfMinus);
            h_wtweighted2->Fill(pdfMinus, pdfMinus);
            
            h_weight3->Fill(pdfFlat);
            h_lowweight3->Fill(pdfFlat);
            h_wtweighted3->Fill(pdfFlat, pdfFlat);
            
// Make pdf choice for weighting of standard histograms.
            double pdf;
            double pdfOther = 0.0;
            if(signValue == -1){
                pdf = pdfMinus;
                pdfOther = pdfPlus;
            }
            if(signValue ==  1){
                pdf = pdfPlus;
                pdfOther = pdfMinus;
            }
            if(signValue ==  0)pdf = pdfFlat;            
            
            wt1sum += pdfPlus;
            wt1sumsq += pdfPlus*pdfPlus;
            wt2sum += pdfMinus;
            wt2sumsq += pdfMinus*pdfMinus;
            wt3sum += pdfFlat;
            wt3sumsq += pdfFlat*pdfFlat;                        
            
            if ( passthrough ){
   // KEEP these random variates and weight according to chosen pdf above.        
                ngenerated++; 
                h_accweight->Fill(pdf);
                
                pdfsum += pdf;
                pdfsumsq += pdf*pdf;
                 
   // Transform back to masses
                double mll  = sqrt(z*mA*mA);
                double ml1N = sqrt(x*mA*mA);
                double ml2N = sqrt(y*mA*mA);
                if(ngenerated <= 5)std::cout << "Accept  " << mll << " " << ml1N << " " << ml2N << std::endl; 
                
   // Set up variables as in Byckling and Kajantie chapter 5. But let's dispense with the confusing s1,s2,s3
                double s   = mA*mA;
                double s12 = mll*mll;          // Was s1
                double s23 = ml2N*ml2N;        // Was s2
                double s31 = ml1N*ml1N;        // Was s3
                double m1  = mf;
                double m2  = mf;
                double m3  = mB;
                
   // Compute the opening angles in the 3-body rest frame using cyclic permutations of equation V.1.4.
                double costh12 = angleFormula(s, m1, m2, s23, s31, s12);
                double costh23 = angleFormula(s, m2, m3, s31, s12, s23);
                double costh31 = angleFormula(s, m3, m1, s12, s23, s31);
   // Compute the energies and momenta in the 3-body rest frame             
                double E1 = (s + m1*m1 - s23)/(2.0*std::sqrt(s));
                double E2 = (s + m2*m2 - s31)/(2.0*std::sqrt(s));
                double E3 = (s + m3*m3 - s12)/(2.0*std::sqrt(s));             
                double p1 = std::sqrt(E1*E1 - m1*m1);
                double p2 = std::sqrt(E2*E2 - m2*m2);
                double p3 = std::sqrt(E3*E3 - m3*m3);
   // Compute the "scattering angle" in the two-body rest frame using equation V.1.9             
                double costh23star = decayAngleFormula(s, m2, m3, m1, s23, s12);   //costh_12^R23
                double costh31star = decayAngleFormula(s, m3, m1, m2, s31, s23);   //costh_23^R31
                double costh12star = decayAngleFormula(s, m1, m2, m3, s12, s31);   //costh_31^R12
                
   // Lepton pt with respect to neutralino3 direction
                double pt = p1*sin(acos(costh31));                                 
                
                h_mll->Fill(mll, pdf);
                h_mllp->Fill(mll, pdfOther);
                h_mllu->Fill(mll, pdfFlat);                                
                h_costh12->Fill(costh12, pdf);
                h_costh23->Fill(costh23, pdf);
                h_costh31->Fill(costh31, pdf);
                h_pt->Fill(pt, pdf);
                h_costh12star->Fill(costh12star, pdf); 
                h_costh23star->Fill(costh23star, pdf); 
                h_costh31star->Fill(costh31star, pdf);
                h_xy->Fill(x,y,pdf);
                h_xz->Fill(z,x,pdf);
                h_pneut->Fill(p3,pdf); 
                
                // reweighting tests under either Odd or Even hypothesis
                /*
                double wt = pdf/pdfave;                
                h_mllR->Fill(mll,wt);
                h_costh12starR->Fill(costh12star,wt);
                h_xzR->Fill(z,x,wt);
                */
            }
        }    
        else{ 
            noutofbounds++;
        }
    }
    
    std::cout << "rBsq " << rBsq << " sign = " << signValue << std::endl;
//    std::cout << "WTMAX set to " << WTMAX << std::endl;
//    std::cout << "Max weight observed of " << maxweight << std::endl;
    std::cout << "Ntrials, NOB, NIB, Ngen, Efficiency, EfficiencyR " << ntrials << " " << noutofbounds << " " << ninbounds 
              << " " << ngenerated << " " << double(ngenerated)/double(ntrials) << " " << double(ngenerated)/double(ninbounds) << std::endl;
    std::cout << "Mean pdf value for generated events " << pdfsum/double(ngenerated) << std::endl;
 
    std::cout.precision(12); 
        
    // Phase-space fraction measurement.
    double pcontained = double(ninbounds)/double(ntrials);
    double dpcontained = std::sqrt(pcontained*(1.0-pcontained)/double(ntrials));
    std::cout << "Dalitz region fraction = " << pcontained << " +- " << dpcontained << std::endl;
    
    long int denominator = ninbounds;
    
    auto pPlus  = Calculations("Plus", denominator, normConstantPlus, wt1sum, wt1sumsq);
    auto pMinus = Calculations("Minus", denominator, normConstantMinus, wt2sum, wt2sumsq);
    auto pFlat  = Calculations("Flat", denominator, normConstantFlat, wt3sum, wt3sumsq);
    
    fout << imA << " " << imB << " " << pPlus.first << " " << pPlus.second << " " << pMinus.first << " " << pMinus.second << " " << pcontained << " " << 100.0*dpcontained/pcontained << std::endl;       

/*        
    double fbar = wt1sum/double(denominator);
    double ffbar = wt1sumsq/double(denominator);
    double sigf = std::sqrt(ffbar - fbar*fbar);
    std::cout << "MC integral " << volume*fbar << " +- " << volume*sigf/std::sqrt(double(denominator)) << std::endl;
    std::cout << "MC integral normalization check " << fbar << " +- " << sigf/std::sqrt(double(denominator)) << std::endl;
    std::cout << "Deviation from unity: " << fbar - 1.0 << std::endl;
    std::cout << "Current normconstant: " << normConstantPlus << " Suggested revised value of normConstantPlus : " 
              << normConstantPlus*fbar << " +- " << normConstantPlus*sigf/std::sqrt(double(denominator)) << std::endl; 
*/
    
    h_mll->Write();
    h_mllu->Write();
    h_mllp->Write();
    h_costh12->Write();
    h_costh23->Write();
    h_costh31->Write();
    h_pt->Write();
    h_costh12star->Write();
    h_costh23star->Write();
    h_costh31star->Write();
    h_xy->Write();
    h_xz->Write();
    h_mllR->Write();
    h_costh12starR->Write();
    h_xzR->Write();
    h_pneut->Write();
    h_weight->Write();
    h_accweight->Write();
    h_wtweighted->Write();
    h_lowweight->Write();
    myFile->Close();
    
    fout.close();
    
}

int main(int argc, char **argv) {

    CLI::App app{"Sample from 3-body decay phase-space with appropriate matrix element"};
    
    int nevents = 1000000;
    app.add_option("-n,--nevents", nevents, "Number of events");    

    unsigned long int seed = 13579L;
    app.add_option("-s,--seed", seed, "Seed");
    
    double mA = 300.0;
    app.add_option("-a,--ma", mA, "Parent neutralino mass (GeV)");
    
    double mB = 270.0;
    app.add_option("-b,--mb", mB, "Child neutralino mass (GeV)");    
    
    double mf = 0.0;
    app.add_option("-f,--mf", mf, "Fermion mass (GeV)");
    
    int signValue = -1;
    app.add_option("-e,--evsign", signValue, "Relative mass eigenvalue sign (+-1)");
    
    double WTMAX = 6.0;
    app.add_option("-w,--weightmax", WTMAX, "Maximum weight");     

//    std::string filename = "CopulaGen.EDAT";
//    app.add_option("-o,--outputfile", filename, "Output copula file");
    
    double pdfave = 1.25391;     //for signValue = -1
//    double pdfave=3.23978;     // for signValue = +1    
    app.add_option("-r,--reweight",pdfave, "Average pdf value for reweighting");
    
    double normConstantPlus = 2.8245;    
    app.add_option("--constPlus",normConstantPlus, "Normalization constant (+) for pdf");  
    
    double normConstantMinus = 0.9039;    
    app.add_option("--constMinus",normConstantMinus, "Normalization constant (-) for pdf"); 
    
    double normConstantFlat = 1.0;    
    app.add_option("--constFlat",normConstantFlat, "Normalization constant (0) for pdf");            
    
    bool passthrough = true;
    app.add_flag("-p,--passthrough", passthrough, "Disable event rejection based on matrix element weight");    
    
    CLI11_PARSE(app, argc, argv);

    std::cout.precision(12);
    std::cout << "nevents      " << nevents << std::endl;
    std::cout << "seed         " << seed << std::endl;
    std::cout << "mA           " << mA << std::endl;
    std::cout << "mB           " << mB << std::endl;
    std::cout << "mf           " << mf << std::endl;
//    double signValue = double(signvalue);
    std::cout << "signValue    " << signValue << std::endl;   
    std::cout << "WTMAX        " << WTMAX << std::endl; 
    std::cout << "pdfave       " << pdfave << std::endl;
    std::cout << "normConstantPlus " << normConstantPlus << std::endl;
    std::cout << "normConstantMinus " << normConstantMinus << std::endl; 
    std::cout << "normConstantFlat " << normConstantFlat << std::endl;       
//    std::cout << "filename     " << filename << std::endl;
    std::cout << "passthrough  " << passthrough << std::endl;    
    
    BiGenerator(nevents, seed, mA, mB, mf, signValue, WTMAX, pdfave, passthrough, normConstantPlus, normConstantMinus, normConstantFlat);    
       
    return 0;
    
}
