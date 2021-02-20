#include "Analysis.h"
#include <vector>

void Analysis::ReweightEvent(FVector* ExtMom, double beta, double* ExtMass, int* ExtPID, int BornNext, double* weight){
    
    // This is a sample Analysis file

    double wgt = 1;

    ///////////////////////////////////////////////////
    //
    // Gets all A and g particles from a channel
    // it cuts on the energy of them below some amount
    //
    ///////////////////////////////////////////////////

    // std::vector<FVector> radiation_mom;
    // for(int i=2;i<BornNext;i++){
    //     if(std::abs(ExtPID[i])==21||std::abs(ExtPID[i])==0)radiation_mom.push_back(ExtMom[i]);
    // }

    // for(FVector p : radiation_mom){
    //     if( p.p0 < 1) wgt = 0;
    // }

    //////////////////////////////////////////////////
    //
    // Gets all final b-quarks from a channel and cuts 
    // on both its eta and transverse momentum
    //
    ///////////////////////////////////////////////////

    std::vector<FVector> final_bs_mom;
    std::vector<int> final_bs_ids;
    for(int i=2;i<BornNext;i++){
        if(std::abs(ExtPID[i])==5)final_bs_mom.push_back(ExtMom[i]);
    }

    for(FVector p : final_bs_mom){
        FMatrix LabFrameBoost = Boost(0.0,0.0,beta);
        double eta = std::abs(Kinematics::PseudoRapidity(LabFrameBoost*p));
        double pt = Kinematics::TransverseMomentum(LabFrameBoost*p);
        if( pt < 25 || eta > 2.5 ) wgt = 0;
    }

    // FVector p = final_bs_mom.at(0) + final_bs_mom.at(1);
    // if (std::abs(sqrt(p*p)-91.1876) < 10) wgt = 0;
    
    *weight = wgt;

}

void Analysis::InitializeHistograms(std::vector<Histogram>* Histograms){
    
    //
    //  First Histogram 
    //

    //
    //  At position 1 we are storing the differential X-Section
    //  dsigma/dmbb
    //

    std::string VarName = "2b-jet Invariant Mass (GeV)";
    std::vector<double> Boundaries;
    for(int i=0;i<13000;i+=100) Boundaries.push_back(i);
    Histogram H1(VarName,Boundaries);
    Histograms->push_back(H1);

    // std::string Name2 = "E_Gamma_energy";
    // Histogram H2(Name2,Boundaries);
    // Histograms->push_back(H2);

}

void Analysis::FillHistograms(FVector* ExtMom, double beta, double* ExtMass, int* ExtPID, int BornNext, double xsec, std::vector<Histogram>* Histo){
    
    //
    //  Here we fill Histograms 
    //

    std::vector<FVector> final_bs_mom;
    for(int i=2;i<BornNext;i++){
        if(std::abs(ExtPID[i])==5)final_bs_mom.push_back(ExtMom[i]);
    }

    FVector p = final_bs_mom.at(0) + final_bs_mom.at(1);
    Histo->at(0).Append(sqrt(p*p),xsec);

    // std::vector<FVector> photon_mom;
    // for(int i=0;i<BornNext;i++){
    //     if(std::abs(ExtPID[i])==22) photon_mom.push_back(ExtMom[i]);
    // }

    // Histo->at(1).Append(photon_mom.at(0).p0,xsec);

}
