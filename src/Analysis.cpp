#include "Analysis.h"
#include <vector>

void Analysis::ReweightEvent(FVector* ExtMom, double* ExtMass, int* ExtPID, int BornNext, double* weight){
    
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

    // std::vector<FVector> final_bs_mom;
    // for(int i=2;i<BornNext;i++){
    //     if(std::abs(ExtPID[i])==5)final_bs_mom.push_back(ExtMom[i]);
    // }

    // for(FVector p : final_bs_mom){
    //     double eta = std::abs(Kinematics::PseudoRapidity(p));
    //     double pt = Kinematics::TransverseMomentum(p);
    //     if( pt < 25 || eta > 2.5 ) wgt = 0;
    // }
    
    *weight = wgt;

}