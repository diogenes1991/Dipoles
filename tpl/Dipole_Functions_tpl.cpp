#include "####SubProcName####.h"

####SubProcName####::####SubProcName####(OLP& Prov, Model& Mod){

    Provider = &Prov;
    model = &Mod;
    nPar = NextR;
    RMomenta = new FVector[NextR];
    RParticles = new Particle*[NextR];

    BMomenta = new FVector[NextR-1];

####SubProcConst####

}

void ####SubProcName####::Subtracted(std::string cp, double sqrts, double* rand, double mu, double* rval){
/*
    double J = 1.0;
    double EWKFac = 4*M_PI*(model->alpha_e);
    double QCDFac = 4*M_PI*(model->alpha_s);
    double MatrixElement;
    FMatrix SpinCorr;
    CMatrix ColorCorr(NextR);
    double DipFac = 1.0;
    int BornIndex;
    *rval = 0.;

    OLP::Arguments Args;
    Args.mu_ren = mu;
    Args.RVal = &MatrixElement;
    Args.NExt = NextR;
    Args.Order = "LO";    
    Args.Momenta = RMomenta;
    Args.SCRVal = &SpinCorr;
    Args.CCRVal = &ColorCorr;

    SGenerate(sqrts,rand,&J);
    
####SubProcSub####

   else{
        std::cout << "Error: Coupling power: "<<cp<<" not found in process"<<std::endl;
        throw "Unrecognized Coupling";
   }

   *rval *= J;
*/
}

void ####SubProcName####::PlusDistribution(std::string cp, double sqrts, double* rand, double mu, double* rval){
/*
    int i;
    double Ix,I1;
    double EWKFac = model->alpha_e/(2*M_PI);
    double QCDFac = model->alpha_s/(2*M_PI);
    double DipFac = 1.0;
    *rval = 0.;

####SubProcPlu####

    else{
        std::cout << "Error: Coupling power: "<<cp<<" not found in process"<<std::endl;
        throw "Unrecognized Coupling";
   }
*/
}

void ####SubProcName####::Endpoint(std::string cp, double sqrts, double* rand, double mu, double* rval){
/*
    int i,j;
    double EWKFac = model->alpha_e/(2*M_PI);
    double QCDFac = model->alpha_s/(2*M_PI);
    double DipFac = 1.0;
    double J = 1.0;
    double acc,Invariant;
    std::vector<FVector> p;
    double pp[5*NextR];
    rval[0] = 0.;
    rval[1] = 0.;
    rval[2] = 0.;
    
####SubProcEnd####

    else{
        std::cout << "Error: Coupling power: "<<cp<<" not found in process"<<std::endl;
        throw "Unrecognized Coupling";
   }
*/
}