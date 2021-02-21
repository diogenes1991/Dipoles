#include "####SubProcName####.h"

####SubProcName####::####SubProcName####(OLP& Prov, Model& Mod){

    Provider = &Prov; 
    model = &Mod;
    nPar = NextV; 
    Particles = new Particle*[NextV]; 
    Momenta   = new FVector[NextV];
    
####SubProcConst####

}

void ####SubProcName####::Born(std::string cp, double sqrts, double* rand, double mu, double* rval){

    double J = 1.0;
    *rval = 0.;

    BGenerate(sqrts,rand,&J);
    
    OLP::Arguments Args;
    Args.mu_ren = mu;
    Args.RVal = rval;
    Args.NExt = NextV;
    Args.Order = "LO";    
    Args.Momenta = Momenta; 
           
    
####SubProcBorn####

    else{
        std::cout << "Error: Coupling power: "<<cp<<" not found in process"<<std::endl;
        throw "Unrecognized Coupling";
    }
    
    *rval *= J;
}

void ####SubProcName####::Virtual(std::string cp, double sqrts, double* rand, double mu, double* rval){

    double J = 1.0;
    *rval = 0.;

    BGenerate(sqrts,rand,&J);

    OLP::Arguments Args;
    Args.mu_ren = mu;
    Args.RVal = rval;
    Args.NExt = NextV;
    Args.Order = "NLO";    
    Args.Momenta = Momenta; 
    
####SubProcVirt####

    else{
        std::cout << "Error: Coupling power: "<<cp<<" not found in process"<<std::endl;
        throw "Unrecognized Coupling";
    }
    
    *rval *= J;
}