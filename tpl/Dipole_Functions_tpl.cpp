#include "####SubProcName####.h"

####SubProcName####::####SubProcName####(Process& process){

    Proc = &process; 
    nPar = NextR;
    Momenta = new FourVector[NextR];
    Masses = new double[NextR];
    PID = new int[NextR];
    BornMomenta = new FourVector[NextR-1];

####SubProcConst####

}

void ####SubProcName####::Subtracted(std::string cp, double sqrts, double* rand, double mu, double* rval){

    int i;
    double radiative[3];
    double EWKFac = 4*M_PI*(Proc->pc.alpha_e);
    double QCDFac = 4*M_PI*(Proc->pc.alpha_s);
    double DipFac = 1.0;
    double J = 1.0;
    double acc;
    *rval = 0.;

    SGenerate(sqrts,rand,&J);
    std::vector<FourVector> p,p_tilde;
    double pp[5*NextR],pp_tilde[5*NextR];
    for(int j=0;j<NextR;j++){
        p.push_back(Momenta[j]);
    }
    
####SubProcSub####

   else{
        std::cout << "Error: Coupling power: "<<cp<<" not found in process"<<std::endl;
        abort();
   }

   *rval *= J;

}

void ####SubProcName####::PlusDistribution(std::string cp, double sqrts, double* rand, double mu, double* rval){

    int i;
    double Ix,I1;
    double EWKFac = Proc->pc.alpha_e/(2*M_PI);
    double QCDFac = Proc->pc.alpha_s/(2*M_PI);
    double DipFac = 1.0;
    *rval = 0.;

####SubProcPlu####

    else{
        std::cout << "Error: Coupling power: "<<cp<<" not found in process"<<std::endl;
        abort();
   }

}

void ####SubProcName####::Endpoint(std::string cp, double sqrts, double* rand, double mu, double* rval){

    int i,j;
    double EWKFac = Proc->pc.alpha_e/(2*M_PI);
    double QCDFac = Proc->pc.alpha_s/(2*M_PI);
    double DipFac = 1.0;
    double J = 1.0;
    double acc,Invariant;
    std::vector<FourVector> p;
    double pp[5*NextR];
    rval[0] = 0.;
    rval[1] = 0.;
    rval[2] = 0.;
    
####SubProcEnd####

    else{
        std::cout << "Error: Coupling power: "<<cp<<" not found in process"<<std::endl;
        abort();
   }
}