#include "####SubProcName####.h"

####SubProcName####::####SubProcName####(Process& process){

    
    Proc = &process; 

####SubProcConst####

}

void ####SubProcName####::Subtracted(std::string cp, std::vector<FourVector> p, double mu, double* rval, double* acc){

    int i;
    double radiative[3];
    std::vector<FourVector> p_tilde = p;
    double EWKFac = 4*M_PI*(Proc->pc.alpha_e);
    double QCDFac = 4*M_PI*(Proc->pc.alpha_s);
    rval[0] = 0.;
    rval[1] = 0.;
    rval[2] = 0.;
    *acc = 0.;

####SubProcSub####

   else{
        std::cout << "Error: Coupling power: "<<cp<<" not found in process"<<std::endl;
        abort();
   }


}

void ####SubProcName####::PlusDistribution(std::string cp, std::vector<FourVector> p, double mu, double* rval, double* acc){

####SubProcPlu####

}

void ####SubProcName####::Endpoint(std::string cp, std::vector<FourVector> p, double mu, double* rval, double* acc){

####SubProcEnd####

}