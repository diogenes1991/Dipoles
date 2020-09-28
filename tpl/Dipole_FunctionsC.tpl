#include "####SubProcName####.h"

####SubProcName####::####SubProcName####(Process& process){

    
    Proc = &process; 

####SubProcConst####

}

void ####SubProcName####::Subtracted(std::string cp, std::vector<FourVector> p, double* rval){

    int i;
    double radiative[3];
    std::vector<FourVector> p_tilde = p;
    double EWKFac = 4*M_PI*(Proc->pc.alpha_e);
    double QCDFac = 4*M_PI*(Proc->pc.alpha_s);
    double DipFac = 1.0;
    double mu = 0.0;
    double acc;
    *rval = 0.;

####SubProcSub####

   else{
        std::cout << "Error: Coupling power: "<<cp<<" not found in process"<<std::endl;
        abort();
   }

}

void ####SubProcName####::PlusDistribution(std::string cp, std::vector<FourVector> p, double mu, double* rval, double* acc){

    int i;
    double radiative[3];
    std::vector<FourVector> p_tilde = p;
    double EWKFac = Proc->pc.alpha_e/(2*M_PI);
    double QCDFac = Proc->pc.alpha_s/(2*M_PI);
    *rval = 0.;

####SubProcPlu####

    else{
        std::cout << "Error: Coupling power: "<<cp<<" not found in process"<<std::endl;
        abort();
   }

}

void ####SubProcName####::Endpoint(std::string cp, std::vector<FourVector> p, double mu, double* rval){

    int i;
    double MEm,MSp;
    double EWKFac = Proc->pc.alpha_e/(2*M_PI);
    double QCDFac = Proc->pc.alpha_s/(2*M_PI);
    double DipFac = 1.0;
    rval[0] = 0.;
    rval[1] = 0.;
    rval[2] = 0.;
    double acc;


####SubProcEnd####

    else{
        std::cout << "Error: Coupling power: "<<cp<<" not found in process"<<std::endl;
        abort();
   }
}