#include "####SubProcName####.h"

####SubProcName####::####SubProcName####(Process& process){

    
    Proc = &process; 

####SubProcConst####

}

void ####SubProcName####::SetECM(double ECM){
    sqrts=ECM;
    double rIn[2] = {0.0,0.0};
    SetInMom(rIn);
}

void ####SubProcName####::SetFiMom(double* rFi, double* J){
    FourVector PCM(sqrts,0,0,0);
    FourVector PFi[Next-2];
    double mFi[Next-2];
    double dummy = 1.0;
    for (int i=0;i<Next-2;i++) mFi[i] = masses[i+2];
    Recursive_PSP(PCM,Next-2,PFi,mFi,rFi,dummy);
    for (int i=0;i<Next-2;i++) momenta[i+2] = PFi[i];
    *J = dummy;    
}

void ####SubProcName####::SetInMom(double* rIn){
    FourVector PCM(sqrts,0,0,0);
    FourVector PIn[2];
    double mIn[2]={masses[0],masses[1]};
    double dummy = 1.0;
    Recursive_PSP(PCM,2,PIn,mIn,rIn,dummy);
    momenta[0] = PIn[0];
    momenta[1] = PIn[1];
}

void ####SubProcName####::Subtracted(std::string cp, double* rand, double* rval){

    int i;
    double radiative[3];
    double EWKFac = 4*M_PI*(Proc->pc.alpha_e);
    double QCDFac = 4*M_PI*(Proc->pc.alpha_s);
    double DipFac = 1.0;
    double mu = 0.0;
    double J = 1.0;
    double acc;
    *rval = 0.;

    SetFiMom(rand,&J);
    std::vector<FourVector> p,p_tilde;
    for(int j=0;j<Next;j++){
        p.push_back(momenta[j]);
        p_tilde.push_back(momenta[j]);
    }
    
####SubProcSub####

   else{
        std::cout << "Error: Coupling power: "<<cp<<" not found in process"<<std::endl;
        abort();
   }

}

void ####SubProcName####::PlusDistribution(std::string cp, double* rand, double mu, double* rval){

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

void ####SubProcName####::Endpoint(std::string cp, double* rand, double mu, double* rval){

    int i;
    double EWKFac = Proc->pc.alpha_e/(2*M_PI);
    double QCDFac = Proc->pc.alpha_s/(2*M_PI);
    double DipFac = 1.0;
    double J = 1.0;
    rval[0] = 0.;
    rval[1] = 0.;
    rval[2] = 0.;
    double acc;
    
    SetFiMom(rand,&J);
    std::vector<FourVector> p;
    for(int j=0;j<Next;j++) p.push_back(momenta[j]);

####SubProcEnd####

    else{
        std::cout << "Error: Coupling power: "<<cp<<" not found in process"<<std::endl;
        abort();
   }
}