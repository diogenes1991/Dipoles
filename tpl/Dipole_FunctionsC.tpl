#include "####SubProcName####.h"

####SubProcName####::####SubProcName####(Process& process){

    
    Proc = &process; 

####SubProcConst####

}

void ####SubProcName####::setECM(double sqrts){
    PCM.p0=sqrts;
    double mIn[2]={masses[0],masses[1]};
    double rIn[2]={0.0,0.0};
    double dummy = 1.0;
    Recursive_PSP(PCM,2,PIn,mIn,rIn,dummy);
}

void ####SubProcName####::Subtracted(std::string cp, double* rand, double* rval){

    int i;
    double radiative[3];
    double EWKFac = 4*M_PI*(Proc->pc.alpha_e);
    double QCDFac = 4*M_PI*(Proc->pc.alpha_s);
    double DipFac = 1.0;
    double mu = 0.0;
    double acc;
    *rval = 0.;

    FourVector pFi[Next-2];
    double mFi[Next-2];
    for(int i=2;i<Next;i++) mFi[i-2] = masses[i];
    Recursive_PSP(PCM,Next-2,pFi,mFi,rand,DipFac);
    std::vector<FourVector> p{PIn[0],PIn[1]};
    for (int i=0;i<Next-2;i++) p.push_back(pFi[i]);
    std::vector<FourVector> p_tilde = p;
    
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
    rval[0] = 0.;
    rval[1] = 0.;
    rval[2] = 0.;
    double acc;
    
    FourVector pFi[Next-2];
    double mFi[Next-2];
    for(int i=2;i<Next-1;i++) mFi[i-2] = masses[i];
    Recursive_PSP(PCM,Next-3,pFi,mFi,rand,DipFac);
    std::vector<FourVector> p{PIn[0],PIn[1]};
    for (int i=0;i<Next-2;i++) p.push_back(pFi[i]);
    std::vector<FourVector> p_tilde = p;
    
####SubProcEnd####

    else{
        std::cout << "Error: Coupling power: "<<cp<<" not found in process"<<std::endl;
        abort();
   }
}