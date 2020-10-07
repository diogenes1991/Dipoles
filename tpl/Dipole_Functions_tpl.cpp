#include "####SubProcName####.h"

####SubProcName####::####SubProcName####(Process& process){

    
    Proc = &process; 

####SubProcConst####

}

void ####SubProcName####::SetECM(double ECM){
    P.p0=ECM;
    double rIn[2] = {0.0,0.0};
    SetInMom(rIn);
}

void ####SubProcName####::SetFiMom(double* rFi, double* J){
    FourVector PFi[Next-2];
    double mFi[Next-2];
    double dummy = 1.0;
    for (int i=0;i<Next-2;i++) mFi[i] = Masses[i+2];
    Recursive_PSP(P,Next-2,PFi,mFi,rFi,dummy);
    for (int i=0;i<Next-2;i++) Momenta[i+2] = PFi[i];
    *J = dummy;    
}

void ####SubProcName####::SetInMom(double* rIn){
    FourVector PIn[2];
    double mIn[2]={Masses[0],Masses[1]};
    double dummy = 1.0;
    Recursive_PSP(P,2,PIn,mIn,rIn,dummy);
    Momenta[0] = PIn[0];
    Momenta[1] = PIn[1];
}

void ####SubProcName####::SetFiMom(int BornNum, double* rFi, double* J){
    FourVector PFi[Next-3];
    double mFi[Next-3];
    double dummy = 1.0;
    for (int i=0;i<Next-1;i++) mFi[i] = BornMasses[BornNum][i+2];
    Recursive_PSP(P,Next-1,PFi,mFi,rFi,dummy);
    for (int i=0;i<Next-1;i++) BornMomenta[i+2] = PFi[i];
    *J = dummy;    
}

void ####SubProcName####::SetInMom(int BornNum){
    FourVector PIn[2];
    double mIn[2]={BornMasses[BornNum][0],BornMasses[BornNum][1]};
    double rIn[2]={0.0,0.0};
    double dummy = 1.0;
    Recursive_PSP(P,2,PIn,mIn,rIn,dummy);
    BornMomenta[0] = PIn[0];
    BornMomenta[1] = PIn[1];
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
        p.push_back(Momenta[j]);
        p_tilde.push_back(Momenta[j]);
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

    int i,j;
    double EWKFac = Proc->pc.alpha_e/(2*M_PI);
    double QCDFac = Proc->pc.alpha_s/(2*M_PI);
    double DipFac = 1.0;
    double J = 1.0;
    double acc,Invariant;
    std::vector<FourVector> p;
    rval[0] = 0.;
    rval[1] = 0.;
    rval[2] = 0.;
    
####SubProcEnd####

    else{
        std::cout << "Error: Coupling power: "<<cp<<" not found in process"<<std::endl;
        abort();
   }
}