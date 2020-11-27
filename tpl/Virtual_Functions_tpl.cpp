#include "####SubProcName####.h"

####SubProcName####::####SubProcName####(Process& process){

    
    Proc = &process; 

####SubProcConst####

}

void ####SubProcName####::SetECM(double ECM){
    P.p0=ECM;
    SetInMom();
}

void ####SubProcName####::SetFiMom(double* rFi, double* J){
    FourVector PFi[NextV-2];
    double mFi[NextV-2];
    double dummy = 1.0;
    for (int i=0;i<NextV-2;i++) mFi[i] = BornMasses[i+2];
    Recursive_PSP(P,NextV-2,PFi,mFi,rFi,dummy);
    for (int i=0;i<NextV-2;i++) BornMomenta[i+2] = PFi[i];
    *J = dummy;    
}

void ####SubProcName####::SetInMom(){
    FourVector PIn[2];
    double mIn[2]={BornMasses[0],BornMasses[1]};
    double rIn[2]={0.0,0.0};
    double dummy = 1.0;
    Recursive_PSP(P,2,PIn,mIn,rIn,dummy);
    BornMomenta[0] = PIn[0];
    BornMomenta[1] = PIn[1];
}

void ####SubProcName####::GetMomenta(FourVector* p){

    for(int i=0;i<NextV;i++)p[i]=BornMomenta[i];
}

void ####SubProcName####::GetMasses(double* m){

    for(int i=0;i<NextV;i++)m[i]=BornMasses[i];
}

void ####SubProcName####::GetPID(int* pid){

    for(int i=0;i<NextV;i++)pid[i]=BornPID[i];
}

void ####SubProcName####::Born(std::string cp, double* rand, double* rval){

    double J = 1.0;
    double born[3];
    double mu = 0.0;
    double acc;
    *rval = 0.;

    SetFiMom(rand,&J);
    double pp[5*NextV];
    for(int i=0;i<NextV;i++){pp[5*i+0]=BornMomenta[i].p0;pp[5*i+1]=BornMomenta[i].p1;pp[5*i+2]=BornMomenta[i].p2;pp[5*i+3]=BornMomenta[i].p3;pp[5*i+4]=0.0;}
    
####SubProcBorn####

    else{
        std::cout << "Error: Coupling power: "<<cp<<" not found in process"<<std::endl;
        abort();
    }
    
    *rval = J*born[2];
}

void ####SubProcName####::Virtual(std::string cp, double* rand, double mu, double* rval){

    double J = 1.0;
    double virt[3];
    double acc;
    *rval = 0.;

    SetFiMom(rand,&J);
    double pp[5*NextV];
    for(int i=0;i<NextV;i++){pp[5*i+0]=BornMomenta[i].p0;pp[5*i+1]=BornMomenta[i].p1;pp[5*i+2]=BornMomenta[i].p2;pp[5*i+3]=BornMomenta[i].p3;pp[5*i+4]=0.0;}
    
####SubProcVirt####

    else{
        std::cout << "Error: Coupling power: "<<cp<<" not found in process"<<std::endl;
        abort();
    }
    
    *rval = J*virt[2];
}