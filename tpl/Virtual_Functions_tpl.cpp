#include "####SubProcName####.h"

####SubProcName####::####SubProcName####(Process& process){

    Proc = &process; 
    nPar = NextV; 
    BornMomenta = new FourVector[NextV];
    BornMasses = new double[NextV];
    BornPID = new int[NextV];

####SubProcConst####

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