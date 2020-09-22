
#include "Integrands.h"

int main(int argc, char* argv[]){

    
    Integrands* Integrand;
    Integrand = new Integrands();

    std::string Integral = string(argv[1]);
    std::string Channel = string(argv[2]);
    std::string Coupling = string(argv[3]);

    double rval[3];
    double acc;
    double mu =500;

    FourVector P(1000,0,0,0);

    FourVector PIn[2];
    double mIn[2] = {0,4.2};
    double rIn[2] = {0.,0.};
    double JIn = 1.;
    Recursive_PSP(P,2,PIn,mIn,rIn,JIn);

    FourVector PFi[3];
    double mFi[3] = {4.2,0,91.1876};
    double rFi[5] = {0.381283921,0.8328481,0.58816431,0.9137812,0.5991822};
    double JFi = 1.;
    Recursive_PSP(P,3,PFi,mFi,rFi,JFi);

    std::vector<FourVector> p;
    p.push_back(PIn[0]);
    p.push_back(PIn[1]);
    p.push_back(PFi[0]);
    p.push_back(PFi[1]);
    p.push_back(PFi[2]);

    std::cout<<"Calling "<<Integral<<" integrand for the channel "<<Channel<<" at coupling power "<<Coupling<<std::endl;

    Integrand->Call(Integral,Channel,Coupling,rval,&acc,p,mu);

    std::cout << "Subtracted Double Pole  = "<<rval[0]<<std::endl;
    std::cout << "Subtracted Single Pole  = "<<rval[1]<<std::endl;
    std::cout << "Subtracted    Finite    = "<<rval[2]<<std::endl;

    return 0;
}