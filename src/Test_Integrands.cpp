
#include "Integrands.h"

// 
// There are three integrand types bundled by
// the integrand class "Sub", "Plu" and "End"
// Currently there are no checks on the arguments 
// passed to each adn so they are assumed to be well
// called with the correct arguments, which are:
// 
// Subtracted: - Coupling power passed as an std::string "asXaeY" for the desired X and Y
//             - Phase-Space point passed as an std::vector<FourVector> a.k.a. std::vector<FourVectorT<double>>
//             - Pointer to store result a single double value is necessary 
// 


int main(int argc, char* argv[]){

    if(argc!=4) abort();

    Integrands* Integrand;
    Integrand = new Integrands();

    std::string Integral = string(argv[1]);
    std::string Channel = string(argv[2]);
    std::string Coupling = string(argv[3]);

    std::cout<<"Calling "<<Integral<<" integrand for the channel "<<Channel<<" at coupling power "<<Coupling<<std::endl;

    double mu =500;

    FourVector P(1000,0,0,0);

    FourVector PIn[2];
    double mIn[2] = {0.0,0.0};
    double rIn[2] = {0.0,0.0};
    double JIn = 1.;
    Recursive_PSP(P,2,PIn,mIn,rIn,JIn);


    if (Integral=="Sub"){
        // Prepare a 2->n+1 psp

    FourVector PFi[3];
    double mFi[3] = {91.1876,91.1876,0};
    double rFi[5] = {0.381283921,0.8328481,0.58816431,0.9137812,0.5991822};
    double JFi = 1.;
    Recursive_PSP(P,3,PFi,mFi,rFi,JFi);

    std::vector<FourVector> p;
    p.push_back(PIn[0]);
    p.push_back(PIn[1]);
    p.push_back(PFi[0]);
    p.push_back(PFi[1]);
    p.push_back(PFi[2]);

    double rval;
    SubArg SArg;
    SArg.psp = p;
    SArg.cp = Coupling;
    SArg.rval = &rval;

    Integrand->Call(Channel,Integral,&SArg);

    std::cout << "Subtracted Double Pole  = "<<0.0<<std::endl;
    std::cout << "Subtracted Single Pole  = "<<0.0<<std::endl;
    std::cout << "Subtracted    Finite    = "<<rval<<std::endl;

    }

    if (Integral=="End"){

        // Prepare a 2->n psp
        FourVector PFi[2];
        double mFi[2] = {91.1876,91.1876};
        double rFi[2] = {0.321837231984612,0.716746127852};
        double JFi = 1.0;
        Recursive_PSP(P,2,PFi,mFi,rFi,JFi);

        std::vector<FourVector> p;
        p.push_back(PIn[0]);
        p.push_back(PIn[1]);
        p.push_back(PFi[0]);
        p.push_back(PFi[1]);

        double rval[3];
        EndArg EArg;
        EArg.psp = p;
        EArg.cp = Coupling;
        EArg.mu = mu;
        EArg.rval = rval;

        Integrand->Call(Channel,Integral,&EArg);

        std::cout << "Endpoint Double Pole  = "<<rval[0]<<std::endl;
        std::cout << "Endpoint Single Pole  = "<<rval[1]<<std::endl;
        std::cout << "Endpoint    Finite    = "<<rval[2]<<std::endl;

        
        char typ[] = "tree_loop";
        // char cp[] = Coupling;
        double acc;


        int iv = 0;
        std::string types = typ;
        std::string cps = Coupling;
        int next = 4;
        // Evaluate
        Integrand->Proc.evaluate_alpha(iv, types, cps, p, next, mu, rval, &acc);

        std::cout << "Virtual Double Pole  = "<<rval[0]<<std::endl;
        std::cout << "Virtual Single Pole  = "<<rval[1]<<std::endl;
        std::cout << "Virtual    Finite    = "<<rval[2]<<std::endl;
        std::cout << "Virtual      acc     = "<<acc<<std::endl;


    }

    

    
    return 0;
}