
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
    Integrand->setECM(2000);    

    std::string Integral = string(argv[1]);
    std::string Channel = string(argv[2]);
    std::string Coupling = string(argv[3]);

    std::cout<<"Calling "<<Integral<<" integrand for the channel "<<Channel<<" at coupling power "<<Coupling<<std::endl;

    double mu =500;


    if (Integral=="Sub"){

    int step = 100;

    for (float i=1;i<step;i+=1){
    
        double rFi[5] = {1.0-i/step,0.8328481,0.58816431,0.9137812,0.5991822};
        double rval;
    
        Integrand->Subtracted(Channel,Coupling,rFi,&rval);

        // std::cout << "Subtracted Double Pole  = "<<0.0<<std::endl;
        // std::cout << "Subtracted Single Pole  = "<<0.0<<std::endl;
        std::cout << "Subtracted    Finite    = "<<rval<<std::endl;
    }

    }

    if (Integral=="End"){

        double rFi[2] = {0.321837231984612,0.716746127852};
        double rval[3];

        Integrand->Endpoint(Channel,Coupling,rFi,mu,rval);

        std::cout << "Endpoint Double Pole  = "<<rval[0]<<std::endl;
        std::cout << "Endpoint Single Pole  = "<<rval[1]<<std::endl;
        std::cout << "Endpoint    Finite    = "<<rval[2]<<std::endl;

        // int iv = 0;
        // std::string types = typ;
        // std::string cps = Coupling;
        // int next = 4;
        // // Evaluate
        // Integrand->Proc.evaluate_alpha(iv, types, cps, p, next, mu, rval, &acc);

        // std::cout << "Virtual Double Pole  = "<<rval[0]<<std::endl;
        // std::cout << "Virtual Single Pole  = "<<rval[1]<<std::endl;
        // std::cout << "Virtual    Finite    = "<<rval[2]<<std::endl;
        // std::cout << "Virtual      acc     = "<<acc<<std::endl;

    }

    

    
    return 0;
}