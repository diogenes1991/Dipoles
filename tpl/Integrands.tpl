#ifndef __INTEGRANDS_H_
#define __INTEGRANDS_H_

#include "nlox_process.h"
#define DEBUG 0
####Include Integrands####

struct SubArg{

    std::vector<FourVector> psp;
    std::string cp;
    double* rval;
};

struct PluArg{

    std::string cp;
    std::vector<FourVector> psp;
    double J;
    std::vector<FourVector> psp_1;
    double J_1;
    std::vector<double> mlist;
    double x;
    double* rval;
};

struct EndArg{

    std::string cp;
    std::vector<FourVector> psp;
    double mu;
    double* rval;
};

class Integrands{
    
    int nChannels;
    DipoleStructure** Channels;
    std::unordered_map<std::string,int> ChannelMap;     
            
    public:

        Process Proc;
        
        Integrands(){

            nChannels = ####nRadiative####;
            Channels = new DipoleStructure* [nChannels];
                        
####Integrand Catalogue####

        }

        ~Integrands(){
            for (int i=0; i<nChannels; i++) delete Channels[i];
            delete [] Channels;
        }

        void Call(std::string CH, std::string IG, void* Arg){

            int Channel;
            try {Channel = ChannelMap.at(CH);}
            catch (const std::out_of_range& oor) {
                std::cerr<<"Error: Channel "<<CH<<" not found in Process"<<std::endl;
                std::cout<<"The available Radiative processes are:"<<std::endl;
                for ( auto& x : ChannelMap ) std::cout<<x.first<<" => "<<x.second<<std::endl;
                abort();
            }

            if(IG=="Sub"){

                SubArg* SubArgPtr = static_cast<SubArg*>(Arg);
                Channels[Channel]->Subtracted(SubArgPtr->cp,SubArgPtr->psp,SubArgPtr->rval);

            }

            else if (IG=="Plu"){

                // Channels[Channel]->PlusDistribution(CP,P,MU,RV,ACC);                

            }

            else if (IG=="End"){

                EndArg* EndArgPtr = static_cast<EndArg*>(Arg);
                Channels[Channel]->Endpoint(EndArgPtr->cp,EndArgPtr->psp,EndArgPtr->mu,EndArgPtr->rval);

            }

            else{
                std::cout<<"Integrand type: "<<IG<<" not supported, the suported types are: "<<std::endl;
                std::cout<<" - Sub for the Subtracted Integrand"<<std::endl;
                std::cout<<" - Plu for the Plus Distribution Integrand"<<std::endl;
                std::cout<<" - End for the Endpoint Integrand"<<std::endl;
                abort();

            }
        }

};

#endif


