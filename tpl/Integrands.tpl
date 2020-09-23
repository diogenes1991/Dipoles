#ifndef __INTEGRANDS_H_
#define __INTEGRANDS_H_

#include "nlox_process.h"
#define DEBUG 0
####Include Integrands####

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

        void Call(std::string IG, std::string CH, std::string CP, double* RV, double* ACC, std::vector<FourVector> P, double MU,  double x = 1.){

            int Channel = ChannelMap.at(CH);

            if(IG=="Sub"){

                Channels[Channel]->Subtracted(CP,P,MU,RV,ACC);

            }

            else if (IG=="Plu"){

                Channels[Channel]->PlusDistribution(CP,P,MU,RV,ACC);                

            }

            else if (IG=="End"){

                Channels[Channel]->Endpoint(CP,P,MU,RV,ACC);

            }

            else{
                std::cout<<"Integrand type: "<<IG<<" not supported, the suported types are: "<<std::endl;
                std::cout<<" - Sub for the Subtracted Integrand"<<std::endl;
                std::cout<<" - Plu for the Plus Distribution Integrand"<<std::endl;
                std::cout<<" - End for the Endpoint Integrand"<<std::endl;

            }
        }

};

#endif


