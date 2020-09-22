#ifndef __INTEGRANDS_H_
#define __INTEGRANDS_H_

#include "nlox_process.h"
####Include Integrands####

class Integrands{
    
    Process Proc;
    int nChannels;
    DipoleStructure** Channels;
    std::unordered_map<std::string,int> ChannelMap;     
            
    public:
        
        Integrands(){

            nChannels = ####nRadiative####;
            Channels = new DipoleStructure* [nChannels];
                        
####Integrand Catalogue####

        }

        ~Integrands(){
            for (int i=0; i<nChannels; i++) delete Channels[i];
            delete [] Channels;
        }

        void Call(std::string IG, std::string CH, std::string CP, double* RV, double* ACC, double* P, double MU,  double x = 1.){

            int Channel = ChannelMap.at(CH);

            if(IG=="Sub"){

                std::cout << "Subtracted function called!" <<std::endl;
                Channels[Channel]->Subtracted(CP,P,MU,RV,ACC);

            }

            else if (IG=="Plu"){

                std::cout << "Plus Distribution function called!" <<std::endl;


            }

            else if (IG=="End"){

                std::cout << "Endpoint function called!" <<std::endl;

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


