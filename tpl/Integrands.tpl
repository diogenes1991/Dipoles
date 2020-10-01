#ifndef __INTEGRANDS_H_
#define __INTEGRANDS_H_

#include "nlox_process.h"
####Include Integrands####

class Integrands{
    
    int nChannels;
    int nExternal;
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

        int ChannelSelect(std::string CH){
            int Channel;
            try {Channel = ChannelMap.at(CH);}
            catch (const std::out_of_range& oor) {
                std::cerr<<"Error: Channel "<<CH<<" not found in Process"<<std::endl;
                std::cout<<"The available Radiative processes are:"<<std::endl;
                for ( auto& x : ChannelMap ) std::cout<<x.first<<" => "<<x.second<<std::endl;
                abort();
            }
            return Channel;
        }

        void setECM(double sqrts){
            for ( int i=0;i<nChannels;i++) Channels[i]->setECM(sqrts);
        }

        void Subtracted(std::string ch, std::string cp, double* rand, double* rval){
            int Channel = ChannelSelect(ch);
            Channels[Channel]->Subtracted(cp,rand,rval);
        }

        void PlusDistribution(std::string ch, std::string cp, double* rand, double mu, double* rval){

        }      

        void Endpoint(std::string ch, std::string cp, double* rand, double mu, double* rval){
            int Channel = ChannelSelect(ch);
            Channels[Channel]->Endpoint(cp,rand,mu,rval);
        }
};

#endif


