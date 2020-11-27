#ifndef __VIRTUAL_H_
#define __VIRTUAL_H_

#include "nlox_process.h"
####Include Virtuals####

class VirtualIntegrands{
    
    int nChannels;
    VirtualStructure** Channels;
    std::unordered_map<std::string,int> ChannelMap;     
            
    public:

        Process * Proc;
        
        VirtualIntegrands(Process * process){

            Proc = process;

            nChannels = ####nChannels####;
            Channels = new VirtualStructure* [nChannels];
            
####Virtual Catalogue####

        }

        ~VirtualIntegrands(){
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
            for ( int i=0;i<nChannels;i++) Channels[i]->SetECM(sqrts);
        }

        void GetMomenta(std::string ch, FourVector* p){
            int Channel = ChannelSelect(ch);
            Channels[Channel]->GetMomenta(p);
        }

        void GetMasses(std::string ch, double* m){
            int Channel = ChannelSelect(ch);
            Channels[Channel]->GetMasses(m);
        }

        void GetPID(std::string ch, int* pid){
            int Channel = ChannelSelect(ch);
            Channels[Channel]->GetPID(pid);
        }

        void Born(std::string ch, std::string cp, double* rand, double* rval){
            int Channel = ChannelSelect(ch);
            Channels[Channel]->Born(cp,rand,rval);
        }

        void Virtual(std::string ch, std::string cp, double* rand, double mu, double* rval){
            int Channel = ChannelSelect(ch);
            Channels[Channel]->Virtual(cp,rand,mu,rval);
        }


};

#endif


