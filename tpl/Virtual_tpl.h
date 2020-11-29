#ifndef __VIRTUAL_H__
#define __VIRTUAL_H__

#include "Integrand.h"
####Include Virtuals####

class VirtualIntegrands : public Integrand{
    
    typedef void(VirtualIntegrands::*MemberFunction)(int,std::string,double*,double,double*);
    std::unordered_map<std::string,MemberFunction> VCatalog;
            
    public:

        VirtualStructure** Channels;
        
        VirtualIntegrands(Process * process){

            Proc = process;
            VCatalog.insert({"Virtual",&VirtualIntegrands::Virtual});
            VCatalog.insert({"Born",&VirtualIntegrands::Born});

            nChannels = ####nChannels####;
            Channels = new VirtualStructure* [nChannels];
            
####Virtual Catalogue####

        }

        ~VirtualIntegrands(){
            for (int i=0; i<nChannels; i++) delete Channels[i];
            delete [] Channels;
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

        void Born(int Channel, std::string cp, double* rand, double mu, double* rval){
            Channels[Channel]->Born(cp,rand,rval);
        }

        void Virtual(int Channel, std::string cp, double* rand, double mu, double* rval){
            Channels[Channel]->Virtual(cp,rand,mu,rval);
        }

        void Call(std::string in, std::string ch, std::string cp, double* rand, double mu, double* rval){
            int Channel = ChannelSelect(ch);
            (this->*VCatalog.at(in))(Channel,cp,rand,mu,rval);
        }

};

#endif


