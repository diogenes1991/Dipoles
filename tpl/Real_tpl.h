#ifndef __REAL_H__
#define __REAL_H__

#include "Integrand.h"
####Include Integrands####

class RealIntegrands : public Integrand{
    
    typedef void(RealIntegrands::*MemberFunction)(int,std::string,double,double*,double,double*);
    std::unordered_map<std::string,MemberFunction> RCatalog;

    public:

        DipoleStructure** Channels;

        RealIntegrands(Process * process){

            Proc = process;
            RCatalog.insert({"Subtracted",&RealIntegrands::Subtracted});
            RCatalog.insert({"PlusDistribution",&RealIntegrands::PlusDistribution});
            RCatalog.insert({"Endpoint",&RealIntegrands::Endpoint});    

            nChannels = ####nRadiative####;
            Channels = new DipoleStructure* [nChannels];
                        
####Integrand Catalogue####

        }

        ~RealIntegrands(){
            for (int i=0; i<nChannels; i++) delete Channels[i];
            delete [] Channels;
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

        void Subtracted(int Channel, std::string cp, double sqrts, double* rand, double mu, double* rval){
            Channels[Channel]->Subtracted(cp,sqrts,rand,mu,rval);
        }

        void PlusDistribution(int Channel, std::string cp, double sqrts, double* rand, double mu, double* rval){
            Channels[Channel]->PlusDistribution(cp,sqrts,rand,mu,rval);
        }      

        void Endpoint(int Channel, std::string cp, double sqrts, double* rand, double mu, double* rval){
            Channels[Channel]->Endpoint(cp,sqrts,rand,mu,rval);
        }

        void Call(std::string in, std::string ch, std::string cp, double sqrts, double* rand, double mu, double* rval){
            int Channel = ChannelSelect(ch);
            (this->*RCatalog.at(in))(Channel,cp,sqrts,rand,mu,rval);
        }

};

#endif


