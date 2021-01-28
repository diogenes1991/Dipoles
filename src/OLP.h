#ifndef __OLP_H__
#define __OLP_H__

#include <unordered_map>
#include "Four_Vector.h"
#include "Model.h"

class OLP{

    public:

        std::unordered_map<std::string,int> ChannelIndex;

        struct Arguments{
            std::string SubProc;
            std::string Order;
            int CouplingPower[2];
            double mu_ren;
            FVector * Momenta;
            int NExt;
            double * Rval;
        };

        int SelectSubProcess(std::string SP){
            int SubProcess;
            try {SubProcess = ChannelIndex.at(SP);}
            catch (const std::out_of_range& oor) {
                std::cerr<<"Error: Channel "<<SP<<" not found in OLP"<<std::endl;
                std::cout<<"The available Channels are:"<<std::endl;
                for ( auto& x : ChannelIndex ) std::cout<<x.first<<" => "<<x.second<<std::endl;
                abort();
            }
            return SubProcess;
        }

        virtual ~OLP(){};
        virtual void Evaluate(Arguments * args) = 0;
        virtual void Evaluate_CC(Arguments * args) = 0;
        virtual void Evaluate_SC(Arguments * args) = 0;

};

#endif
