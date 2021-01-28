#ifndef __NLOX_OLP_H__
#define __NLOX_OLP_H__

#include "OLP.h"
#include "nlox_process.h"

class NLOX_OLP : public OLP{

    // void evaluate_alpha(int i, std::string type, std::string cp, double* pp, int next, double mu, double* rval2, double* acc)
    // void evaluate_alpha_sc(int i, std::string type, std::string cp, double* pp, int next, int pn, double** Bmunu)
    // void evaluate_alpha_cc(int i, std::string type, std::string cp, double* pp, int next, double* Bij)

    public:

        Process * Proc;
        std::unordered_map<std::string,std::string> OrderMap({{"LO","tree_tree"},{"NLO","tree_loop"}});

        NLOX_OLP(){
            Proc = new Process();

####Define Channels####
        }

        ~NLOX_OLP(){
            delete Proc;
        }

        void Evaluate(Arguments * arg){
            int Channel = SelectSubProcess(arg->SubProc);
            std::string Type = OrderMap.at(arg->Order);
            std::string Coupling = "as"+std::to_string(arg->CouplingPower[0])+"ae"+std::to_string(arg->CouplingPower[1]);
            int NExt = arg->NExt;
            double Momenta[5*NExt];
            double Mu = arg->mu_ren;
            double acc;
            double RVal[3];
            for(int i=0;i<NExt;i++){
                Momenta[5*i+0]=arg->Momenta[i].p0;
                Momenta[5*i+1]=arg->Momenta[i].p1;
                Momenta[5*i+2]=arg->Momenta[i].p2;
                Momenta[5*i+3]=arg->Momenta[i].p3;
                Momenta[5*i+4]=0.0;
            } 

            Proc->evaluate_alpha(Channel,Type,Coupling,Momenta,NExt,Mu,RVal,&acc);
            *arg->RVal = RVal[2];
            
        }


        void Evaluate_CC(Arguments * arg){
            
        }

        void Evaluate_SC(Arguments * arg){
            
        }
};

#endif
