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
        std::unordered_map<std::string,std::string> OrderMap = {{"LO","tree_tree"},{"NLO","tree_loop"}};

        NLOX_OLP(Model * Mod){
            
            model = Mod;
            Proc = new Process();
            UpdateParameters();

####Define Channels####
        }

        ~NLOX_OLP(){
            delete Proc;
        }

        void UpdateParameters(){

            //
            //  Scheme
            //

            Proc->pc.nlf       = model->NLF;
            Proc->pc.nhf       = model->NHF;

            Proc->pc.cmScheme  = model->UseCMScheme;
            Proc->pc.GmuScheme = model->UseGMuScheme;

            //
            //  Masses: Only t and b quarks are allowed 
            //  massive by NLOX
            //

            Proc->pc.set_param("mb",model->b.Mass);
            Proc->pc.set_param("wb",model->b.Width);
            Proc->update_mass("mb");
            Proc->update_mass("wb");

            Proc->pc.set_param("mt",model->t.Mass);
            Proc->pc.set_param("wt",model->t.Width);
            Proc->update_mass("mt");
            Proc->update_mass("wt");

            Proc->pc.set_param("mW",model->Wp.Mass);
            Proc->pc.set_param("wW",model->Wp.Width);
            Proc->update_mass("mW");
            Proc->update_mass("wW");

            Proc->pc.set_param("mZ",model->Z.Mass);
            Proc->pc.set_param("wZ",model->Z.Width);
            Proc->update_mass("mZ");
            Proc->update_mass("wZ");

            Proc->pc.set_param("mH",model->h.Mass);
            Proc->pc.set_param("wH",model->h.Width);
            Proc->update_mass("mH");
            Proc->update_mass("wH");

            //
            //  Couplings
            //

            Proc->pc.set_param("alpha_s",model->alpha_s);
            Proc->pc.set_param("alpha_e",model->alpha_e);
            Proc->pc.set_param("GF",model->GFermi);
            Proc->pc.compute(true,true);

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
