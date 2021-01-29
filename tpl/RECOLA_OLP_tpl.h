#ifndef __RECOLA_OLP_H__
#define __RECOLA_OLP_H__

#include "OLP.h"
#include "recola.hpp"

class RECOLA_OLP : public OLP{

    // Recola::compute_process_rcl(int i,double* p,std::string Order);
    // Recola::get_squared_amplitude_rcl(int i,int pow,std::string Order,double RVal);

    public:

        RECOLA_OLP(Model * Mod){

            model = Mod;
####Define Channels####
            Recola::generate_processes_rcl();
        }

        ~RECOLA_OLP(){
            Recola::reset_recola_rcl();
        }

        void UpdateParameters(){
            
        }

        void Evaluate(Arguments * arg){
            int Channel = SelectSubProcess(arg->SubProc);
            std::string Order = arg->Order;
            int NExt = arg->NExt;
            double Momenta[NExt][4];
            double RVal;
            int Pow = arg->CouplingPower[0];
            for(int i=0;i<NExt;i++){
                Momenta[i][0]=arg->Momenta[i].p0;
                Momenta[i][1]=arg->Momenta[i].p1;
                Momenta[i][2]=arg->Momenta[i].p2;
                Momenta[i][3]=arg->Momenta[i].p3;
            }
            
            Recola::compute_process_rcl(Channel,Momenta,Order);
            Recola::get_squared_amplitude_rcl(Channel,Pow,Order,RVal);

            *arg->RVal = RVal;
            
        }


        void Evaluate_CC(Arguments * arg){
            
        }

        void Evaluate_SC(Arguments * arg){
            
        }
};

#endif
