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
            UpdateParameters();
            Recola::generate_processes_rcl();
        }

        ~RECOLA_OLP(){
            Recola::reset_recola_rcl();
        }

        void UpdateParameters(){

            //
            //  Masses: RECOLA
            //  
            
            if(model->UseCMScheme) Recola::set_complex_mass_scheme_rcl();

            if(model->d.Mass){
                Recola::set_pole_mass_down_rcl(model->d.Mass);
                Recola::unset_light_down_rcl();
            }
            else Recola::set_light_down_rcl();
            
            if(model->u.Mass){
                Recola::set_pole_mass_up_rcl(model->u.Mass);
                Recola::unset_light_up_rcl();
            }
            else Recola::set_light_up_rcl();

            if(model->s.Mass){
                Recola::set_pole_mass_strange_rcl(model->s.Mass);
                Recola::unset_light_strange_rcl();
            }
            else Recola::set_light_strange_rcl();
            
            if(model->c.Mass){
                Recola::set_pole_mass_charm_rcl(model->c.Mass,model->c.Width);
                Recola::unset_light_charm_rcl();
            }
            else Recola::set_light_charm_rcl();
            
            if(model->b.Mass){
                Recola::set_pole_mass_bottom_rcl(model->b.Mass,model->b.Width);
                Recola::unset_light_bottom_rcl();
            }
            else Recola::set_light_bottom_rcl();
            
            Recola::set_pole_mass_top_rcl(model->t.Mass,model->t.Width);
            Recola::set_pole_mass_w_rcl(model->Wp.Mass,model->Wp.Width);
            Recola::set_pole_mass_z_rcl(model->Z.Mass,model->Z.Width);
            Recola::set_pole_mass_h_rcl(model->h.Mass,model->h.Width);
            
            //
            //  QCD Scheme
            //

            Recola::set_alphas_rcl(model->alpha_s,1000,model->NLF);
            
            //
            //  EWK Scheme
            //
            
            if(model->UseGMuScheme) Recola::use_gfermi_scheme_and_set_gfermi_rcl(model->GFermi);
            else{
                Recola::use_alpha0_scheme_rcl(model->alpha_e);
                Recola::set_delta_uv_rcl(1.0);
            }

            Recola::set_dynamic_settings_rcl(1);
            Recola::set_mu_uv_rcl(1000);
            Recola::set_mu_ir_rcl(1000);
            Recola::set_delta_ir_rcl(0,0);
               
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
