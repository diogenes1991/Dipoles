#ifndef __GSL_INTEGRATOR_H__
#define __GSL_INTEGRATOR_H__

#include "Montecarlo_Integrator.h"

#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

class GSL_Integrator : public Montecarlo_Integrator{
    
    const gsl_rng_type * t = gsl_rng_default;
    gsl_rng * R = gsl_rng_alloc(t);
    
    double *max,*min;
    const std::string Available[3] = {"Plain","Vegas","Miser"};
    
    double(*Integrand)(double*,size_t,void*);

    public:

        GSL_Integrator(double(*integrand)(double*,size_t,void*), size_t dim){
            Integrand = integrand;
            Dimension = dim;
            min = new double[Dimension];
            max = new double[Dimension];
            for(size_t i =0;i<Dimension;i++){max[i]=1.0;min[i]=0.0;}
            gsl_rng_env_setup ();
        }

        ~GSL_Integrator(){
            delete [] max;
            delete [] min;
            gsl_rng_free(R);
        }

        void Integrate(Specifications * mc_specs, double* result, double* error){
            
            double res,err;
            gsl_monte_function F = {Integrand,Dimension,mc_specs->Params};
            std::string METHOD;
            if(mc_specs->Method=="") METHOD = "Plain";
            else METHOD = mc_specs->Method;

            if(METHOD=="Plain"){
                std::cout<<"###################################################################"<<std::endl;
                std::cout<<"                                                                   "<<std::endl;
                std::cout<<"                   GSL-Plain integration initiated                 "<<std::endl;
                std::cout<<"                                                                   "<<std::endl;
                gsl_monte_plain_state * s = gsl_monte_plain_alloc(Dimension);
                size_t Calls = 0;
                size_t Evals = mc_specs->NStart;
                size_t NCalls = mc_specs->MaxEval;
                size_t Iter = 0;
                while(Calls<NCalls){
                    gsl_monte_plain_integrate (&F,min,max,Dimension,Evals,R,s,&res,&err);
                    Iter += 1;
                    Calls += Evals;
                    Evals += mc_specs->NIncrease;
                    std::cout<<"Iteration["<<Iter<<"] "<<Calls<<" evaluations so far"<<std::endl;
                    std::cout<<res<<" +/- "<<err<<std::endl<<std::endl;
                    if(std::abs(err/res)<mc_specs->RelErr) break;
                }
                gsl_monte_plain_free(s);
                std::cout<<"                                                                   "<<std::endl;
                std::cout<<"###################################################################"<<std::endl;
            }
            
            else if(METHOD=="Vegas"){
                std::cout<<"###################################################################"<<std::endl;
                std::cout<<"                                                                   "<<std::endl;
                std::cout<<"                   GSL-Vegas integration initiated                 "<<std::endl;
                std::cout<<"                                                                   "<<std::endl;
                gsl_monte_vegas_state * s = gsl_monte_vegas_alloc(Dimension);
                size_t Calls = 0;
                size_t Evals = mc_specs->NStart;
                size_t NCalls = mc_specs->MaxEval;
                size_t Iter = 0;
                while(Calls<NCalls){
                    gsl_monte_vegas_integrate (&F,min,max,Dimension,Evals,R,s,&res,&err);
                    Iter += 1;
                    Calls += Evals;
                    Evals += mc_specs->NIncrease;
                    std::cout<<"Iteration["<<Iter<<"] "<<Calls<<" evaluations so far"<<std::endl;
                    std::cout<<res<<" +/- "<<err<<" ("<<gsl_monte_vegas_chisq(s)<<" chisq/dof)"<<std::endl<<std::endl;
                    if(std::abs(err/res)<mc_specs->RelErr) break;
                }
                gsl_monte_vegas_free(s);
                std::cout<<"                                                                   "<<std::endl;
                std::cout<<"###################################################################"<<std::endl;
            }
            
            else if(METHOD=="Miser"){
                std::cout<<"###################################################################"<<std::endl;
                std::cout<<"                                                                   "<<std::endl;
                std::cout<<"                   GSL-Miser integration initiated                 "<<std::endl;
                std::cout<<"                                                                   "<<std::endl;
                gsl_monte_miser_state *s = gsl_monte_miser_alloc(Dimension);
                size_t Calls = 0;
                size_t Evals = mc_specs->NStart;
                size_t NCalls = mc_specs->MaxEval;
                size_t Iter = 0;
                while(Calls<NCalls){
                    gsl_monte_miser_integrate (&F,min,max,Dimension,Evals,R,s,&res,&err);
                    Iter += 1;
                    Calls += Evals;
                    Evals += mc_specs->NIncrease;
                    std::cout<<"Iteration["<<Iter<<"] "<<Calls<<" evaluations so far"<<std::endl;
                    std::cout<<res<<" +/- "<<err<<std::endl<<std::endl;
                    if(std::abs(err/res)<mc_specs->RelErr) break;
                }
                gsl_monte_miser_free(s);
                std::cout<<"                                                                   "<<std::endl;
                std::cout<<"###################################################################"<<std::endl;
            }
            
            else{
                std::cout<<"Error: The integrator routine you requested "<<METHOD<<" is not available, the available GSL routines are:"<<std::endl;
                for(auto a : Available) std::cout<<"    - "<<a<<std::endl;
                throw "Unrecognized Montecarlo Routine";
            }

            *result = res;
            *error  = err; 
        }     
};

#endif
