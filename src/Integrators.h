#ifndef _INTEGRATORS_H
#define _INTEGRATORS_H

#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

template <class T>
class GSL_Integrator{
    const gsl_rng_type * t = gsl_rng_default;
    gsl_rng * R = gsl_rng_alloc(t);
    
    T* IntArgs;
    double *max,*min;

    const std::string Available[3] = {"Plain","Vegas","Miser"};

    public:

        GSL_Integrator(T& IntArgs0){
            IntArgs = &IntArgs0;
            max = new double[IntArgs->Dimension];
            for(int i =0;i<IntArgs->Dimension;i++)max[i]=1.0;
            min = new double[IntArgs->Dimension];
            gsl_rng_env_setup ();
        }

        ~GSL_Integrator(){
            delete [] max;
            delete [] min;
            gsl_rng_free(R);
        }

        void Evaluate(const std::string METHOD = "Plain"){
            double res,err;
            size_t Bound = IntArgs->NCalls/IntArgs->NIterations;
            gsl_monte_function F = {IntArgs->Function,IntArgs->Dimension,IntArgs->Parameters};

            if(METHOD=="Plain"){

                gsl_monte_plain_state *s = gsl_monte_plain_alloc(IntArgs->Dimension);
                for(int i=0;i<IntArgs->NIterations;i++){
                    gsl_monte_plain_integrate (&F,min,max,IntArgs->Dimension,Bound,R,s,&res,&err);
                    std::cout<<"Iteration["<<i<<"] "<<(i+1)*Bound<<" evaluations so far"<<std::endl;
                    std::cout<<res<<" +/- "<<err<<std::endl;
                }
                gsl_monte_plain_free(s);
            }

            else if(METHOD=="Vegas"){

                gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(IntArgs->Dimension);
                for(int i=0;i<IntArgs->NIterations;i++){
                    gsl_monte_vegas_integrate (&F,min,max,IntArgs->Dimension,Bound,R,s,&res,&err);
                    std::cout<<"Iteration["<<i<<"] "<<(i+1)*Bound<<" evaluations so far"<<std::endl;
                    std::cout<<res<<" +/- "<<err<<" ("<<gsl_monte_vegas_chisq(s)<<" chisq/dof)"<<std::endl;
                }
                gsl_monte_vegas_free(s);

            }
            else if(METHOD=="Miser"){

                gsl_monte_miser_state *s = gsl_monte_miser_alloc(IntArgs->Dimension);
                for(int i=0;i<IntArgs->NIterations;i++){
                    gsl_monte_miser_integrate (&F,min,max,IntArgs->Dimension,Bound,R,s,&res,&err);
                    std::cout<<"Iteration["<<i<<"] "<<(i+1)*Bound<<" evaluations so far"<<std::endl;
                    std::cout<<res<<" +/- "<<err<<std::endl;
                }
                gsl_monte_miser_free(s);

            }
            else{
                std::cout<<"Error: The integrator routine you requested "<<METHOD<<" is not available, the available routines are:"<<std::endl;
                for(auto a : Available) std::cout<<"    - "<<a<<std::endl;
                abort();
            }

        }     
};

#endif
