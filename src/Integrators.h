#ifndef _INTEGRATORS_H
#define _INTEGRATORS_H

#include <iostream>

#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

#include <cuba.h>

struct montecarlo_specs{
    int RNGSeed = 1;
    void *Params = NULL;
    const int NVec = 1;
    size_t MaxEval = 1000000;
    size_t NStart = 1000;
    size_t NIncrease = 1000;
    double RelErr = 1E-4;
    double AbsErr = 1E-18;
};

class GSL_Integrator{
    
    const gsl_rng_type * t = gsl_rng_default;
    gsl_rng * R = gsl_rng_alloc(t);
    
    double *max,*min;
    const std::string Available[3] = {"Plain","Vegas","Miser"};
    
    size_t Dimension;
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

        void Integrate(montecarlo_specs * mc_specs, const std::string METHOD = "Plain"){
            
            double res,err;
            gsl_monte_function F = {Integrand,Dimension,mc_specs->Params};

            if(METHOD=="Plain"){
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
                    std::cout<<"Iteration["<<Iter<<"] "<<Calls<<" evaluations so far"<<std::endl<<std::endl;
                    std::cout<<res<<" +/- "<<err<<std::endl;
                    if(std::abs(err/res)<mc_specs->RelErr) break;
                }
                gsl_monte_plain_free(s);
            }
            
            else if(METHOD=="Vegas"){
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
            }
            
            else if(METHOD=="Miser"){
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
                    std::cout<<"Iteration["<<Iter<<"] "<<Calls<<" evaluations so far"<<std::endl<<std::endl;
                    std::cout<<res<<" +/- "<<err<<std::endl;
                    if(std::abs(err/res)<mc_specs->RelErr) break;
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

class CUBA_Integrator{

    // void Vegas(const int ndim, const int ncomp,
    // integrand_t integrand, void *userdata, const int nvec,
    // const double epsrel, const double epsabs,
    // const int flags, const int seed,
    // const int mineval, const int maxeval,
    // const int nstart, const int nincrease, const int nbatch,
    // const int gridno, const char *statefile, void *spin,
    // int *neval, int *fail,
    // double integral[], double error[], double prob[])

    // void Suave(const int ndim, const int ncomp,
    // integrand_t integrand, void *userdata, const int nvec,
    // const double epsrel, const double epsabs,
    // const int flags, const int seed,
    // const int mineval, const int maxeval,
    // const int nnew, const int nmin,
    // const double flatness, const char *statefile, void *spin,
    // int *nregions, int *neval, int *fail,
    // double integral[], double error[], double prob[]);
    
    // void Divonne(const int ndim, const int ncomp,
    // integrand_t integrand, void *userdata, const int nvec,
    // const double epsrel, const double epsabs,
    // const int flags, const int seed,
    // const int mineval, const int maxeval,
    // const int key1, const int key2, const int key3,
    // const int maxpass, const double border,
    // const double maxchisq, const double mindeviation,
    // const int ngiven, const int ldxgiven, double xgiven[],
    // const int nextra, peakfinder_t peakfinder,
    // const char *statefile, void *spin,
    // int *nregions, int *neval, int *fail,
    // double integral[], double error[], double prob[]);
    
    // void Cuhre(const int ndim, const int ncomp,
    // integrand_t integrand, void *userdata, const int nvec,
    // const double epsrel, const double epsabs,
    // const int flags,
    // const int mineval, const int maxeval,
    // const int key, const char *statefile, void *spin,
    // int *nregions, int *neval, int *fail,
    // double integral[], double error[], double prob[]);

    const std::string Available[4]={"Vegas","Suave","Divonne","Cuhre"};
    integrand_t Integrand;
    size_t Dimension;

    public:
        CUBA_Integrator(integrand_t integrand, size_t dim){
            Integrand = integrand;
            Dimension = dim;
        }

        ~CUBA_Integrator(){};

        void Integrate(montecarlo_specs * mc_specs, std::string METHOD = "Vegas"){
            
            // Common Arguments 

            void* Params = mc_specs->Params;
            const int NVec = mc_specs->NVec;

            const int NStart = mc_specs->NStart;
            const int NIncrease = mc_specs->NIncrease;
            const int MinEval = NStart;
            const int MaxEval = mc_specs->MaxEval;

            const double RelErr = mc_specs->RelErr;
            const double AbsErr = mc_specs->AbsErr;

            const int Flag = 1;
            const int Seed = mc_specs->RNGSeed;

            const int Grid = 0;
            const char* State = NULL;
            void* Spin = NULL;

            int NEval,Fail;
            double Integral[NVec],Error[NVec],Prob[NVec];

            // Vegas Specific

            const int NBatch = 1000;

            // Suave Specific

            const int NNew = mc_specs->NIncrease; 
            const int NMin = 50;
            const double Flatness = 2;
            int NRegions = 1;


            if(METHOD=="Vegas"){
                Vegas(Dimension,1,Integrand,Params,
                      NVec,RelErr,AbsErr,Flag,Seed,
                      MinEval,MaxEval,NStart,NIncrease,NBatch,
                      Grid, State, &Spin,
                      &NEval,&Fail,
                      Integral, Error, Prob);
            }

            else if(METHOD=="Suave"){
                Suave(Dimension,1,Integrand,Params, 
                      NVec,RelErr,AbsErr,Flag,Seed,
                      MinEval,MaxEval,NNew,NMin,
                      Flatness,State,&Spin,
                      &NRegions,&NEval,&Fail,
                      Integral,Error,Prob);
            }

            else{
                std::cout<<"Error: The integrator routine you requested "<<METHOD<<" is not available, the available routines are:"<<std::endl;
                for(auto a : Available) std::cout<<"    - "<<a<<std::endl;
                abort();
            }
        }
};

#endif
