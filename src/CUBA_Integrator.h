#ifndef __CUBA_INTEGRATOR_H__
#define __CUBA_INTEGRATOR_H__

#include "Montecarlo_Integrator.h"
#include <cuba.h>

class CUBA_Integrator : public Montecarlo_Integrator{

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

    public:
        
        CUBA_Integrator(integrand_t integrand, size_t dim){
            Integrand = integrand;
            Dimension = dim;
        }

        ~CUBA_Integrator(){};

        void Integrate(Specifications * mc_specs){
            
            std::string METHOD;
            if(mc_specs->Method=="")METHOD="Vegas";
            else METHOD = mc_specs->Method;

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
                std::cout<<"###################################################################"<<std::endl;
                std::cout<<"                                                                   "<<std::endl;
                std::cout<<"                  CUBA-Vegas integration initiated                 "<<std::endl;
                Vegas(Dimension,1,Integrand,Params,
                      NVec,RelErr,AbsErr,Flag,Seed,
                      MinEval,MaxEval,NStart,NIncrease,NBatch,
                      Grid, State, &Spin,
                      &NEval,&Fail,
                      Integral, Error, Prob);
                std::cout<<"                                                                   "<<std::endl;
                std::cout<<"###################################################################"<<std::endl;
            }
            
            

            else if(METHOD=="Suave"){
                std::cout<<"###################################################################"<<std::endl;
                std::cout<<"                                                                   "<<std::endl;
                std::cout<<"                  CUBA-Suave integration initiated                 "<<std::endl;
                Suave(Dimension,1,Integrand,Params, 
                      NVec,RelErr,AbsErr,Flag,Seed,
                      MinEval,MaxEval,NNew,NMin,
                      Flatness,State,&Spin,
                      &NRegions,&NEval,&Fail,
                      Integral,Error,Prob);
                std::cout<<"                                                                   "<<std::endl;
                std::cout<<"###################################################################"<<std::endl;
            }

            else{
                std::cout<<"Error: The integrator routine you requested "<<METHOD<<" is not available, the available CUBA routines are:"<<std::endl;
                for(auto a : Available) std::cout<<"    - "<<a<<std::endl;
                abort();
            }
        }
};

#endif
