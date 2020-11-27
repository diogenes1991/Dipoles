#ifndef __MONTECARLOINTEGRATOR_H__
#define __MONTECARLOINTEGRATOR_H__

#include <iostream>

class Montecarlo_Integrator{

    public:

        struct Specifications{
            int RNGSeed = 1;
            void *Params = NULL;
            const int NVec = 1;
            size_t MaxEval = 1000000;
            size_t NStart = 1000;
            size_t NIncrease = 1000;
            double RelErr = 1E-4;
            double AbsErr = 1E-18;
        };

        virtual ~Montecarlo_Integrator(){};
        virtual void Integrate(Specifications * mc_specs, const std::string method) = 0;

};

#endif
