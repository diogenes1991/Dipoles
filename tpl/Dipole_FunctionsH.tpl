#ifndef __####SubProcHeader####_H_
#define __####SubProcHeader####_H_

#include "Dipole_Structure.h"

class ####SubProcName#### : pubic DipoleStructure{

    public:
        ####SubProcName#### :: ####SubProcName####();

        void Subtracted(std::string cp, std::vector<FourVector> p, double mu, double* rval, double* acc);
        void PlusDistribution(double* p, double* rval);
        void Endpoint(double* p, double* rval);

    
    };

#endif


