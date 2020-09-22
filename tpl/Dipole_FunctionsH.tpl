#ifndef __####SubProcHeader####_H_
#define __####SubProcHeader####_H_

#include "Dipole_Structure.h"

class ####SubProcName#### : public DipoleStructure{

    public:
        ####SubProcName####();

        void Subtracted(std::string cp, double* p, double mu, double* rval, double* acc);
        void PlusDistribution(std::string cp, double* p, double mu, double* rval, double* acc);
        void Endpoint(std::string cp, double* p, double mu, double* rval, double* acc);

    };

#endif