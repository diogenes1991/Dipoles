#ifndef __DIPOLE_STRUCTURE_H_
#define __DIPOLE_STRUCTURE_H_

#include "nlox_process.h"
#include "Dipole_Definitions.h"

class DipoleStructure {

    Process* Proc;
    
    public:

        virtual ~DipoleStructure() {};

        virtual void Subtracted(std::string cp, std::vector<FourVector> p, double mu, double* rval, double* acc) = 0;
        virtual void PlusDistribution(double*,double*) = 0;
        virtual void Endpoint(double*,double*) = 0;
    
    };

#endif


