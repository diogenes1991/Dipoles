#ifndef __DIPOLE_STRUCTURE_H_
#define __DIPOLE_STRUCTURE_H_

#include "nlox_process.h"
#include "Dipole_Definitions.h"
#include <string>
#include <unordered_map>
#include <iostream>

#define DEBUG 0

class DipoleStructure {

    public:

        Process* Proc;

        virtual ~DipoleStructure() {};

        virtual void SetECM(double sqrts) = 0;
        virtual void SetInMom(double* rand) = 0;
        virtual void SetFiMom(double* rand, double* J) = 0;

        virtual void Subtracted(std::string cp, double* rand, double* rval) = 0;
        virtual void PlusDistribution(std::string cp, double* rand, double mu, double* rval) = 0;
        virtual void Endpoint(std::string cp, double* rand, double mu, double* rval) = 0;

    };

#endif


