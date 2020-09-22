#ifndef __DIPOLE_STRUCTURE_H_
#define __DIPOLE_STRUCTURE_H_

#include "nlox_process.h"
#include <string>
#include <unordered_map>
#include <iostream>
#include "Dipole_Definitions.h"

class DipoleStructure {

    public:

        Process* Proc;

        virtual ~DipoleStructure() {};

        virtual void Subtracted(std::string cp, double* p, double mu, double* rval, double* acc) = 0;
        virtual void PlusDistribution(std::string cp, double* p, double mu, double* rval, double* acc) = 0;
        virtual void Endpoint(std::string cp, double* p, double mu, double* rval, double* acc) = 0;

    };

#endif


