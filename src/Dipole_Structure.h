#ifndef __DIPOLE_STRUCTURE_H_
#define __DIPOLE_STRUCTURE_H_

#include "nlox_process.h"
#include "Dipole_Definitions.h"
#include <string>
#include <unordered_map>
#include <iostream>

class DipoleStructure {

    public:

        Process* Proc;

        virtual ~DipoleStructure() {};

        virtual void Subtracted(std::string cp, std::vector<FourVector> p, double* rval) = 0;
        virtual void PlusDistribution(std::string cp, std::vector<FourVector> p, double mu, double* rval, double* acc) = 0;
        virtual void Endpoint(std::string cp, std::vector<FourVector> p, double mu, double* rval) = 0;

    };

#endif


