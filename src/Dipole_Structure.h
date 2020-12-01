
#ifndef __DIPOLE_STRUCTURE_H__
#define __DIPOLE_STRUCTURE_H__

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

        virtual void GetMomenta(FourVector* p) = 0;
        virtual void GetMasses(double* m) = 0;
        virtual void GetPID(int* pid) = 0;

        virtual void SetInMom(int BornNum) = 0;
        virtual void SetFiMom(int BornNum, double* rand, double* J) = 0;

        virtual void Subtracted(std::string cp, double* rand, double mu, double* rval) = 0;
        virtual void PlusDistribution(std::string cp, double* rand, double mu, double* rval) = 0;
        virtual void Endpoint(std::string cp, double* rand, double mu, double* rval) = 0;

    };

#endif


