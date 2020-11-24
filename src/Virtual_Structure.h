#ifndef __VIRTUAL_STRUCTURE_H_
#define __VIRTUAL_STRUCTURE_H_

#include "nlox_process.h"
#include <string>
#include <unordered_map>


class VirtualStructure {

    public:

        Process* Proc;

        virtual ~VirtualStructure() {};

        virtual void SetECM(double sqrts) = 0;
        virtual void SetInMom() = 0;
        virtual void SetFiMom(double* rand, double* J) = 0;

        virtual void Born(std::string cp, double* rand, double* rval) = 0;
        virtual void Virtual(std::string cp, double* rand, double mu, double* rval) = 0;

    };

#endif


