#ifndef __####SubProcHeader####_H__
#define __####SubProcHeader####_H__

#include "PSP_Generator.h"
#include "Virtual_Structure.h"
#include "Utilities.h"

#define NextV ####Next####

class ####SubProcName#### : public VirtualStructure{

    public:
        
        ####SubProcName####(OLP& Prov, Model& Mod);
        
        void Born(std::string cp, double sqrts, double* rand, double mu, double* rval);
        void Virtual(std::string cp, double sqrts, double* rand, double mu, double* rval);

    };

#endif