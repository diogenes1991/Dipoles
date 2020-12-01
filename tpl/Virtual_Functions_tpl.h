#ifndef __####SubProcHeader####_H__
#define __####SubProcHeader####_H__

#include "PSP_Generator.h"
#include "Virtual_Structure.h"
#include "Utilities.h"

#define NextV ####Next####

class ####SubProcName#### : public VirtualStructure{

    public:
        
        ####SubProcName####(Process& process);
        
        void Born(std::string cp, double* rand, double* rval);
        void Virtual(std::string cp, double* rand, double mu, double* rval);

    };

#endif