#ifndef __####SubProcHeader####_H__
#define __####SubProcHeader####_H__

#include "PSP_Generator.h"
#include "Dipole_Structure.h"
#include "Utilities.h"

#define NextR ####Next####

class ####SubProcName#### : public DipoleStructure{

    public:
        
        ####SubProcName####(OLP& Prov, Model& Mod);
        
        void Subtracted(std::string cp, double sqrts, double* rand, double mu, double* rval);
        void PlusDistribution(std::string cp, double sqrts, double* rand, double mu, double* rval);
        void Endpoint(std::string cp, double sqrts, double* rand, double mu, double* rval);

};

#endif