#ifndef __####SubProcHeader####_H_
#define __####SubProcHeader####_H_

#include "../PSP_Generator.h"
#include "../Dipole_Structure.h"

#define Next ####Next####

class ####SubProcName#### : public DipoleStructure{

    double sqrts;
    double masses[Next];
    FourVector momenta[Next];

    public:
        
        ####SubProcName####(Process& process);
        
        void SetECM(double sqrts);
        void SetInMom(double* rand);
        void SetFiMom(double* rand, double* J);

        void Subtracted(std::string cp, double* rand, double* rval);
        void PlusDistribution(std::string cp, double* rand, double mu, double* rval);
        void Endpoint(std::string cp, double* rand, double mu, double* rval);

    };

#endif