#ifndef __####SubProcHeader####_H_
#define __####SubProcHeader####_H_

#include "../PSP_Generator.h"
#include "../Dipole_Structure.h"

class ####SubProcName#### : public DipoleStructure{

    int Next;
    std::vector<double> masses;
    FourVector PCM;
    FourVector PIn[2];

    public:
        
        ####SubProcName####(Process& process);
        
        void setECM(double sqrts);

        void Subtracted(std::string cp, double* rand, double* rval);
        void PlusDistribution(std::string cp, double* rand, double mu, double* rval);
        void Endpoint(std::string cp, double* rand, double mu, double* rval);

    };

#endif