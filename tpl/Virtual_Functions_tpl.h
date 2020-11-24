#ifndef __####SubProcHeader####_H_
#define __####SubProcHeader####_H_

#include "PSP_Generator.h"
#include "Virtual_Structure.h"
#include "Utilities.h"

#define NextV ####Next####

class ####SubProcName#### : public VirtualStructure{

    FourVector P;
    double* BornMasses;
    FourVector BornMomenta[NextV];

    public:
        
        ####SubProcName####(Process& process);
        ~####SubProcName####(){
            delete [] BornMasses;
        }
        
        void SetECM(double sqrts);
        void SetInMom();
        void SetFiMom(double* rand, double* J);

        void Born(std::string cp, double* rand, double* rval);
        void Virtual(std::string cp, double* rand, double mu, double* rval);

    };

#endif