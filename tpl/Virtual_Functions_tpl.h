#ifndef __####SubProcHeader####_H__
#define __####SubProcHeader####_H__

#include "PSP_Generator.h"
#include "Virtual_Structure.h"
#include "Utilities.h"

#define NextV ####Next####

class ####SubProcName#### : public VirtualStructure{

    FourVector P;
    double* BornMasses;
    int* BornPID;
    FourVector BornMomenta[NextV];

    public:
        
        ####SubProcName####(Process& process);
        ~####SubProcName####(){
            delete [] BornMasses;
            delete [] BornPID;
        }
        
        void SetECM(double sqrts);
        void SetInMom();
        void SetFiMom(double* rand, double* J);

        void GetMomenta(FourVector* p);
        void GetMasses(double* m);
        void GetPID(int* pid);

        void Born(std::string cp, double* rand, double* rval);
        void Virtual(std::string cp, double* rand, double mu, double* rval);

    };

#endif