#ifndef __####SubProcHeader####_H_
#define __####SubProcHeader####_H_

#include "../PSP_Generator.h"
#include "../Dipole_Structure.h"
#include "../Utilities.h"

#define Next ####Next####

class ####SubProcName#### : public DipoleStructure{

    FourVector P;
    double Masses[Next];
    FourVector Momenta[Next];

    int nBorn;
    double** BornMasses;
    std::unordered_map<std::string,int> BornMap;
    FourVector BornMomenta[Next-1];

    public:
        
        ####SubProcName####(Process& process);
        ~####SubProcName####(){
            for (int i=0; i<nBorn; i++) delete BornMasses[i];
            delete [] BornMasses;
        }
        
        void SetECM(double sqrts);
        void SetInMom(double* rand);
        void SetFiMom(double* rand, double* J);

        void SetInMom(int BornNum);
        void SetFiMom(int BornNum, double* rand, double* J);

        void Subtracted(std::string cp, double* rand, double* rval);
        void PlusDistribution(std::string cp, double* rand, double mu, double* rval);
        void Endpoint(std::string cp, double* rand, double mu, double* rval);

    };

#endif