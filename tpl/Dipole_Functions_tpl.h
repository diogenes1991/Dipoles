#ifndef __####SubProcHeader####_H__
#define __####SubProcHeader####_H__

#include "PSP_Generator.h"
#include "Dipole_Structure.h"
#include "Utilities.h"

#define NextR ####Next####

class ####SubProcName#### : public DipoleStructure{

    FourVector P;
    double Masses[NextR];
    FourVector Momenta[NextR];
    int PID[NextR];

    int nBorn;
    double** BornMasses;
    int** BornPID;
    std::unordered_map<std::string,int> BornMap;
    FourVector BornMomenta[NextR-1];

    public:
        
        ####SubProcName####(Process& process);
        ~####SubProcName####(){
            for (int i=0; i<nBorn; i++){
                delete BornMasses[i];
                delete BornPID[i];
            }
            delete [] BornMasses;
            delete [] BornPID;
        }
        
        void SetECM(double sqrts);
        void SetInMom(double* rand);
        void SetFiMom(double* rand, double* J);

        void GetMomenta(FourVector* p);
        void GetMasses(double* m);
        void GetPID(int* pid);

        void SetInMom(int BornNum);
        void SetFiMom(int BornNum, double* rand, double* J);

        void Subtracted(std::string cp, double* rand, double mu, double* rval);
        void PlusDistribution(std::string cp, double* rand, double mu, double* rval);
        void Endpoint(std::string cp, double* rand, double mu, double* rval);

    };

#endif