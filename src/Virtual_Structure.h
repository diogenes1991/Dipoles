#ifndef __VIRTUAL_STRUCTURE_H__
#define __VIRTUAL_STRUCTURE_H__

#include "OLP.h"
#include "Model.h"
#include <string>
#include <unordered_map>

class VirtualStructure {

    public:

        OLP * Provider;
        Model * model;
        int nPar;
        double* BornMasses;
        int* BornPID;
        FVector * BornMomenta;

        virtual ~VirtualStructure(){
            delete [] BornMomenta;
            delete [] BornMasses;
            delete [] BornPID;
        }
        
        void BGenerate(double sqrts, double* rand, double* Jac){
            
            FVector P(sqrts,0,0,0);
            double dummy;

            FVector PIn[2];
            double mIn[2]={BornMasses[0],BornMasses[1]};
            double rIn[2]={0.0,0.0};
            dummy = 1.0;
            Recursive_PSP(P,2,PIn,mIn,rIn,dummy);
            BornMomenta[0] = PIn[0];
            BornMomenta[1] = PIn[1];
            
            FVector PFi[nPar-2];
            double mFi[nPar-2];
            for (int i=0;i<nPar-2;i++) mFi[i] = BornMasses[i+2];
            dummy = 1.0;
            Recursive_PSP(P,nPar-2,PFi,mFi,rand,dummy);
            for (int i=0;i<nPar-2;i++) BornMomenta[i+2] = PFi[i];
            *Jac = dummy;   

        }

        void GetMomenta(FVector* p){

            for(int i=0;i<nPar;i++)p[i]=BornMomenta[i];
        }

        void GetMasses(double* m){

            for(int i=0;i<nPar;i++)m[i]=BornMasses[i];
        }

        void GetPID(int* pid){

            for(int i=0;i<nPar;i++)pid[i]=BornPID[i];
        }

        virtual void Born(std::string cp, double sqrts,  double* rand, double mu, double* rval) = 0;
        virtual void Virtual(std::string cp, double sqrts,  double* rand, double mu, double* rval) = 0;

    };

#endif


