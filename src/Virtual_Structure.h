#ifndef __VIRTUAL_STRUCTURE_H__
#define __VIRTUAL_STRUCTURE_H__

#include "nlox_process.h"
#include <string>
#include <unordered_map>


class VirtualStructure {

    public:

        Process* Proc;
        FourVector P;
        int nPar;
        double* BornMasses;
        int* BornPID;
        FourVector * BornMomenta;

        virtual ~VirtualStructure(){
            delete [] BornMomenta;
            delete [] BornMasses;
            delete [] BornPID;
        }
        
        void SetECM(double ECM){
            P.p0=ECM;
            SetInMom();
        }

        void SetFiMom(double* rFi, double* J){
            FourVector PFi[nPar-2];
            double mFi[nPar-2];
            double dummy = 1.0;
            for (int i=0;i<nPar-2;i++) mFi[i] = BornMasses[i+2];
            Recursive_PSP(P,nPar-2,PFi,mFi,rFi,dummy);
            for (int i=0;i<nPar-2;i++) BornMomenta[i+2] = PFi[i];
            *J = dummy;    
        }

        void SetInMom(){
            FourVector PIn[2];
            double mIn[2]={BornMasses[0],BornMasses[1]};
            double rIn[2]={0.0,0.0};
            double dummy = 1.0;
            Recursive_PSP(P,2,PIn,mIn,rIn,dummy);
            BornMomenta[0] = PIn[0];
            BornMomenta[1] = PIn[1];
        }

        void GetMomenta(FourVector* p){

            for(int i=0;i<nPar;i++)p[i]=BornMomenta[i];
        }

        void GetMasses(double* m){

            for(int i=0;i<nPar;i++)m[i]=BornMasses[i];
        }

        void GetPID(int* pid){

            for(int i=0;i<nPar;i++)pid[i]=BornPID[i];
        }

        virtual void Born(std::string cp, double* rand, double* rval) = 0;
        virtual void Virtual(std::string cp, double* rand, double mu, double* rval) = 0;

    };

#endif


