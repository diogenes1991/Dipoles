#ifndef __DIPOLE_STRUCTURE_H__
#define __DIPOLE_STRUCTURE_H__

#include "nlox_process.h"
#include "Dipole_Definitions.h"
#include <string>
#include <unordered_map>
#include <iostream>

#define DEBUG 0

class DipoleStructure {

    public:

        Process * Proc;
        std::unordered_map<std::string,int> BornMap;

        FourVector P;

        FourVector* Momenta;
        double* Masses;
        int* PID;

        int nBorn,nPar;
        FourVector* BornMomenta;
        double** BornMasses;
        int** BornPID;
        
        virtual ~DipoleStructure(){
            for (int i=0; i<nBorn; i++){
                delete BornMasses[i];
                delete BornPID[i];
            }
            delete [] BornMasses;
            delete [] BornPID;

            delete BornMomenta;
            delete PID;
            delete Masses;
            delete Momenta;
        }

        void SetECM(double sqrts){
            P.p0=sqrts;
            double rIn[2] = {0.0,0.0};
            SetInMom(rIn);
        }

        void SetFiMom(double* rFi, double* J){
            FourVector PFi[nPar-2];
            double mFi[nPar-2];
            double dummy = 1.0;
            for (int i=0;i<nPar-2;i++) mFi[i] = Masses[i+2];
            Recursive_PSP(P,nPar-2,PFi,mFi,rFi,dummy);
            for (int i=0;i<nPar-2;i++) Momenta[i+2] = PFi[i];
            *J = dummy;    
        }

        void SetInMom(double* rIn){
            FourVector PIn[2];
            double mIn[2]={Masses[0],Masses[1]};
            double dummy = 1.0;
            Recursive_PSP(P,2,PIn,mIn,rIn,dummy);
            Momenta[0] = PIn[0];
            Momenta[1] = PIn[1];
        }

        void SetFiMom(int BornNum, double* rFi, double* J){
            FourVector PFi[nPar-3];
            double mFi[nPar-3];
            double dummy = 1.0;
            for (int i=0;i<nPar-3;i++) mFi[i] = BornMasses[BornNum][i+2];
            Recursive_PSP(P,nPar-1,PFi,mFi,rFi,dummy);
            for (int i=0;i<nPar-3;i++) BornMomenta[i+2] = PFi[i];
            *J = dummy;    
        }

        void SetInMom(int BornNum){
            FourVector PIn[2];
            double mIn[2]={BornMasses[BornNum][0],BornMasses[BornNum][1]};
            double rIn[2]={0.0,0.0};
            double dummy = 1.0;
            Recursive_PSP(P,2,PIn,mIn,rIn,dummy);
            BornMomenta[0] = PIn[0];
            BornMomenta[1] = PIn[1];
        }

        void GetMomenta(FourVector* p){

            for(int i=0;i<nPar;i++)p[i]=Momenta[i];
        }

        void GetMasses(double* m){

            for(int i=0;i<nPar;i++)m[i]=Masses[i];
        }

        void GetPID(int* pid){

            for(int i=0;i<nPar;i++)pid[i]=PID[i];
        }

        virtual void Subtracted(std::string cp, double* rand, double mu, double* rval) = 0;
        virtual void PlusDistribution(std::string cp, double* rand, double mu, double* rval) = 0;
        virtual void Endpoint(std::string cp, double* rand, double mu, double* rval) = 0;

};

#endif


