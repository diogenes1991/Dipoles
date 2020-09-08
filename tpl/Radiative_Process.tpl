#ifndef __RADIATIVE_PROCESS_H_
#define __RADIATIVE_PROCESS_H_

#include "Dipole_Definitions.h"

####Include Born####
####Include Radi####

class Process {
    
    Subprocess** RadiSubProc;
    Subprocess** BornSubProc;
    int nBorn,nRadi;
    std::map<std::string,int> BornMap;
    std::map<std::string,int> RadiMap;
    
    public:
        ProcessConst pc;
    
        Process() {
        nBorn = ####nBorn####;
        nRadi = ####nRadi####;
        BornSubProc = new Subprocess* [nBorn];
        RadiSubProc = new Subprocess* [nRadi];
        ####Construct Born####
        ####Construct Radi####
        }

        ~Process() {
        for (int i = 0; i != nBorn; i++) delete BornSubProc[i];
        for (int i = 0; i != nRadi; i++) delete RadiSubProc[i];
        delete [] BornSubProc;
        delete [] RadiSubProc;


        void call(std::string subprocess){

        }



  }

};

#endif


