#ifndef __INTEGRANDS_H_
#define __INTEGRANDS_H_

#include "Radiative_Process.h"

class Integrands{
    
    Process* Proc = NULL;
    std::unordered_map<std::string,std::vector<std::string>> AmpMap;

    public:
        
        Integrands Integrands(){

            Proc = new Process();
            // ####AmpMap####

        }

        Integrands ~Integrands(){
            delete Proc;
        }

        // Now that all the subprocess have been initialized we can 
        // build the integrands. These are the naming convetions:
        //     - The SubProcesses are called by born tag 
        //     - The CP is the radiative CP
        //     - We will build the dependency tree here at construction
        //       it will have an array of names of the radiative processes 
        //       that need to be included. 

        void Subtracted(std::string SP, std::string CP, std::vector<FourVector> P, double* rval){

            
        }

  }

};

#endif


