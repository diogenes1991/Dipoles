#ifndef __####SubProcHeader####_H_
#define __####SubProcHeader####_H_

#include "DipoleStructure.h"

####SubProcMat####

class ####SubProcName#### : pubic DipoleStructure{

    public:
        ####SubProcName#### :: ####SubProcName####();

        void Subtracted(double* p, double* rval);
        void PlusDistribution(double* p, double* rval);
        void Endpoint(double* p, double* rval);

    
    };

#endif


