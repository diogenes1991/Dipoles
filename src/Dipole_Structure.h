#ifndef __DIPOLE_STRUCTURE_H_
#define __DIPOLE_STRUCTURE_H_

#include <vector>
#include <string>
#include "tred.h"
#include "poly3.h"
#include "four_vector.h"
#include "processconst.h"
#include "Dipole_Definitions.h"


class DipoleStructure {
    
    public:

        virtual ~DipoleStructure() {};

        virtual void Subtracted(double*,double*) = 0;
        virtual void PlusDistribution(double*,double*) = 0;
        virtual void Endpoint(double*,double*) = 0;
    
    };

#endif


