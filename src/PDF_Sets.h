#ifndef _PDF_SETS_H
#define _PDF_SETS_H

#include "LHAPDF/LHAPDF.h"

class LHAPDF_Set{

    public:

        LHAPDF_Set(std::string NAME){
            LHAPDF::initPDFSetByName(NAME);
        }

        ~LHAPDF_Set(){};

        double Evaluate(int PID, double x, double MuFact){
            return LHAPDF::xfxphoton(x,MuFact,PID)/x;
        }

};

#endif
