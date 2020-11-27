#ifndef __PDF_SETS_H__
#define __PDF_SETS_H__

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
