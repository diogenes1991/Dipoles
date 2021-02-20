#ifndef __LHA_PDF_H__
#define __LHA_PDF_H__

#include "PDF_Set.h"
#include "LHAPDF/LHAPDF.h"

class LHA_PDF : public PDF_Set{

    public:

        LHA_PDF(std::string NAME){
            LHAPDF::initPDFSetByName(NAME);
        }

        ~LHA_PDF(){};

        double Evaluate(int PID, double x, double MuFact){
            return LHAPDF::xfxphoton(x,MuFact,PID)/x;
        }

};

#endif
