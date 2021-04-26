#ifndef __DUMMY_PDF_H__
#define __DUMMY_PDF_H__

#include "PDF_Set.h"

class Dummy_PDF : public PDF_Set{

    public:

        Dummy_PDF(){}

        ~Dummy_PDF(){}

        double Evaluate(int PID, double x, double MuFact) {return 1.0;}
        double Alpha_S(double Q) {return 0.118;}

};

#endif
