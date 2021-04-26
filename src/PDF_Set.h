#ifndef __PDF_SET_H__
#define __PDF_SET_H__

class PDF_Set{

    public:

        virtual ~PDF_Set(){};
        virtual double Evaluate(int PID, double x, double MuFact) = 0;
        virtual double Alpha_S(double Q) = 0;
};

#endif
