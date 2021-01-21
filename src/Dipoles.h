#ifndef __DIPOLE_H__
#define __DIPOLE_H__

#include "Phase_Space_Tools.h"
#include "Constants.h"

class Dipole{

    double m_emit,m_spec;
    double mu_rn;
    double cp_sq;

    public:

        Dipole(double m_em, double m_sp){
            m_emit = m_em;
            m_spec = m_sp;
        }
        ~Dipole(){}

        static void setCoupling(double c) cp_sq = c;
        static void setCoupling(double c) mu_rn = c;

        virtual void BuildTildeMomenta(FourVector * RadiativeMomenta, FourVector * TildeMomenta) = 0;

        virtual void gFunction(FourVector p_em, FourVector p_sp, FourVector p_rd, double* rval) = 0;    
        virtual void CurlyGFunction(double mandelstam, double x, double Ix, double I1, double* rval) = 0;
        virtual void GFunction(double mandelstam, double* rval) = 0;
};

class EWK_Dipole : public Dipole{

};

class QCD_Dipole : public Dipole{
    
};

#endif
