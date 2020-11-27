#ifndef __XSECTION_H__
#define __XSECTION_H__

#include "Real.h"
#include "Virtual.h"
#include "PDF_Sets.h"
#include "Analysis.h"
#include "Constants.h"

class XSection{

    bool usingpdfs = false;
    bool usingjetalg = false;

    double sqrts;
    double mu_ren;
    double mu_fac;

    Process * Proc;
    RealIntegrands * Reals;
    VirtualIntegrands * Virtuals;
    LHAPDF_Set * PDF;

    public:

        XSection(std::string pdfset = "none"){
            if(pdfset!="none"){
                PDF = new LHAPDF_Set(pdfset);
                usingpdfs = true;
            }
            Proc = new Process();
            Reals  = new RealIntegrands(Proc);
            Virtuals = new VirtualIntegrands(Proc);
        }

        ~XSection(){
            if(usingpdfs) delete PDF;
            delete Virtuals;
            delete Reals;
            delete Proc;
        }

        bool ispdfset(){return usingpdfs;}

        void SetScales(double sqrts0, double mu_ren0, double mu_fac0){
            sqrts = sqrts0;
            mu_ren = mu_ren0;
            mu_fac = mu_fac0;
            Virtuals->setECM(sqrts);
            Reals->setECM(sqrts);
        }

        void SetXSection(std::string Integrand, std::string Channel, std::string Coupling, double* x, double* xsec){
            
            double sqrtshat;
            double prefactor=1.0;
            int nfin = NextV-2;
            
            int PID[NextV];
            Virtuals->GetPID(Channel,PID);

            double mass[NextV];
            Virtuals->GetMasses(Channel,mass);
            
            double sqrtshat_min = 0;
            
            if(usingpdfs){
                for(int i=2;i<NextV;i++)sqrtshat_min+=mass[i];
                double x1_min = mass[0]/sqrts;
                double x2_min = mass[1]/sqrts;
                double x1_max = 0.5 * (1.0+sqrt(1.0-4.0*mass[0]*mass[0]/sqrts/sqrts));
                double x2_max = 0.5 * (1.0+sqrt(1.0-4.0*mass[1]*mass[1]/sqrts/sqrts));
                double x1 = x1_min + x[3*nfin-4] * (x1_max-x1_min);
                double x2 = x2_min + x[3*nfin-3] * (x2_max-x2_min);
                sqrtshat = sqrt(mass[0]*mass[0]+mass[1]*mass[1]+mass[0]*mass[0]*mass[1]*mass[1]/(sqrts*sqrts*x1*x2)+sqrts*sqrts*x1*x2);
                prefactor = PDF->Evaluate(PID[0],x1,mu_fac)*PDF->Evaluate(PID[1],x2,mu_fac)+PDF->Evaluate(PID[1],x1,mu_fac)*PDF->Evaluate(PID[0],x2,mu_fac);
                }
            else{
                sqrtshat = sqrts;
            }

            if(sqrtshat > sqrtshat_min){
                
                prefactor *= 1.0/(2.0*sqrt(lambda(sqrtshat*sqrtshat,mass[0]*mass[0],mass[1]*mass[1])));
                
                double partxsec;
                FourVector p[NextV];
                Virtuals->setECM(sqrtshat);
                if(Integrand=="Born")Virtuals->Born(Channel,Coupling,x,&partxsec);
                else if (Integrand=="Virtual")Virtuals->Virtual(Channel,Coupling,x,mu_ren,&partxsec);
                else partxsec = 0;

                double reweight;
                Virtuals->GetMomenta(Channel,p);
                Analysis::ReweightEvent(p,mass,PID,NextV,&reweight);
                prefactor *= reweight;

                *xsec = prefactor*partxsec;
            }
            else *xsec = 0;

            *xsec *= GeVtoPB;
        }
};

#endif
