#ifndef __XSECTION_H__
#define __XSECTION_H__

#include "Real.h"
#include "Virtual.h"
#include "Integrators.h"
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
            
            int PID[4];
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

XSection * X1;
void XSectionStart(std::string pdfname);
void XSectionEnd();

class XSecCalc{

    size_t NVarB = 3*NextV-4;
    size_t NVarP = 3*NextV-3;
    size_t NVarS = 3*NextR-4;

    bool usinggsl = false;
    GSL_Integrator *GBIntegrator,*GSIntegrator,*GPIntegrator;
    
    bool usingcuba = false;
    CUBA_Integrator *CBIntegrator,*CSIntegrator,*CPIntegrator;

    public:

        struct xsec_sel{
        std::string Channel;
        std::string Integrand;
        std::string Coupling;
        };

        XSecCalc(std::string Integrator = "All"){
            
            if(X1->ispdfset()){
                NVarB+=2;
                NVarS+=2;
                NVarP+=2;
            }

            if(Integrator=="GSL"||Integrator=="All"){
                usinggsl = true;
                GBIntegrator = new GSL_Integrator(GSL_Integrand,NVarB);
                GSIntegrator = new GSL_Integrator(GSL_Integrand,NVarS);
                GPIntegrator = new GSL_Integrator(GSL_Integrand,NVarP);
            }
            if(Integrator=="CUBA"||Integrator=="All"){
                usingcuba = true;
                CBIntegrator = new CUBA_Integrator(CUBA_Integrand,NVarB);
                CSIntegrator = new CUBA_Integrator(CUBA_Integrand,NVarS);
                CPIntegrator = new CUBA_Integrator(CUBA_Integrand,NVarP);
            }
        }

        ~XSecCalc(){
            if(usinggsl){
                delete GBIntegrator;
                delete GSIntegrator;
                delete GPIntegrator;
            }
            if(usingcuba){
                delete CBIntegrator;
                delete CSIntegrator;
                delete CPIntegrator;
            }
        }
            
        static double GSL_Integrand(double *x, size_t dim, void* param){
            double rval;
            xsec_sel xs;
            xs = *(xsec_sel*)(param);
            X1->SetXSection(xs.Integrand,xs.Channel,xs.Coupling,x,&rval);
            return rval;
        }

        static int CUBA_Integrand(const int *ndim, const double x[], const int *ncomp, double f[], void* param){
            double rval;
            xsec_sel xs;
            xs = *(xsec_sel*)(param);
            size_t dim = *ndim;
            double y[dim];
            for(size_t i=0;i<dim;i++)y[i]=x[i];
            X1->SetXSection(xs.Integrand,xs.Channel,xs.Coupling,y,&rval);
            f[0] = rval;
            return 0;
        }

        void ComputeXSections(std::string Integrand, std::string Channel, std::string Coupling, std::string Integrator, std::string Method){
            montecarlo_specs mc;
            xsec_sel xs;
            xs.Channel = Channel;
            xs.Coupling = Coupling;
            xs.Integrand = Integrand;
            mc.NStart = 5000;
            mc.NIncrease = 750;
            mc.MaxEval = 100000;
            mc.RelErr = 1E-10;
            mc.Params = &xs;
            if(Integrator=="CUBA")CBIntegrator->Integrate(&mc,Method);
            else if(Integrator=="GSL")GBIntegrator->Integrate(&mc,Method);
            else{
                std::cout<<"Error: Unimplemented integrator routine: "<<Integrator<<std::endl;
                abort();
            }
        }
};

#endif
