#ifndef __XSECTION_INTEGRATOR_H__
#define __XSECTION_INTEGRATOR_H__

#include "XSection.h"
#include "GSL_Integrator.h"
#include "CUBA_Integrator.h"

XSection * X1;
void XSectionStart(std::string pdfname);
void XSectionEnd();

class XSection_Integrator{

    size_t NVarB = 3*NextV-4;
    size_t NVarP = 3*NextV-3;
    size_t NVarS = 3*NextR-4;

    Montecarlo_Integrator *BIntegrator, *SIntegrator, *PIntegrator;

    public:

        struct XSection_Selector{
        std::string Channel;
        std::string Integrand;
        std::string Coupling;
        };

        XSection_Integrator(std::string Integrator){
            
            if(X1->ispdfset()){
                NVarB+=2;
                NVarS+=2;
                NVarP+=2;
            }

            if(Integrator=="GSL"){
                BIntegrator = new GSL_Integrator(GSL_Integrand,NVarB);
                SIntegrator = new GSL_Integrator(GSL_Integrand,NVarS);
                PIntegrator = new GSL_Integrator(GSL_Integrand,NVarP);
            }
            else if(Integrator=="CUBA"){
                BIntegrator = new CUBA_Integrator(CUBA_Integrand,NVarB);
                SIntegrator = new CUBA_Integrator(CUBA_Integrand,NVarS);
                PIntegrator = new CUBA_Integrator(CUBA_Integrand,NVarP);
            }
            else{
                std::cout<<"Error: Integrator interface not supported:"<<Integrator<<std::endl;
                abort();
            }
        }

        ~XSection_Integrator(){
            delete BIntegrator;
            delete SIntegrator;
            delete PIntegrator;
        }
            
        static double GSL_Integrand(double *x, size_t dim, void* param){
            double rval;
            XSection_Selector xs;
            xs = *(XSection_Selector*)(param);
            X1->SetXSection(xs.Integrand,xs.Channel,xs.Coupling,x,&rval);
            return rval;
        }

        static int CUBA_Integrand(const int *ndim, const double x[], const int *ncomp, double f[], void* param){
            double rval;
            XSection_Selector xs;
            xs = *(XSection_Selector*)(param);
            size_t dim = *ndim;
            double y[dim];
            for(size_t i=0;i<dim;i++)y[i]=x[i];
            X1->SetXSection(xs.Integrand,xs.Channel,xs.Coupling,y,&rval);
            f[0] = rval;
            return 0;
        }

        void ComputeXSection(std::string Integrand, std::string Channel, std::string Coupling, std::string Integrator, std::string Method){
            Montecarlo_Integrator::Specifications mc;
            XSection_Selector xs;
            xs.Channel = Channel;
            xs.Coupling = Coupling;
            xs.Integrand = Integrand;
            mc.NStart = 5000;
            mc.NIncrease = 750;
            mc.MaxEval = 100000;
            mc.RelErr = 1E-10;
            mc.Params = &xs;
            BIntegrator->Integrate(&mc,Method);
        }
};

#endif
