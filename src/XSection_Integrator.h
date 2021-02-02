#ifndef __XSECTION_INTEGRATOR_H__
#define __XSECTION_INTEGRATOR_H__

#include "XSection.h"
#include "GSL_Integrator.h"
#include "CUBA_Integrator.h"

class XSection_Integrator{

    size_t NVarB = 3*NextV-4;
    size_t NVarP = 3*NextV-3;
    size_t NVarS = 3*NextR-4;

    Montecarlo_Integrator *BIntegrator, *SIntegrator, *PIntegrator;

    public:

        static XSection * XSec;

        struct XSection_Selector{
        std::string Channel;
        std::string Integrand;
        std::string Coupling;
        std::string Catalog;
        size_t NVars;
        };

        XSection_Integrator(std::string olp, std::string pdfname, std::string Integrator){
            
            std::cout<<"###################################################################"<<std::endl;
            std::cout<<"                                                                   "<<std::endl;
            std::cout<<"         XSection Integrator Instance Created                      "<<std::endl;
            std::cout<<"                                                                   "<<std::endl;
            std::cout<<"             1-Loop Provider     : "  <<   olp                      <<std::endl;
            std::cout<<"             PDF Set             : "  <<   pdfname                  <<std::endl;
            std::cout<<"             Integration Routine : "  <<   Integrator               <<std::endl;
            std::cout<<"                                                                   "<<std::endl;
            std::cout<<"###################################################################"<<std::endl;
            
            XSec = new XSection(olp,pdfname);

            if(pdfname!=""){
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
            delete XSec;
            delete BIntegrator;
            delete SIntegrator;
            delete PIntegrator;
        }
            
        static double GSL_Integrand(double *x, size_t dim, void* param){
            double rval;
            XSection_Selector xs;
            xs = *(XSection_Selector*)(param);
            XSec->SetXSection(xs.Catalog,xs.Integrand,xs.Channel,xs.Coupling,x,xs.NVars,&rval);
            return rval;
        }

        static int CUBA_Integrand(const int *ndim, const double x[], const int *ncomp, double f[], void* param){
            double rval;
            XSection_Selector xs;
            xs = *(XSection_Selector*)(param);
            size_t dim = *ndim;
            double y[dim];
            for(size_t i=0;i<dim;i++)y[i]=x[i];
            XSec->SetXSection(xs.Catalog,xs.Integrand,xs.Channel,xs.Coupling,y,xs.NVars,&rval);
            f[0] = rval;
            return 0;
        }

        void ComputeXSection(XSection_Selector XS, Montecarlo_Integrator::Specifications MC){
            Montecarlo_Integrator::Specifications mc;
            mc.Method = MC.Method;
            mc.NStart = MC.NStart;
            mc.NIncrease = MC.NIncrease;
            mc.MaxEval = MC.MaxEval;
            mc.RelErr = MC.RelErr;

            Montecarlo_Integrator * MCI;

            Clock C1();

            if (XS.Integrand=="Born"||XS.Integrand=="Virtual"){
                MCI = BIntegrator;
                XS.Catalog = "Virtuals";
                XS.NVars = MCI->Dimension;
            }
            else if (XS.Integrand=="Endpoint"){
                MCI = BIntegrator;
                XS.Catalog = "Reals";
                XS.NVars = MCI->Dimension;
            }
            else if (XS.Integrand=="PlusDistribution"){
                MCI = PIntegrator;
                XS.Catalog = "Reals";
                XS.NVars = MCI->Dimension;
            }
            else if (XS.Integrand=="Subtracted"){
                MCI = SIntegrator;
                XS.Catalog = "Reals";
                XS.NVars = MCI->Dimension;
            }
            else{
                std::cout<<"Error: Unavailable integrand "<<XS.Integrand<<" for "<<XS.Channel<<" channel"<<std::endl;
                abort();
            }

            mc.Params = &XS;
            MCI->Integrate(&mc);
            std::cout<<"The integration took: ";
            C1.ShowTime(); 
        }

};

XSection * XSection_Integrator::XSec = NULL;

#endif
