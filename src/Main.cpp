#include "Input.h"
#include "XSection_Integrator.h"

class Madisqe : public Input{

    XSection_Integrator * XSec_Int;

    Montecarlo_Integrator::Specifications MC_SP;
    XSection_Integrator::XSection_Selector XS;

    public:

        Madisqe(std::string InputFileName){
            
            LoadInput(InputFileName,InputFile);
            
            std::cout<<"#############################################################"<<std::endl;
            std::cout<<"                                                             "<<std::endl;
            std::cout<<"   ███╗░░░███╗░█████╗░██████╗░██╗░██████╗░██████╗░███████╗   "<<std::endl;
            std::cout<<"   ████╗░████║██╔══██╗██╔══██╗██║██╔════╝██╔═══██╗██╔════╝   "<<std::endl;
            std::cout<<"   ██╔████╔██║███████║██║░░██║██║╚█████╗░██║██╗██║█████╗░░   "<<std::endl;
            std::cout<<"   ██║╚██╔╝██║██╔══██║██║░░██║██║░╚═══██╗╚██████╔╝██╔══╝░░   "<<std::endl;
            std::cout<<"   ██║░╚═╝░██║██║░░██║██████╔╝██║██████╔╝░╚═██╔═╝░███████╗   "<<std::endl;
            std::cout<<"   ╚═╝░░░░░╚═╝╚═╝░░╚═╝╚═════╝░╚═╝╚═════╝░░░░╚═╝░░░╚══════╝   "<<std::endl;
            std::cout<<"                                                             "<<std::endl;
            std::cout<<"                        Version 1.0.0                        "<<std::endl;
            std::cout<<"                     Author: D.Figueroa.                     "<<std::endl;
            std::cout<<"                                                             "<<std::endl;
            std::cout<<"#############################################################"<<std::endl;
            
            std::cout<<"#############################################################"<<std::endl;
            std::cout<<"                                                             "<<std::endl;
            std::cout<<"        MADISQE Instance Created                             "<<std::endl;
            std::cout<<"                                                             "<<std::endl;
            
            for (auto Setting : InputFile){
                std::cout<<"           "<< Setting.first << " : " << Setting.second   <<std::endl;
            }

            std::cout<<"                                                             "<<std::endl;
            std::cout<<"#############################################################"<<std::endl;
            
            std::string PDFSet     = InputFile.at("LHAPDFSet");
            std::string Provider   = InputFile.at("OLP");
            std::string Integrator = InputFile.at("Integrator");
            
            MC_SP.Method    = InputFile.at("Method");
            MC_SP.MaxEval   = stoi(InputFile.at("NEvaluations"));
            MC_SP.NStart    = stoi(InputFile.at("NStart"));
            MC_SP.NIncrease = stoi(InputFile.at("NIncrease"));
            
            XS.Integrand = InputFile.at("Integrand");
            XS.Channel   = InputFile.at("Channel");
            XS.Coupling  = InputFile.at("Coupling");
                
            double sqrts = stod(InputFile.at("sqrts"));
            double muRen = stod(InputFile.at("muRen"));
            double muFac = stod(InputFile.at("muFac"));

            XSec_Int = new XSection_Integrator(Provider,PDFSet,Integrator);
            XSec_Int->XSec->SetScales(sqrts,muRen,muFac);     
            
        }

        ~Madisqe(){
            delete XSec_Int;
        }

        void Run(){
            XSec_Int->ComputeXSection(XS,MC_SP);
        }
};

int main(int argc, char* argv[]){

    if( argc < 2){
        std::cout<<"Error: No input file specified"<<std::endl;
        abort();
    }
    if( argc > 2){
        std::cout<<"Error: Too many arguments"<<std::endl;
        abort();
    }

    Madisqe E1(argv[1]);
    E1.Run();

    return 0;
}