#include "XSection_Integrator.h"

class Environment{

    std::unordered_map<std::string,std::string> InputFile;

    static void LoadInput(const std::string & filename, std::unordered_map<std::string,std::string>& settings){
        std::ifstream data;
        std::string linebuf;
        data.open(filename.data());
            if (!data.fail()) {
                while(!std::getline(data,linebuf).eof()) {
                    if(linebuf=="") continue;
                    std::istringstream ss(linebuf);
                    std::string name;
                    char equal;
                    std::string val;
                    ss >> name;
                    ss >> equal;
                    if (!ss.fail() && equal == '=') {
                        ss >> val;
                        settings.insert({name,val});
                    }
                    else{
                        std::cout<<"Warning: Malformed line at input"<<linebuf<<std::endl;
                    }
                }
            data.close();
            }
            else {
                std::cout << "Error: No input file found" << std::endl;
                abort();
            }
    }

    XSection_Integrator * XSec_Int;

    Montecarlo_Integrator::Specifications MC_SP;
    XSection_Integrator::XSection_Selector XS;

    public:

        Environment(std::string InputFileName){
            LoadInput(InputFileName,InputFile);
            std::string PDFSet     = InputFile.at("LHAPDFSet");
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
            
            XSec_Int = new XSection_Integrator(PDFSet,Integrator);
            XSec_Int->XSec->SetScales(sqrts,muRen,muFac);
        }

        ~Environment(){
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
    if( argc >2){
        std::cout<<"Error: Too many arguments"<<std::endl;
        abort();
    }

    Environment E1(argv[1]);
    E1.Run();

    return 0;
}