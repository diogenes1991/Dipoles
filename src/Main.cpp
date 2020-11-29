#include "XSection_Integrator.h"

void XSectionStart(std::string pdfname){

    X1 = new XSection(pdfname);
}

void XSectionEnd(){

    delete X1;
}

void LoadInput(const std::string & filename, std::unordered_map<std::string,std::string>& settings){
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

int main(int argc, char* argv[]){

    if( argc < 2){
        std::cout<<"Error: No input file specified"<<std::endl;
        abort();
    }
    if( argc >2){
        std::cout<<"Error: Too many arguments"<<std::endl;
        abort();
    }

    std::unordered_map<std::string,std::string> Settings;
    LoadInput(argv[1],Settings);
            
    std::string PDFSet = Settings.at("LHAPDFSet");

    std::string Integrator = Settings.at("Integrator");
    std::string Method = Settings.at("Method");

    Montecarlo_Integrator::Specifications MC_SP;
    MC_SP.MaxEval = stoi(Settings.at("NEvaluations"));
    MC_SP.NStart = stoi(Settings.at("NStart"));
    MC_SP.NIncrease = stoi(Settings.at("NIncrease"));

    XSection_Integrator::XSection_Selector XS;
    XS.Integrand = Settings.at("Integrand");
    XS.Channel = Settings.at("Channel");
    XS.Coupling = Settings.at("Coupling");

    XSectionStart(PDFSet);
    X1->SetScales(13000,91.1876,91.1876);
    XSection_Integrator C1(Integrator);
    C1.ComputeXSection(XS,MC_SP,Method);
    XSectionEnd();

    return 0;
}