#include "XSection.h"

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

    std::unordered_map<std::string,std::string> Settings;

    LoadInput(argv[1],Settings);
    
    std::string PDFSet = Settings.at("LHAPDFSet");

    std::string Integrator = Settings.at("Integrator");
    std::string Method = Settings.at("Method");

    XSecCalc::xsec_sel XSC;

    std::string Integrand = Settings.at("Integrand");
    std::string Channel = Settings.at("Channel");
    std::string Coupling = Settings.at("Coupling");
    

    XSectionStart(PDFSet);
    X1->SetScales(13000,91.1876,91.1876);
    XSecCalc C1(Integrator); // Careful one can set up for seg-fault by not initializong an integrator but calling it anyways
    C1.ComputeXSections(Integrand,Channel,Coupling,Integrator,Method);
    XSectionEnd();

    return 0;
}