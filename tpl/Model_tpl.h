#ifndef __MODEL_H__
#define __MODEL_H__

#include <string>

class Particle{

    public:

        std::string Name;
        double Mass;
        double Width;
        double Charge;
        int PID;

        Particle(std::string N, double m, double w, double c, int p){
            Name = N;
            Mass = m;
            PID  = p;
            Charge = c;
            Width = w;
        }

        ~Particle(){

        }

};

class Model {

    const int NParticles = ####NParticles####;

    public:

        int NLF = 4;
        int NHF = 2;

        bool UseCMScheme = true;
        bool UseGMuScheme = true;

        double alpha_e = 1.0;
        double alpha_s = 1.0;
        double GFermi  = 1.16637E-5;
    
####Build Model####

        Model();

        ~Model(){}

};


#endif