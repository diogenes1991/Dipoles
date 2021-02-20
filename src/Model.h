#ifndef __MODEL_H__
#define __MODEL_H__

#include <string>

class Particle{

    public:

        std::string Name;
        double Mass;
        double Width;
        int PID;

        Particle(std::string N, double m, double w, int p){
            Name = N;
            Mass = m;
            PID  = p;
            Width = w;
        }

        ~Particle(){

        }

};

class Model{

    const int NParticles = 30;

    public:

        int NLF = 4;
        int NHF = 2;

        bool UseCMScheme = true;
        bool UseGMuScheme = true;

        double alpha_e = 1.0;
        double alpha_s = 1.0;
        double GFermi  = 1.16637E-5;
        
        Particle em = Particle("em",0.0,0.0,11);
        Particle sbar = Particle("sbar",0.0,0.0,-3);
        Particle ntbar = Particle("ntbar",0.0,0.0,-16);
        Particle dbar = Particle("dbar",0.0,0.0,-1);
        Particle nmbar = Particle("nmbar",0.0,0.0,-14);
        Particle ep = Particle("ep",0.0,0.0,-11);
        Particle nm = Particle("nm",0.0,0.0,14);
        Particle tbar = Particle("tbar",0.0,0.0,-6);
        Particle Wm = Particle("Wm",0.0,0.0,-24);
        Particle ne = Particle("ne",0.0,0.0,12);
        Particle bbar = Particle("bbar",0.0,0.0,-5);
        Particle Wp = Particle("Wp",0.0,0.0,24);
        Particle mum = Particle("mum",0.0,0.0,13);
        Particle u = Particle("u",0.0,0.0,2);
        Particle nt = Particle("nt",0.0,0.0,16);
        Particle A = Particle("A",0.0,0.0,7);
        Particle s = Particle("s",0.0,0.0,3);
        Particle Z = Particle("Z",0.0,0.0,23);
        Particle mup = Particle("mup",0.0,0.0,-13);
        Particle c = Particle("c",0.0,0.0,4);
        Particle b = Particle("b",0.0,0.0,5);
        Particle d = Particle("d",0.0,0.0,1);
        Particle g = Particle("g",0.0,0.0,0);
        Particle h = Particle("h",0.0,0.0,25);
        Particle taum = Particle("taum",0.0,0.0,15);
        Particle taup = Particle("taup",0.0,0.0,-15);
        Particle nebar = Particle("nebar",0.0,0.0,-12);
        Particle ubar = Particle("ubar",0.0,0.0,-2);
        Particle t = Particle("t",0.0,0.0,6);
        Particle cbar = Particle("cbar",0.0,0.0,-4);

        Model();

        ~Model(){
            
        }

};

#endif