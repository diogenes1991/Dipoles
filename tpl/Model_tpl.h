#ifndef __MODEL_H__
#define __MODEL_H__

class Particle{

    public:

        std::string Name;
        double Mass;
        int PID;

        Particle(std::string N, double m, int p){
            Name = N;
            Mass = m;
            PID  = p;
        }

        ~Particle(){

        }

};

class Model{

    public:

        double alpha_e = 1.0;
        double alpha_s = 1.0;
        
        const int NParticles = ####NParticles####;
    
####Build Model####

        Model(){

        }

        ~Model(){
            
        }

};

#endif