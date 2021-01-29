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

class Model : public Input{

    public:

        double alpha_e = 1.0;
        double alpha_s = 1.0;
        
        const int NParticles = ####NParticles####;
    
####Build Model####

        Model(std::string InputFileName){

            LoadInput(InputFileName,InputFile);
            std::cout<<"Model Mass Environment Initialized"<<std::endl;
            
            for (auto Setting : InputFile){
                std::cout<<Setting.first<<" = "<<Setting.second<<std::endl;
            }

            
            
        }

        ~Model(){
            
        }

};

#endif