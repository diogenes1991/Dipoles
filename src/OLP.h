#ifndef __OLP_H__
#define __OLP_H__

class OLP{

    std::unordered_map<std::string,int> SubProc;

    public:

        struct Arguments{
            std::string SubProc;
            std::string Order;
            std::string Coupling;
            double mu_ren;
            FourVector * P;
            double * Rval;
        };

        int SelectSubProcess(std::string SP){
            int SubProcess;
            try {SubProcess = ChannelMap.at(SP);}
            catch (const std::out_of_range& oor) {
                std::cerr<<"Error: Channel "<<SP<<" not found in Process"<<std::endl;
                std::cout<<"The available Radiative processes are:"<<std::endl;
                for ( auto& x : SubProc ) std::cout<<x.first<<" => "<<x.second<<std::endl;
                abort();
            }
            return SubProcess;
        }

        virtual ~OLP(){};
        virtual void Evaluate(Arguments * args) = 0;

};

#endif
