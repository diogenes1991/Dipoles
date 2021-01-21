#ifndef __NLOX_OLP_H__
#define __NLOX_OLP_H__

#include "nlox_process.h"

class NLOX_OLP : public OLP{

    public:

        Process * Proc;
        std::unordered_map<std::string,std::string> ChannelMap;

        NLOX_OLP(){
            Proc = new Process();

            ####BuildChannels####
        }

        ~NLOX_OLP(){
            delete Proc;
        }

        void Evaluate(Arguments * arg){
            std::string subproc = SelectSubProc(arg.SubProc);
            evaluate(subroc,arg)
        }
};

#endif
