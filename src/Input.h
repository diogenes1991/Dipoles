#ifndef __INPUT_H__
#define __INPUT_H__

#include <unordered_map>
#include <iostream>
#include <fstream>
#include <sstream>

class Input{

    public:

        static void LoadInput(const std::string & filename, std::unordered_map<std::string,std::string>& settings){
            
            // There is a bug in here were the last line is not read 
            // for now we get around it by adding an extra inert line 
            // at the end of the Input file

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

        std::unordered_map<std::string,std::string> InputFile;
        virtual ~Input() {};
};

#endif