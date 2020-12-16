#ifndef __UTILITIES_H__
#define __UTILITIES_H__

#include <time.h>

template <class T>
class ColorAndSpinMatrixT{

    public:
        T*** ccsc;
        int n;
        ColorAndSpinMatrixT(int NCC){
            n = NCC;
            ccsc = new T**[n];
            for(int cc=0;cc<n;cc++){
                ccsc[cc] = new T*[4];
                for(int mu=0;mu<4;mu++){
                    ccsc[cc][mu] = new T[4];
                }
            }
        }
        ~ColorAndSpinMatrixT(){
            for(int cc=0;cc<n;cc++){
                for(int mu=0;mu<4;mu++){
                    delete [] ccsc[cc][mu];
                }
                delete [] ccsc[cc];
            }
            delete [] ccsc;
        }
        void Show(){
            for(int cc=0;cc<n;cc++){
                FourMatrixT<T> Bmunu(ccsc[cc]);
                std::cout<<"B("<<cc<<")= \n";
                std::cout<<Bmunu<<std::endl;
            }
        }
};

typedef ColorAndSpinMatrixT<double> ColorAndSpinMatrix;

class Clock{

    size_t start;

    public:

        Clock(){
            start = time(NULL);
        }
        ~Clock(){}
        
        void Reset(){
            start = time(NULL);
        }

        void ShowTime(){

            size_t t = time(NULL)-start;

            size_t sec = 1;
            size_t min = 60*sec;
            size_t hour = 60*min;
            size_t day = 24*hour;

            size_t t_sequence[4] = {day,hour,min,sec};
            std::string n_sequence[4] = {"day","hour","min","sec"};
            
            for(int i=0;i<4;i++){
                int t_i = t/t_sequence[i];
                if(t_i)std::cout<<t_i<<" "<<n_sequence[i]<<(t_i == 1 ? " " : "s ");
                t -= t_i*t_sequence[i];
            }

            std::cout<<std::endl;
        }

            
};

#endif
