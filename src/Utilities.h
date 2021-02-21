#ifndef __UTILITIES_H__
#define __UTILITIES_H__

#include <time.h>
#include <chrono>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <vector>
#include <cmath>

#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/wait.h>
#include <unistd.h>

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
                FMatrixT<T> Bmunu(ccsc[cc]);
                std::cout<<"B("<<cc<<")= \n";
                std::cout<<Bmunu<<std::endl;
            }
        }
};

typedef ColorAndSpinMatrixT<double> ColorAndSpinMatrix;

class Clock{

    uint64_t start;

    public:

        uint64_t musec   = 1;
        uint64_t msec    = 1000*musec;
        uint64_t sec     = 1000*msec;
        uint64_t min     = 60*sec;
        uint64_t hour    = 60*min;
        uint64_t day     = 24*hour;

        uint64_t t_sequence[6] = {day,hour,min,sec,msec,musec};
        std::string n_sequence[6] = {"day","hour","min","sec","milisec","microsec"};

        Clock(){
            start = Time();
        }
        ~Clock(){}
        
        uint64_t Time() {
            return std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
        }


        void Reset(){
            start = Time();
        }

        void ShowTime(std::ostream& out = std::cout){

            uint64_t t = Time()-start;

            for(int i=0;i<6;i++){
                int t_i = t/t_sequence[i];
                if(t_i)std::cout<<t_i<<" "<<n_sequence[i]<<(t_i == 1 ? " " : "s ");
                t -= t_i*t_sequence[i];
            }

            out<<std::endl;
        }

        std::string GetTime(){

            uint64_t t = Time()-start;
            std::string rval;
            for(int i=0;i<6;i++){
                int t_i = t/t_sequence[i];
                if(t_i) rval = rval+std::to_string(t_i)+" "+n_sequence[i]+(t_i == 1 ? " " : "s ");
                t -= t_i*t_sequence[i];
            }
            return rval;
        }          
};

template<class T>
class PoleVector{
    
    public:
        
        T DoublePole;
        T SinglePole;
        T FinitePart;

        PoleVector(){


        }

        ~PoleVector(){

        }

};

template<class T>
class HistogramT{

    struct Bin{
        T LOW;
        T HIGH;
        T ACUM;
        size_t COUNT;
    };

    std::string Name;
    std::vector<Bin> Bins;
    
    public:

        static unsigned int HeavisideTheta(T Value){
            return not std::signbit(Value);
        }

        HistogramT(const std::string N, const std::vector<T> BinBoundaries = {}){
            Name = N;
            for(unsigned int i=0;i<(BinBoundaries.size()-1);i++){
                Bin NewBin;
                NewBin.LOW = (BinBoundaries.at(i));
                NewBin.HIGH = (BinBoundaries.at(i+1));
                NewBin.ACUM = 0;
                NewBin.COUNT = 0;
                Bins.push_back(NewBin);
            }

        }

        ~HistogramT(){
            Print();
        };

        void Append(const T Value, const T Result){
            unsigned int Belongs;
            for (auto& Bin : Bins){
                Belongs = HeavisideTheta(Bin.HIGH-Value)*HeavisideTheta(Value-Bin.LOW);
                Bin.ACUM  += Belongs*Result;
                Bin.COUNT += Belongs;
            }
        }

        void Print(){
            std::string HName = Name;
            std::ofstream out(HName);
            int GLOBAL = 0;
            for(auto Bin : Bins){
                out<<(Bin.LOW+Bin.HIGH)/2.0<<" "<<(Bin.COUNT==0?0:Bin.ACUM/Bin.COUNT)<<std::endl;
                GLOBAL += Bin.COUNT;
            }
            out<<"# Total number of recorded events : "<<GLOBAL<<std::endl;
        }

        void Show(){
            for(auto Bin : Bins){
                std::cout<<(Bin.LOW+Bin.HIGH)/2.0<<" "<<(Bin.COUNT==0?0:Bin.ACUM/Bin.COUNT)<<std::endl;
            }
        }
};

#define DIV(a,b) ((b)==0?0:((a)/(b)))

template<class T>
class SharedHistogramT{

    static const size_t _max_nbins = 300;

    static unsigned int HeavisideTheta(double Value){
            return not std::signbit(Value);
    }
    
    public:

        T Boundaries[_max_nbins];
        struct Bins{
            T AVG[_max_nbins];
            T VAR[_max_nbins];
            size_t CNT;
        };

        std::string XAxis;
        size_t NBins;
        Bins * PTR;
        key_t Key;;
        int ID;
        int PID = 1245;

        SharedHistogramT(const std::string VarName ,const std::vector<T> BinBoundaries){
            Key = ftok(VarName.c_str(),PID);
            ID = shmget(Key,sizeof(struct Bins),IPC_CREAT|0777);
            PTR = (struct Bins*)shmat(ID,NULL,0);
            XAxis = VarName;
            NBins = BinBoundaries.size()-1;
            if(NBins<=_max_nbins) PTR = new Bins[NBins];
            else{
                std::cout<<"Error: Histogram is too large, the current limit on the number of bins is:"<<_max_nbins<<std::endl;
                throw "Histogram too large";
            }
            for(unsigned int i=0;i<BinBoundaries.size();i++) Boundaries[i] = BinBoundaries[i];
        }

        ~SharedHistogramT(){
            shmctl(ID,IPC_RMID,NULL);
        }

        void Append(T Value, T Result){
            PTR = (struct Bins*)shmat(ID,NULL,0);
            unsigned int Belongs;
            PTR->CNT += 1;
            T Delta;
            for (size_t i=0;i<NBins;i++){
                Delta = Result-PTR->AVG[i];
                Belongs = HeavisideTheta(Boundaries[i+1]-Value)*HeavisideTheta(Value-Boundaries[i]);
                PTR->AVG[i] += Belongs*DIV(Delta,PTR->CNT);
                PTR->VAR[i] += Belongs*DIV(Delta*Delta,(PTR->CNT*(PTR->CNT-1)));
            }
            shmdt(PTR); 
        }

        void Write(std::string Title){
            std::string Path = Title+XAxis+".grc";
            std::ofstream out(Title);
            std::cout<<"Writing Histogram for: "<<Title<<" -> "<<XAxis<<std::endl;
            out<<"# Grace project file"<<std::endl;
            out<<"@    title \""<<Title<<"\""<<std::endl;
            out<<"@    title font 10"<<std::endl;
            out<<"@    xaxis  label \""<<XAxis<<"\""<<std::endl;
            out<<"@    xaxis  label font 10"<<std::endl;
            out<<"@    yaxes scale Logarithmic"<<std::endl;
            out<<"@    yaxis  label \"X-Section (pb)\""<<std::endl;
            out<<"@    yaxis  label font 10"<<std::endl;
            out<<"@    s0 line type 3"<<std::endl;
            out<<"@target G0.S0"<<std::endl;
            out<<"@type xydy"<<std::endl;
            PTR = (struct Bins*)shmat(ID,NULL,0);
            T XSec = 0;
            for(size_t i=0;i<NBins;i++){
                XSec += PTR->AVG[i];
                out<<(Boundaries[i]+Boundaries[i+1])/2.0<<" "<<PTR->AVG[i]<<" "<<sqrt(PTR->VAR[i])<<std::endl;
            }
            out<<"&"<<std::endl;
            out<<"Accumulated Xsec = "<<XSec<<std::endl;
            shmdt(PTR); 
        }
};

typedef SharedHistogramT<double> Histogram;

#endif
