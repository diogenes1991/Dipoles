
import sys
import os
from Dipole_Headers import *

def find_substring(string, substring):
    substring_length = len(substring)    
    def recurse(locations_found, start):
        location = string.find(substring, start)
        if location != -1:
            return recurse(locations_found + [location], location+substring_length)
        else:
            return locations_found
    return recurse([], 0)


try:
    fconfig = open(sys.argv[1],'r')
except IndexError:
    print "Error: No input file specified"
    quit()
except IOError: 
    print "Error: No input file",sys.argv[1],"found"
    quit()






#   PARTIICLE IDs FETCHING ROUTINE  #

mystr = fconfig.readline()
INIDS = []
mystr = mystr.strip()
start = mystr.find('=')
ends = find_substring(mystr,",")
ends = [start] + ends
ends.append(len(mystr))

for i in range(0,len(ends)-1):
    toappend = ""
    for j in range(ends[i]+1,ends[i+1]):
        toappend = toappend + mystr[j]
    INIDS.append(toappend)
INIDS = map(int,INIDS)


mystr = fconfig.readline()
FIIDS = []
mystr = mystr.strip()
start = mystr.find('=')
ends = find_substring(mystr,",")
ends = [start] + ends
ends.append(len(mystr))

for i in range(0,len(ends)-1):
    toappend = ""
    for j in range(ends[i]+1,ends[i+1]):
        toappend = toappend + mystr[j]
    FIIDS.append(toappend)
FIIDS = map(int,FIIDS)

mystr = fconfig.readline()
if mystr.find("0")>0: VERBOSE = 0
elif mystr.find("1")>0: VERBOSE = 1
elif mystr.find("2")>0: VERBOSE = 2




#   C++ INTEGRAND FILES CREATION  #



EXTIDS = INIDS + FIIDS + [22]
NEXT = len(EXTIDS)

PROCESSNAME = PID_to_Name(EXTIDS[0])+PID_to_Name(EXTIDS[1])+'_'
for ip in range(2,len(EXTIDS)-1):
    PROCESSNAME = PROCESSNAME + PID_to_Name(EXTIDS[ip])


PROCESSNAME_SUB = PROCESSNAME+"_Subtracted.cpp"
FSUB = open (PROCESSNAME_SUB,"w+")
PROCESSNAME_PLU = PROCESSNAME+"_Plus.cpp"
FPLU = open (PROCESSNAME_PLU,"w+")
PROCESSNAME_END = PROCESSNAME+"_Endpoint.cpp"
FEND = open (PROCESSNAME_END,"w+")


def_vec_args = ""
cal_vec_args = ""
for i in range(1,NEXT+1):
    def_vec_args = def_vec_args +"Vector p"+str(i)
    cal_vec_args = cal_vec_args +"p"+str(i) 
    if i!=NEXT: 
        def_vec_args = def_vec_args + ","
        cal_vec_args = cal_vec_args + ","

line = "double "+PROCESSNAME+"_Subtracted"+"("+def_vec_args+"){\n"
FSUB.write(line)
line = "double out = 0;\n"
FSUB.write(line)
line = "Vector p_data["+str(NEXT)+"]={"+cal_vec_args+"};\n"
FSUB.write(line)
line = "Vector p_tilde["+str(NEXT-1)+"];\n"
FSUB.write(line)


line = "double "+PROCESSNAME+"_Part_Int"+"("+def_vec_args+"){\n"
FPLU.write(line)
line = "double out = 0;\n"
FPLU.write(line)

line = "double "+PROCESSNAME+"_Full_Int"+"("+def_vec_args+"){\n"
FEND.write(line)
line = "double out = 0;\n"
FEND.write(line)


if VERBOSE:
    for i in range(0,NEXT):
        if i==0: 
            print "QED dipole routine triggered for the process ",
        if i==2: print " ---> ",
        if i==NEXT-1: print '[',
        line = PID_to_Name(EXTIDS[i])+"("+str(i+1)+")"
        print line,
        if i==NEXT-1: print ']',
        if i!=NEXT-1 and i!=1: print " + ",
        elif i!=1: print
            

for i in range(0,NEXT):
        for j in range(0,NEXT):
            ch = Charge(abs(EXTIDS[i]))*Charge(abs(EXTIDS[j]))*signo(EXTIDS[i])*signo(EXTIDS[j])
            sym_ch = "("+str(signo(EXTIDS[i])*signo(EXTIDS[j]))+")*("+Symb_Charge(abs(EXTIDS[i]))+")*("+Symb_Charge(abs(EXTIDS[j]))+")"
            if ch!=0 and j!=i:
                conf = 0;
                if i<2 and j<2: conf = 0
                if i<2 and j>=2: conf = 1
                if i>=2 and j<2: conf = 2
                if i>=2 and j>=2: conf = 3
                dipconf = ["II","IF","FI","FF"]
                dipsubconf = ["fermion","boson"]
                if abs(EXTIDS[i])!=24:
                    subconf = 0
                else:
                    subconf = 1
            
                if(VERBOSE):
                    print "Including ", dipconf[conf] , " dipole with ",
                    line = PID_to_Name(EXTIDS[i]) + "(" + str(i+1) + ") as emitter and "
                    print line,
                    line = PID_to_Name(EXTIDS[j]) +"(" + str(j+1) + ") as spectator"
                    print line
                SUBNAMES = ["g_ab","g_ai","g_ia","g_ij"]
                PLUNAMES = ["CurlyG_ab","CurlyG_ai","CurlyG_ia","CurlyG_ij"]
                ENDNAMES = ["G_ab","G_ai","G_ia","G_ij"]
                
                
                line = "Build_"+dipconf[conf]+"_Momenta("+str(NEXT-1)+","+"p_data,p_tilde,"+str(i+1)+","+str(j+1)+");\n"
                FSUB.write(line)
                line = "out = out + ("+sym_ch+")*"+SUBNAMES[conf] + "_" + dipsubconf[subconf] +  "(p"+str(i+1)+",p"+str(j+1)+",p"+str(NEXT)+")"+"*"+PROCESSNAME+"(p_tilde);\n"
                FSUB.write(line)
                
                
                
                
                line = "out = out + ("+sym_ch+")*"+PLUNAMES[conf] + "_" + dipsubconf[subconf] + "(p"+str(i+1)+",p"+str(j+1)+",x"+"_"+str(i+1)+str(j+1)+");\n"
                FPLU.write(line)
                line = "out = out + ("+sym_ch+")*"+ENDNAMES[conf] + "_" + dipsubconf[subconf] + "(p"+str(i+1)+",p"+str(j+1)+");\n"
                FEND.write(line)
                
                
line = "return " + PROCESSNAME + "A(p_data) + out;\n"        
FSUB.write(line)            
FSUB.write("}\n\n")
FSUB.close()
print 'Subtracted integrand built at',PROCESSNAME_SUB


FPLU.write("return out;\n}")            




FPLU.close()
print 'Plus Distribution integrand built at',PROCESSNAME_PLU




FEND.write("return out;\n")            
FEND.write("}\n")
line = "double "+PROCESSNAME+"_Endpoint"+"("+def_vec_args+"){\n"
FEND.write(line)
line = "double out = "+PROCESSNAME+'_Full_Int('+cal_vec_args+')'+";\n"
FEND.write(line)
line = "out = out*"+PROCESSNAME+"("+cal_vec_args+");\n"
FEND.write(line)
FEND.write('return out;\n}                                \n')
FEND.close()
print 'Endpoint integrand built at',PROCESSNAME_END

print 'The matix elements must be provided as c++ callable objects',
line = PROCESSNAME+'.o'
print line,
line = 'and '+PROCESSNAME+'A.o'
print line




#   MAIN C++ FILE CREATION   #

MAIN_NAME = PROCESSNAME+"_Main.cpp"
FMAIN = open(MAIN_NAME,"w+")

FMAIN.write('#include <iostream>\n')
FMAIN.write('#include <fstream>\n')
FMAIN.write('#include <cmath>\n')
FMAIN.write('#include "Standard_Model.h"\n')
FMAIN.write('#include "Vector_Real.h"\n')
FMAIN.write('#include "Phase_Spaces_Real.h"\n')
FMAIN.write('#include "Dipole_Definitions.h"\n')
line = '#include "'+PROCESSNAME_SUB+'"\n'
FMAIN.write(line)
#line = '#include "'+PROCESSNAME_PLU+'"\n'
#FMAIN.write(line)
#line = '#include "'+PROCESSNAME_END+'"\n'
#FMAIN.write(line)

FMAIN.write('\n\nmain(int argc, char* argv[]){        \n\n')        
FMAIN.write('Vector p1 (5,3,4,0);                       \n')
FMAIN.write('Vector p2 (13,12,0,5);                     \n')
FMAIN.write('Vector p3 (12,11,5,1);                     \n')
FMAIN.write('Vector p4 (10,8,1,1);                      \n')
FMAIN.write('Vector p5 (17,8,15,0);                     \n')
FMAIN.write('Vector p6 (7,6,1,1);                       \n')
FMAIN.write('                                           \n')
FMAIN.write('cout.precision(16);                        \n')
FMAIN.write('                                           \n')
line = 'cout << ' + PROCESSNAME + '_Subtracted('
for i in range(0,NEXT): 
    line = line + 'p'+str(i+1)
    if i!=NEXT-1: line = line + ','
line = line + ') << endl;\n'
FMAIN.write(line)
FMAIN.write('\n\n}                                      \n')

FMAIN.close()




# NLOX Configuration files creation               #
# NOTE: THIS CREATES THE LEADING ALPHA CORRECTION #


NLOXC = PROCESSNAME+'_NLOX.dat'
FNLOX = open(NLOXC,'w+')
FTEMP = open('NLOX_Configuration.Template','r') 
        
for j in FTEMP:
    i = ''
    for k in range(0,len(j)-1):
        i = i+j[k]
    if i.find('initialState')==0:
        i = i + PID_to_Name(EXTIDS[0])+','+PID_to_Name(EXTIDS[1])
    if i.find('finalState')==0: 
        for ip in range(2,len(EXTIDS)-1):
            i = i + PID_to_Name(EXTIDS[ip])
            if ip!=len(EXTIDS)-2: i = i+','
    if i.find('alphaQCD')==0:
        i = i + '0'
    if i.find('alphaEW')==0:
        i = i + str(2+(len(EXTIDS)-4))
    if i.find('processPath')==0:
        i = i + PROCESSNAME
    i = i+'\n'
    FNLOX.write(i)
FTEMP.close()
FNLOX.close()
line = 'NLOX configuration file generated at '+NLOXC
print line


#  MAKEFILE CREATION # 


MKFLN = PROCESSNAME+'_Makefile'
FMAKE = open(MKFLN,'w+')
FTEMP = open('Makefile.Template','r') 
        
for j in FTEMP:
    i = ''
    for k in range(0,len(j)-1):
        i = i+j[k]
    #print i
    if i.find('PROCESS')==0:
        i = i + PROCESSNAME
    if i.find('SOURCE')==0: 
        i = i + PROCESSNAME+'_Main'
    if i.find('OUT')==0: 
        i = i + PROCESSNAME +'_Main_Test'
    if i.find('@nlox.py')>0: 
        i = i +' '+NLOXC+' > '+PROCESSNAME+'_NLOX.log'
    i = i+'\n'
    FMAKE.write(i)
FTEMP.close()
FMAKE.close()
line = 'Makefile generated at '+MKFLN
print line
#line = 'make -f '+MKFLN
#os.system(line)



#  FILE ACCUMULATION AND COMPILATION # 


line  = 'mkdir '+PROCESSNAME
os.system(line)
line = 'mv '+PROCESSNAME_END+' '+PROCESSNAME_PLU+' '+PROCESSNAME_SUB+' '
line = line +' '+MAIN_NAME+' '+NLOXC+' '+MKFLN+' '+ PROCESSNAME
os.system(line)
line = 'cp ' + sys.argv[1]+' '+PROCESSNAME
os.system(line)
















