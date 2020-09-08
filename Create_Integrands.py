#!/usr/bin/env python

from Particle_Engine import *
from Utilities import *

class Amplitude:
    def __init__(self,CONFIGFILE):
        self.ReadConfig(CONFIGFILE)
        self.LoadConfig()
        self.LinkProcs()
        self.ReflectLink()
        self.Build()

    def ReadConfig(self,CONFIGFILE):
        CONFIG = {'INITIAL STATE': [0],\
                  'FINAL STATE'  : [0],\
                  'RADIATIVE'    : [0],\
                  'MODELFILE'    : ['StandardModel'],\
                  'PATH'         : [os.getcwd()],\
                  'NLOX PATH'    : [0] ,\
                  'VERBOSE'      : [0] }
        
        REQUIRED = set(['INITIAL STATE','FINAL STATE','NLOX PATH'])

        COMMCHAR = '#'

        f = open(CONFIGFILE,'r')
        configfile= []
        for l in f:
            if l.find(COMMCHAR)==0 or len(l)==0 or l=='\n':
                continue
            configfile.append([False,l])

        for configuration in CONFIG.keys():
            warn = True
            for line in configfile:
                if configuration in line[1]:
                    thisconfig = line[1].strip().split('=')
                    CONFIG[thisconfig[0].strip()] = (thisconfig[1].strip()).split(',')
                    warn = False
                    line[0] = True
            if warn:
                if configuration not in REQUIRED:
                    print "\33[95mWarning\33[0m: Configuartion",configuration,"not found in configuration file, using default:",CONFIG[configuration][0]
                else:
                    print "\33[31mError\33[0m: Configuration",configuration,"not found in configuration, this is a required configuration."
                    sys.exit()

        
        for line in configfile:
            if not line[0]:
                print "\33[95mWarning\33[0m: Configuartion",(((line[1].strip()).split('='))[0]).strip(),"is not a valid configuration setting, please check the syntax"
        f.close()
        self.config = CONFIG

    def LoadConfig(self):
        
        self.verbose = int(self.config['VERBOSE'][0])
        self.nlox = self.config['NLOX PATH'][0]+'/nlox.py'

        try:
            LoadedModel = __import__(self.config['MODELFILE'][0])
        except:
            print "\33[31mError\33[0m: Model file",self.config['MODELFILE'][0],"not found."
            sys.exit()
                
        self.Model =  LoadedModel.WorkingModel
        self.ParticleContent = self.Model.ParticleContent

        self.ini = []
        for particle in self.config['INITIAL STATE']:
            try:
                self.ini.append(self.ParticleContent[particle])
            except:
                print "\33[31mError\33[0m: Particle class",particle,"not defined."
                sys.exit()
                    
        self.fin = []
        self.finrad = []
        for particle in self.config['FINAL STATE']:
            try:
                self.fin.append(self.ParticleContent[particle])
                self.finrad.append(self.ParticleContent[particle])
            except:
                print "\33[31mError\33[0m: Particle class",particle,"not defined."
                sys.exit()
        
        try:    
            self.finrad.append(self.ParticleContent[self.config['RADIATIVE'][0]])
        except:
            print "\33[31mError\33[0m: Particle class",self.config['RADIATIVE'][0],"not defined"
            sys.exit()

        self.BornProc = Process(self.ini,self.fin)
        self.RadiProc = Process(self.ini,self.finrad)
        
        self.path = self.config['PATH'][0]+"/"+str(self.BornProc)
        if os.path.exists(self.path):
            print "\33[31mError\33[0m: Source directory",self.path,"already exists, please remove it first"
            sys.exit()

    def LinkProcs(self):
        self.link = {}
        for subproc in self.RadiProc.subproc:
            self.link[subproc] = []
            counter = -1
            for particle in self.RadiProc.subproc[subproc]:
                counter += 1
                
                #   Basically I need to append a sublist of the process, I should 
                #   make use of the body of the function to copy the subprocess 
                #   only after striping the undesired particle 
                #
                #   If the Boson was found in the initial state we loop over all
                #   possible swappings that send it to the final state and then 
                #   remove it

                if particle.typ != 'Boson':
                    continue

                WHICH = [counter]

                # Otherwise we eliminate the Boson from the process and build a newlist 

                newlist = []
                for index in range(len(self.RadiProc.subproc[subproc])):
                    if index == counter:
                        continue
                    newlist.append(self.RadiProc.subproc[subproc][index])
                
                if counter >= len(self.RadiProc.ini):
                    # print "   The Boson",str(particle),"will be eliminated from the final state"
                    tentborn = self.BornProc.BuildString(newlist)
                    if tentborn in self.BornProc.subproc.keys() and tentborn not in self.link[subproc]:
                        self.link[subproc].append(tentborn)
                else:
                    # print "   The Boson",str(particle),"will be eliminated from the initial state"
                    for index in range(len(self.RadiProc.ini),len(self.RadiProc.ini)+len(self.RadiProc.fin)-1):
                        # Copy 
                        swaped = [particle for particle in self.RadiProc.subproc[subproc]]
                        
                        # Swap 
                        # print "      Swapping",swaped[counter],"<->",swaped[index],"(",AMPLITUDE.Model.CrossDictionary[swaped[index].nam],")"
                        # print "      Poping",swaped[counter]
                        swaped[index] = self.ParticleContent[self.Model.CrossDictionary[swaped[index].nam]]
                        swaped[counter],swaped[index]=swaped[index],swaped[counter]
                        
                        # Copy, but poped
                        newlist = []
                        for i in range(len(self.RadiProc.subproc[subproc])):
                            if i == index:
                                continue
                            newlist.append(swaped[i])
                        
                        # Resort the initials
                        sortedin = []
                        for i in range(len(self.RadiProc.ini)):
                            sortedin.append(newlist[i])
                        def MyF(p):
                            return -p.pid
                        sortedin.sort(key=MyF)
                        for i in range(len(self.RadiProc.ini)):
                            newlist[i] = sortedin[i]

                        tentborn = self.RadiProc.BuildString(newlist)
                        if tentborn in self.BornProc.subproc.keys() and tentborn not in self.link[subproc]:
                            self.link[subproc].append(tentborn)
                            # print "      The subprocess",tentborn,"\33[95mwas\33[0m found in the Born"
                        # else:
                            # print "      The subprocess",tentborn,"\33[95mwas not\33[0m found in the Born"            

    def ReflectLink(self):
        self.rlink = {}
        for rad in self.link:
            for born in self.link[rad]:
                if born in self.rlink:
                    self.rlink[born].append(rad)
                else:
                    self.rlink[born] = [rad]

    def Build(self):

        ###
        ###   This is the part of the code that will 
        ###   handle templates, creates and moves files
        ### 


        self.TplDir = 'tpl/'
        self.SrcDir = self.path
        self.MatDir = self.SrcDir+'/Matrix_Elements'

        MakeDir(self.SrcDir)
        self.BuildNLOXProcess()  

        self.BuildSeeds()
        # self.BuildMatrixElements()
        self.BuildSubProc()
        self.BuildDipoles()

    def BuildNLOXProcess(self):
        RadiDict = { 'Include Born'   : '' ,\
                     'Include Radi'   : '' ,\
                     'nBorn'          : str(len(self.BornProc.subproc)) ,\
                     'nRadi'          : str(len(self.RadiProc.subproc)) ,\
                     'Construct Born' : '\n' ,\
                     'Construct Radi' : '\n'  }
        TAB3 = '   '
        TAB6 = TAB3+TAB3
        TAB9 = TAB6+TAB3
        
        count = 0
        for subproc in self.BornProc.scatters:
            SubProc = self.BornProc.BuildString(subproc)
            RadiDict['Include Born'] += '#include "../'+SubProc+'/code/Sub_'+SubProc+'.h"\n'
            RadiDict['Construct Born'] += TAB9+'subproc['+str(count)+'] = '+'new Sub_'+SubProc+'(pc);\n'
            # RadiDict['Construct Born'] += TAB9+'BornMap.insert({"'+SubProc+'",'+str(count)+'});\n'
            count += 1

        for subproc in self.RadiProc.scatters:
            SubProc = self.RadiProc.BuildString(subproc)
            RadiDict['Include Radi'] += '#include "../'+SubProc+'/code/Sub_'+SubProc+'.h"\n'
            RadiDict['Construct Radi'] += TAB9+'subproc['+str(count)+'] = '+'new Sub_'+SubProc+'(pc);\n'
            # RadiDict['Construct Radi'] += TAB9+'RadiMap.insert({"'+SubProc+'",'+str(count)+'});\n'
            count += 1
        
        ProcessHeader = seek_and_destroy(self.TplDir+"nlox_process.tpl",RadiDict)
        WriteFile(self.SrcDir+"/nlox_process.h",ProcessHeader)

    def BuildSeeds(self,DetectCP=False):
        cutpoint = len(self.BornProc.ini)
        for subproc in self.BornProc.subproc:
            DICT = {'Initial State'   : '' ,\
                    'Final State'     : '' ,\
                    'Tree CP'         : 'bornAlphaQCD = ' ,\
                    'Loop CP'         : 'virtAlphaQCD = ' ,\
                    'Path'            : '' }
            Ini = [ s.nam for s in self.BornProc.subproc[subproc][:cutpoint]]
            Fin = [ s.nam for s in self.BornProc.subproc[subproc][cutpoint:]]
            count = 1
            for parname in Ini:
                DICT['Initial State'] += parname
                if count != len(Ini):
                    DICT['Initial State'] += ','
                count += 1
            count = 1
            for parname in Fin:
                DICT['Final State'] += parname
                if count != len(Fin):
                    DICT['Final State'] += ','
                count += 1
            BornPow = len(Ini)+len(Fin)-1
            for i in range(BornPow+1):
                if i!=0:
                    DICT['Tree CP'] += str(BornPow-i)
                    if i != BornPow:
                        DICT['Tree CP'] += ','
                DICT['Loop CP'] += str(BornPow-i)
                if i != BornPow:
                    DICT['Loop CP'] += ','
            DICT['Path'] +=  self.MatDir
            SubProcSeed = seek_and_destroy(self.TplDir+"NLOX_seed.tpl",DICT)
            WriteFile(self.SrcDir+'/'+subproc+'.in',SubProcSeed)  


        cutpoint = len(self.RadiProc.ini)
        for subproc in self.RadiProc.subproc:
            DICT = {'Initial State'   : '' ,\
                    'Final State'     : '' ,\
                    'Tree CP'         : 'bornAlphaQCD = ' ,\
                    'Loop CP'         : 'virtAlphaQCD = ' ,\
                    'Path'            : '' }
            Ini = [ s.nam for s in self.RadiProc.subproc[subproc][:cutpoint]]
            Fin = [ s.nam for s in self.RadiProc.subproc[subproc][cutpoint:]]
            count = 1
            for parname in Ini:
                DICT['Initial State'] += parname
                if count != len(Ini):
                    DICT['Initial State'] += ','
                count += 1
            count = 1
            for parname in Fin:
                DICT['Final State'] += parname
                if count != len(Fin):
                    DICT['Final State'] += ','
                count += 1
            RadiPow = len(Ini)+len(Fin)
            for i in range(RadiPow+1):
                if i!=0:
                    DICT['Tree CP'] += str(RadiPow-i)
                    if i != RadiPow:
                        DICT['Tree CP'] += ','
            DICT['Loop CP'] = '##virtAlphaQCD'
            DICT['Path'] +=  self.MatDir
            SubProcSeed = seek_and_destroy(self.TplDir+"NLOX_seed.tpl",DICT)
            WriteFile(self.SrcDir+'/'+subproc+'.in',SubProcSeed)          

    def BuildMatrixElements(self):
        Here = os.getcwd()
        os.chdir(self.SrcDir)
        
        for subproc in self.BornProc.subproc:
            print 'Generating born and virtual process:',subproc
            cmd = self.nlox+' '+subproc +'.in > '+subproc+'.log'
            os.system(cmd)
            cmd = 'mv '+subproc+'.log '+self.MatDir+'/'+subproc
            os.system(cmd)

        for subproc in self.RadiProc.subproc:
            print 'Generating real radiative process:',subproc
            cmd = self.nlox+' '+subproc +'.in > '+subproc+'.log'
            os.system(cmd)
            cmd = 'mv '+subproc+'.log '+self.MatDir+'/'+subproc
            os.system(cmd)

        os.chdir(Here)

    def BuildSubProc(self):
        for subproc in self.RadiProc.subproc:
            DICT={'SubProcHeader'    : '' ,\
                  'SubProcName'      : '' ,\
                  'SubProcMat'       : '' ,\
                  'SubProcSub'       : '' ,\
                  'SubProcPlu'       : '//empty for now' ,\
                  'SubProcEnd'       : '//empty for now'}
            DICT['SubProcHeader'] = subproc.upper()
            DICT['SubProcName'] = subproc+'_Dipoles'
            DICT['SubProcMat'] = '#include "Matrix_Elements/'+subproc+'/code/Sub_'+subproc+'.h" \n'
            
            # Include the necesary Borns
            for linked in self.link[subproc]:
                DICT['SubProcMat'] += '#include "Matrix_Elements/'+linked+'/code/Sub_'+linked+'.h"\n'

            # Construction of the EWK Subtracted Function
            # This uses auto coupling-power detection
            
            EWKc = 0
            QCDc = 0

            for particle in self.RadiProc.subproc[subproc]:
                if particle in self.Model.EWKPars:
                    EWKc += 1
                if particle in self.Model.QCDPars:
                    QCDc += 1
            # print QCDc,'QCD particles and',EWKc,'EWK particles'

            MaxQCD = QCDc - 2
            MinQCD = self.RadiProc.nex - EWKc   

            MaxEWK = EWKc - 2
            MinEWK = self.RadiProc.nex - QCDc

            # print 'QCD range = [',MinQCD,',',MaxQCD,']'
            # print 'EWK range = [',MinEWK,',',MaxEWK,']'

            CPSQCD = []
            for val in range(MinQCD,MaxQCD+1):
                CPSQCD.append('as'+str(val)+'ae'+str(self.RadiProc.nex-val-2))

            CPSEWK = []
            for val in range(MinEWK,MaxEWK+1):
                CPSEWK.append('as'+str(self.RadiProc.nex-val-2)+'ae'+str(val))
            # print CPSQCD
            # print CPSEWK

            for CP in CPSQCD:
                DICT['SubProcSub'] += 'Subprocess* Radiative = proc->call(RadiProc('+CP+'));\n'



            
            ce = 1
            for emitter in self.RadiProc.subproc[subproc]:
                cs = 1
                for spectator in self.RadiProc.subproc[subproc]:
                    if 'Charge' in emitter.sym and 'Charge' in spectator.sym and cs!=ce:
                        CONF = ('I' if (ce<len(self.RadiProc.ini)) else 'F' ) + ('I' if (cs<len(self.RadiProc.ini)) else 'F' )
                        DICT['SubProcSub'] += 'Add '+emitter.nam+'('+str(ce)+') as emitter and '+spectator.nam+'('+str(cs)+') as spectator '+CONF+' Dipole\n'
                    cs += 1
                ce += 1


            


            
                # PlusDistrb += subproc+'(p)-'+linked+'(Maptype(p));\n'
                # Endpoint   += subproc+'()'


            SubFunctionH = seek_and_destroy(self.TplDir+'/Dipole_FunctionsH.tpl',DICT)
            SubFunctionC = seek_and_destroy(self.TplDir+'/Dipole_FunctionsC.tpl',DICT)
            WriteFile(self.SrcDir+'/'+subproc+'_Function.cpp',SubFunctionC)
            WriteFile(self.SrcDir+'/'+subproc+'_Function.h',SubFunctionH)

    def BuildDipoles(self):
        for born in self.rlink:
            print 'For',born,'we have',len(self.rlink[born]),'corrections'
            for rad in self.rlink[born]:
                fir = self.RadiProc.subproc[rad][self.RadiProc.lni:]
                inr = self.RadiProc.subproc[rad][:self.RadiProc.lni]

                inb = self.BornProc.subproc[born][:self.BornProc.lni]
                fib = self.BornProc.subproc[born][self.BornProc.lni:]

                if inb==inr:
                    print rad,'and',born,'have the same initial state'
                    extra = 0

                else:
                    print rad,'and',born,'have different initial states'

def main(CONFIGFILE):
    
    AMPLITUDE = Amplitude(CONFIGFILE)

if __name__ == "__main__":
    main(sys.argv[1])