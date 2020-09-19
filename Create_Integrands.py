#!/usr/bin/env python

from Particle_Engine import *
from Utilities import *

class Amplitude:
    def __init__(self,CONFIGFILE):
        self.ReadConfig(CONFIGFILE)
        self.LoadConfig()
        self.VetoProcs()
        self.Build()

    #
    # Configfile reading and loading into the amplitude class
    #

    def ReadConfig(self,CONFIGFILE):
        CONFIG = {'INITIAL STATE': [0],\
                  'FINAL STATE'  : [0],\
                  'RADIATIVE'    : [0],\
                  'MODEL FILE'    : ['StandardModel'],\
                  'PATH'         : [os.getcwd()],\
                  'NLOX PATH'    : [0] ,\
                  'SEED TEMP'    : ['tpl/NLOX_seed.tpl'] ,\
                  'VERBOSE'      : [0] ,\
                  'STAGE'        : [0] }
        
        REQUIRED = set(['INITIAL STATE','FINAL STATE','NLOX PATH','RADIATIVE']) 

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
            LoadedModel = __import__(self.config['MODEL FILE'][0])
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
                self.finrad.append(self.ParticleContent[particle])
                if self.ParticleContent[particle].nam in self.Model.Fundamentals:
                    self.fin.append(self.ParticleContent[particle])
            except:
                print "\33[31mError\33[0m: Particle class",particle,"not defined."
                sys.exit()
        
        try:    
            self.finrad.append(self.ParticleContent[self.config['RADIATIVE'][0]])
        except:
            print "\33[31mError\33[0m: Particle class",self.config['RADIATIVE'][0],"not defined"
            sys.exit()

        self.RadiProc = Process(self.ini,self.finrad)
        self.path = self.config['PATH'][0]+"/"+str(self.RadiProc)

        self.seed_tpl = self.config['SEED TEMP'][0]
        self.stage = int(self.config['STAGE'][0])

        if self.verbose:
            print 'Configuartion loaded succesfuly'
            for conf in self.config:
                print conf,'=',self.config[conf]

    def VetoProcs(self):
        ToDelete = []
        for subproc in self.RadiProc.subproc:
            CC = []
            self.Model.ConnectedComponents(self.RadiProc.subproc[subproc],CC)
            if len(CC) != 1:
                ToDelete.append(subproc)
        for proc in ToDelete:
            del self.RadiProc.subproc[proc]

    def Build(self):

        self.TplDir = 'tpl/'
        self.IntDir = 'src/'
        self.SrcDir = self.path
        self.MatDir = self.SrcDir+'/Matrix_Elements'

        if self.stage:
            MakeDir(self.SrcDir)
            self.BuildDipoleTree()
            self.BuildSeeds()
            self.stage -= 1

        if self.stage: 
            # self.BuildMatrixElements()
            # self.BuildProcess()
            self.stage -= 1

        if self.stage:
            self.WriteDipoles()

    #
    # Dipole Tree creation from the Radiative process
    #

    def WriteDipoles(self):
        for Radiative in self.RadiProc.subproc:

            MakeDir(self.SrcDir+'/'+Radiative)

            TAB3 = '    '
            TAB6 = TAB3+TAB3

            DICT={'SubProcHeader'    : '' ,\
                  'SubProcConst'     : '' ,\
                  'NCC'              : '' ,\
                  'SubProcName'      : '' ,\
                  'SubProcMat'       : '' ,\
                  'SubProcSub'       : '' ,\
                  'SubProcPlu'       : '//empty for now' ,\
                  'SubProcEnd'       : '//empty for now'}
            
            DICT['SubProcHeader'] = Radiative.upper()
            DICT['SubProcName'] = Radiative+'_Dipoles'
            DICT['SubProcConst'] += TAB3+'std::unordered_map<std::string,int> CPMap;\n'
            DICT['SubProcConst'] += TAB3+'Proc = new Process();\n' 
            
            RadiCPS = self.Model.GetCPS(self.RadiProc.subproc[Radiative])
            Next = self.RadiProc.nex
        
            ##
            ##   Split all dipoles into EWK and QCD
            ## 

            EWKDipoles = []
            QCDDipoles = []
            for dipole in self.DipoleTree[Radiative]:
                if dipole['DIPTYP']=='EWK':
                    EWKDipoles.append(dipole)
                elif dipole['DIPTYP']=='QCD':
                    QCDDipoles.append(dipole)

            DIPDIC = {'II':'ab','IF':'ai','FI':'ia','FF':'ij'}

            for CP in RadiCPS:

                DICT['SubProcSub'] += TAB3+'if (cp =="'+CP+'"){\n'
                DICT['SubProcSub'] += TAB6+'double born[3];\n'
                DICT['SubProcSub'] += TAB6+'double aux = 0.;\n'
                DICT['SubProcSub'] += TAB6+'double DipFac = 1;\n'
                DICT['SubProcSub'] += TAB6+'i = Proc->AmpMap.at("'+Radiative+'");\n'
                DICT['SubProcSub'] += TAB6+'Proc->evaluate_alpha(i,"tree_tree","'+CP+'",p,'+str(Next)+',mu,radiative,acc);\n'
                
                ##
                ##   EWK Dipoles 
                ##

                for Emitter in EWKDipoles:
                    for Spectator in EWKDipoles:
                        if Emitter == Spectator:
                            continue
                        CPB = 'as'+str(RadiCPS[CP]['as'])+'ae'+str(RadiCPS[CP]['ae']-1)
                        try:
                            K = [ i for i in Emitter['IJ'] if i in Spectator['IJ']][0]
                            I = [ i for i in Emitter['IJ'] if i != K ][0]
                            J = [ i for i in Spectator['IJ'] if i != K ][0]

                            SigmaI = ( 1 if (( I<self.RadiProc.lni and self.RadiProc.subproc[Radiative][I].pid>0)\
                                             or  ( I>self.RadiProc.lni and self.RadiProc.subproc[Radiative][I].pid<0)) else -1 )
                            SigmaJ = ( 1 if (( J<self.RadiProc.lni and self.RadiProc.subproc[Radiative][J].pid>0)\
                                             or  ( J>self.RadiProc.lni and self.RadiProc.subproc[Radiative][J].pid<0)) else -1 )
                            
                            QI = self.RadiProc.subproc[Radiative][I].sym['Charge']
                            QJ = self.RadiProc.subproc[Radiative][J].sym['Charge']

                            ##
                            ## There is a factor of 3 between each charge in StandardModel.py vs RealWorld!
                            ## Hence the (1/9) for squared charges. 
                            ## StandardModel uses Charges*3 to enforce charge conservation using integers.
                            ##

                            ## 
                            ##  Subtracted Function 
                            ## 

                            PREFIX = Emitter['SUBTYP']+Spectator['SUBTYP']
                            DIPFUN = DIPDIC[PREFIX]

                            DICT['SubProcSub'] += TAB6+'DipFac *= '+str(SigmaI*SigmaJ*QI*QJ)+'./9;\n'
                            DICT['SubProcSub'] += TAB6+'i = Proc->AmpMap.at("'+Emitter['BORNTAG']+'");\n'
                            DICT['SubProcSub'] += TAB6+'Build_'+PREFIX+'_Momenta('+str(Next-1)+\
                                                       ',p,p_tilde,'+str(I)+','+str(J)+','+str(K)+');\n'
                            DICT['SubProcSub'] += TAB6+'Proc->evaluate_alpha(i,"tree_tree","'+CPB+'",p_tilde,'+str(Next-1)+',mu,born,acc);\n'
                            DICT['SubProcSub'] += TAB6+'aux -= ((Proc->pc.alpha_e)/2*M_PI)*g_'+DIPFUN+'_fermion(p['+str(I)+'],p['+str(J)+'],p['+str(K)+'])*born[2];\n'
                        except:
                            DICT['SubProcSub'] += str(Emitter)+'\n'+str(Spectator)+'\n'
                
                ##
                ##  QCD Dipoles
                ##

                for Emitter in QCDDipoles:
                    for Spectator in QCDDipoles:
                        if Emitter == Spectator:
                            continue
                        CPB = 'as'+str(RadiCPS[CP]['as']-1)+'ae'+str(RadiCPS[CP]['ae'])
                        try:
                            K = [ i for i in Emitter['IJ'] if i in Spectator['IJ']][0]
                            I = [ i for i in Emitter['IJ'] if i != K ][0]
                            J = [ i for i in Spectator['IJ'] if i != K ][0]
                            
                            ## 
                            ##  Subtracted Function 
                            ## 

                            DICT['SubProcSub'] += TAB6+'i = Proc->AmpMap.at("'+Emitter['BORNTAG']+'");\n'
                            DICT['SubProcSub'] += TAB6+'QCD_Build_'+Emitter['SUBTYP']+Spectator['SUBTYP']+'_Momenta('+str(Next-1)+\
                                                       ',p,p_tilde,'+str(I)+','+str(J)+','+str(K)+');\n'
                            DICT['SubProcSub'] += TAB6+'Proc->evaluate_alpha(i,"tree_tree","'+CPB+'",p_tilde,'+str(Next-1)+',mu,born,acc);\n'
                            DICT['SubProcSub'] += TAB6+'aux -= ((Proc->pc.alpha_e)/2*M_PI)*born[2];\n\n'
                        except:
                            DICT['SubProcSub'] += str(Emitter)+'\n'+str(Spectator)+'\n'

                DICT['SubProcSub'] += TAB3+'}\n\n'
            
            SubFunctionH = seek_and_destroy(self.TplDir+'/Dipole_FunctionsH.tpl',DICT)
            SubFunctionC = seek_and_destroy(self.TplDir+'/Dipole_FunctionsC.tpl',DICT)
            WriteFile(self.SrcDir+'/'+Radiative+'_Function.cpp',SubFunctionC)
            WriteFile(self.SrcDir+'/'+Radiative+'_Function.h',SubFunctionH)
    
    def BuildDipoleTree(self):
        # This method is in charged of fetching all 
        # dipole combinations and collect the unique
        # borns needed.
        self.DipoleTree = {}
        self.Borns = {}
        for radi in self.RadiProc.subproc:
            aux = {}
            self.GetDipoles(radi,aux)
            self.DipoleTree[radi] = aux[radi]
            for born in self.DipoleTree[radi]:
                self.Borns[born['BORNTAG']] = born['SORTEDPARS']

    def GetDipoles(self,SUBPROCESS,LINKED):
        # The rules to build the allowed dipoles are hard-coded Standard Model
        # This needs to be generalized or autogenerated from the model file

        RADIATIVE = self.RadiProc.subproc[SUBPROCESS]
        LINKED[SUBPROCESS] = []
        for id1 in range(len(RADIATIVE)):
            for id2 in range(len(RADIATIVE)):
                if id2 <= id1 or (id1 < self.RadiProc.lni and id2 < self.RadiProc.lni):
                    continue
                for newparticle in self.Model.Fundamentals:
                    dipole = []
                    SUBTYP = ''
                    if id1 < self.RadiProc.lni:
                        dipole = [RADIATIVE[id1],self.Model.Fundamentals[newparticle],RADIATIVE[id2]]
                        SUBTYP = 'I'
                    else:
                        dipole = [self.Model.Fundamentals[newparticle],RADIATIVE[id1],RADIATIVE[id2]]
                        SUBTYP = 'F'
                    dipole_nam = dipole[0].nam+'('+str(id1)+') -> '+dipole[1].nam+'('+str(id1)+') + '+dipole[2].nam+'('+str(id2)+')'
                    if self.Model.ParticleContent['h'] in dipole or self.Model.ParticleContent['Z'] in dipole:
                        continue
                    if self.Model.ParticleContent['g'] in dipole and self.Model.ParticleContent['A'] in dipole:
                        continue
                    if (self.Model.ParticleContent['Wp'] in dipole or self.Model.ParticleContent['Wm'] in dipole)\
                    and self.Model.ParticleContent['A'] not in dipole:
                        continue
                    if dipole[0] == dipole[1] == dipole[2] == self.Model.ParticleContent['A']:
                        continue
                    
                    syms = {}
                    first = True
                    for par in dipole:
                        sign = -1
                        if first:
                            sign = +1
                            first = False
                        for sym in par.sym:
                            if sym in syms:
                                syms[sym] += sign*par.sym[sym]
                            else:
                                syms[sym] = sign*par.sym[sym]
                    
                    reject = False
                    for sym in syms:
                        if syms[sym] != 0:
                            reject = True

                    if not reject:
                        AUX = []
                        for idc in range(len(RADIATIVE)):
                            if idc == id1:
                                AUX.append(self.Model.Fundamentals[newparticle])
                                continue
                            if idc == id2:
                                continue
                            AUX.append(RADIATIVE[idc])
                    
                        Skip = False
                        for finfun in self.fin:
                            if finfun not in AUX[self.RadiProc.lni:]:
                                Skip = True
                        
                        if not Skip:
                            def MyF(p):
                                return -p.pid
                            
                            INI = AUX[:self.RadiProc.lni]
                            FIN = AUX[self.RadiProc.lni:]
                            INI.sort(key=MyF)
                            FIN.sort(key=MyF)
                            SAUX = INI + FIN

                            CC = []
                            self.Model.ConnectedComponents(AUX,CC)
                            TYP = ''
                            if self.Model.ParticleContent['g'] in dipole:
                                # print '|M('+SUBPROCESS+')|^2 - QCD Dipole('+dipole_nam+')'+self.RadiProc.BuildString(SAUX)
                                TYP = 'QCD'
                            elif self.Model.ParticleContent['A'] in dipole:
                                # print '|M('+SUBPROCESS+')|^2 - EWK Dipole('+dipole_nam+')'+self.RadiProc.BuildString(SAUX)
                                TYP = 'EWK'
                            else:
                                print '\33[31mError\33[0m: Requested Dipole not found',dipole_nam
                                sys.exit()

                            if len(CC)==1:
                                LINKED[SUBPROCESS].append({'BORNTAG':self.RadiProc.BuildString(SAUX),'UNSORTEDPARS':AUX,'SORTEDPARS':SAUX,'DIPTYP':TYP,'IJ':[id1,id2],'SUBTYP':SUBTYP})
    
    #
    #  NLOX: Seeeds, Matrix elements and Process class
    #

    def BuildSeeds(self):

        for subproc in self.RadiProc.subproc:
            DICT = {'Initial State'   : '' ,\
                    'Final State'     : '' ,\
                    'Tree CP'         : 'bornAlphaQCD = ' ,\
                    'Loop CP'         : 'virtAlphaQCD = ' ,\
                    'enableTerms'      : '4' ,\
                    'Path'            : '' }
            Ini = [ s.nam for s in self.RadiProc.subproc[subproc][:self.RadiProc.lni]]
            Fin = [ s.nam for s in self.RadiProc.subproc[subproc][self.RadiProc.lni:]]

            EWKc = 0
            QCDc = 0
            for particle in self.RadiProc.subproc[subproc]:
                if particle in self.Model.EWKPars:
                    EWKc += 1
                if particle in self.Model.QCDPars:
                    QCDc += 1
            MaxQCD = QCDc - 2
            MinQCD = self.RadiProc.nex - EWKc 

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

            DICT['Tree CP'] += CSS(MaxQCD,MinQCD)
            DICT['Loop CP'] = '##virtAlphaQCD'
            DICT['Path'] +=  self.MatDir
            SubProcSeed = seek_and_destroy(self.seed_tpl,DICT)
            WriteFile(self.SrcDir+'/'+subproc+'.in',SubProcSeed)

        for born in self.Borns:
            DICT = {'Initial State'   : '' ,\
                    'Final State'     : '' ,\
                    'Tree CP'         : 'bornAlphaQCD = ' ,\
                    'Loop CP'         : 'virtAlphaQCD = ' ,\
                    'enableTerms'      : '7' ,\
                    'Path'            : '' }
            Ini = [ s.nam for s in self.Borns[born][:self.RadiProc.lni]]
            Fin = [ s.nam for s in self.Borns[born][self.RadiProc.lni:]]

            EWKc = 0
            QCDc = 0
            for particle in self.Borns[born]:
                if particle in self.Model.EWKPars:
                    EWKc += 1
                if particle in self.Model.QCDPars:
                    QCDc += 1
            MaxQCD = QCDc - 2
            MinQCD = len(self.Borns[born]) - EWKc

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
            
            DICT['Tree CP'] += CSS(MaxQCD,MinQCD)
            DICT['Loop CP'] += CSS(MaxQCD+1,MinQCD)
            DICT['Path'] +=  self.MatDir
            SubProcSeed = seek_and_destroy(self.seed_tpl,DICT)
            WriteFile(self.SrcDir+'/'+born+'.in',SubProcSeed)

    def BuildProcess(self):
        DICT = {    'Include Born'    : '' ,\
                     'Include Radi'   : '' ,\
                     'NSubProcesses'  : str(len(self.Borns)+len(self.RadiProc.subproc)) ,\
                     'Construct Born' : '\n' ,\
                     'Construct Radi' : '\n' ,\
                     'Amplitude Map'  : '\n'}
        TAB3 = '    '
        TAB6 = TAB3+TAB3
        TAB9 = TAB6+TAB3
        
        count = 0
        for subproc in self.Borns:
            DICT['Include Born'] += '#include "../'+subproc+'/code/Sub_'+subproc+'.h"\n'
            DICT['Construct Born'] += TAB9+'subproc['+str(count)+'] = '+'new Sub_'+subproc+'(pc);\n'
            DICT['Amplitude Map'] += TAB9+'AmpMap.insert({"'+subproc+'",'+str(count)+'});\n'
            count += 1

        for subproc in self.RadiProc.subproc:
            DICT['Include Radi'] += '#include "../'+subproc+'/code/Sub_'+subproc+'.h"\n'
            DICT['Construct Radi'] += TAB9+'subproc['+str(count)+'] = '+'new Sub_'+subproc+'(pc);\n'
            DICT['Amplitude Map'] += TAB9+'AmpMap.insert({"'+subproc+'",'+str(count)+'});\n'
            count += 1
        
        ProcessHeader = seek_and_destroy(self.TplDir+"nlox_process.tpl",DICT)
        WriteFile(self.MatDir+"/code/nlox_process.h",ProcessHeader)
        CopyFile(self.IntDir+'/nlox_olp.cc',self.MatDir+'/code/nlox_olp.cc')
        CopyFile(self.IntDir+'/nlox_olp.h',self.MatDir+'/code/nlox_olp.h')

    def BuildMatrixElements(self):
        Here = os.getcwd()
        os.chdir(self.SrcDir)
        
        count = 1
        tot = str(len(self.Borns))
        for subproc in self.Borns:
            if self.verbose:
                print 'Generating born and virtual process: ',subproc,'('+str(count)+'/'+tot+')'
            cmd = self.nlox+' '+subproc +'.in > '+subproc+'.log'
            os.system(cmd)
            cmd = 'mv '+subproc+'.log '+self.MatDir+'/'+subproc
            os.system(cmd)
            cmd = 'rm '+subproc+'.in '
            os.system(cmd)
            count += 1
        
        count = 1
        tot = str(len(self.RadiProc.subproc))
        for subproc in self.RadiProc.subproc:
            if self.verbose:
                print 'Generating real radiative process:',subproc,'('+str(count)+'/'+tot+')'
            cmd = self.nlox+' '+subproc +'.in > '+subproc+'.log'
            os.system(cmd)
            cmd = 'mv '+subproc+'.log '+self.MatDir+'/'+subproc
            os.system(cmd)
            cmd = 'rm '+subproc+'.in '
            os.system(cmd)
            count += 1

        os.chdir(Here)
    
    #
    #  Dipoles: Overhead Interface, Integrands and CUBA Targets 
    #
    

def main(CONFIGFILE):
    
    AMPLITUDE = Amplitude(CONFIGFILE)

if __name__ == "__main__":
    if len(sys.argv)<2:
        print "\33[31mError\33[0m: No input file specified"
    elif len(sys.argv)==2:
        main(sys.argv[1])
    else:
        print "\33[31mError\33[0m: Too many arguments"  
        
    