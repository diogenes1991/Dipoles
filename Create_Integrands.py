#!/usr/bin/env python

from Particle_Engine import *
from Utilities import *

class Amplitude:
    def __init__(self,CONFIGFILE):
        self.ReadConfig(CONFIGFILE)
        self.LoadConfig()
        self.VetoProcs()
        self.Build()

    ##
    ##  Configfile reading and loading into the amplitude class
    ##  We need to get rid of the extra wrapping the configs get,
    ##  all real configs are at self.config['Key'][0]
    ##

    def ReadConfig(self,CONFIGFILE):
        CONFIG = {'INITIAL STATE': [0],\
                  'FINAL STATE'  : [0],\
                  'RADIATIVE'    : [0],\
                  'MODEL FILE'    : ['StandardModel'],\
                  'PATH'         : [os.getcwd()],\
                  'NLOX PATH'    : [0] ,\
                  'SEED TEMP'    : ['tpl/NLOX_seed_tpl.in'] ,\
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
        if not os.path.exists(self.seed_tpl):
            print '\33[31mError\33[0m: NLOX input seed',self.seed_tpl,'not found'
            sys.exit()
        
        self.nlox = self.config['NLOX PATH'][0]+'/nlox.py'
        if not os.path.exists(self.nlox):
            print '\33[31mError\33[0m: NLOX executable',self.nlox,'not found'
            sys.exit()

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

        self.TplDir = 'tpl'
        self.IntDir = 'src'
        self.SrcDir = self.path
        self.MatDir = self.SrcDir+'/NLOX_Process'

        if self.stage:
            MakeDir(self.SrcDir)
            self.BuildDipoleTree()
            self.BuildSeeds()
            self.stage -= 1

        if self.stage: 
            self.BuildMatrixElements()
            self.BuildProcess()
            self.stage -= 1

        if self.stage:
            self.BuildDipoleTree()
            self.CopyInterface()
            self.BuildDipoleStructures()
            
    ##
    ##  Dipole Tree creation from the Radiative process
    ##

    def BuildDipoleTree(self):
        # This method is in charged of fetching all 
        # dipole combinations and collect the unique
        # borns needed.
        self.DipoleTree = {}
        self.Borns = {}
        self.Linked = {}
        for radi in self.RadiProc.subproc:
            self.Linked[radi] = {}
            aux = {}
            self.GetDipoles(radi,aux)
            self.DipoleTree[radi] = aux[radi]
            for born in self.DipoleTree[radi]:
                self.Borns[born['SBORNTAG']] = born['SPARS']
                self.Linked[radi][born['SBORNTAG']] = born['SPARS']

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
                                LINKED[SUBPROCESS].append({'SPARS':SAUX,'SBORNTAG':self.RadiProc.BuildString(SAUX),\
                                                           'UPARS':AUX,'UBORNTAG':self.RadiProc.BuildString(AUX),\
                                                           'DIPTYP':TYP,'IJ':[id1,id2],'SUBTYP':SUBTYP,\
                                                           'NAME':dipole_nam})
    
    ##
    ##  NLOX: Seeeds, Matrix elements and Process class
    ##

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
        DICT = {     'Include Born'    : '' ,\
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
        
        ## 
        ##  NLOX Interface  
        ##
        
        NLOX_PROCESS_H = seek_and_destroy(self.TplDir+"nlox_process_tpl.h",DICT)
        WriteFile(self.MatDir+"/code/nlox_process.h",NLOX_PROCESS_H)
        CopyFile(self.IntDir+'/nlox_olp.h',self.MatDir+'/code/nlox_olp.h')
        CopyFile(self.IntDir+'/nlox_olp.cc',self.MatDir+'/code/nlox_olp.cc')
        CopyFile(self.IntDir+'/nlox_olp_fortran.h',self.MatDir+'/code/nlox_olp_fortran.h')
        CopyFile(self.IntDir+'/nlox_olp_fortran.cc',self.MatDir+'/code/nlox_olp_fortran.cc')
        CopyFile(self.IntDir+'/nlox_fortran_interface.f90',self.MatDir+'/code/nlox_fortran_interface.f90')

        ##
        ##  NLOX Test Code
        ## 

        # CopyFile(self.IntDir+'/ftest_process.f90',self.MatDir+'/test_process.f90')
        # CopyFile(self.IntDir+'/test_process.cc',self.MatDir+'/test_process.cc')

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
    
    ##
    ##  Dipoles: Overhead Interface, Integrands and CUBA Targets 
    ##
    
    def CopyInterface(self):
        INTEGRANDFILES = ['Dipole_Structure.h','Dipole_Definitions.h','PSP_Generator.h']
        TOPLAYERFILES  = ['Test_Process.cpp','Test_Integrands.cpp']

        MakeDir(self.SrcDir+'/Integrands')

        for file in INTEGRANDFILES:
            CopyFile(self.IntDir+'/'+file,self.SrcDir+'/Integrands/'+file)
        
        for file in TOPLAYERFILES:
            CopyFile(self.IntDir+'/'+file,self.SrcDir+'/'+file)

        MAKEFILEDICT = {'NLOX PATH':self.config['NLOX PATH'][0],'GSL PATH':self.config['GSL PATH'][0]}
        MAKEFILE = seek_and_destroy(self.TplDir+'/makefile',MAKEFILEDICT)
        WriteFile(self.SrcDir+'/makefile',MAKEFILE)
        
    def BuildDipoleStructures(self):

        def ProcessConstMassName(particle):

            ##  This is a NLOX-specific function, it translates the names from StandardModel.py
            ##  into the ones in the NLOX's Processconst class
            
            massvalue = 'Proc->pc.m'
            if (particle in self.Model.quarks):
                if abs(particle.pid) == 2:
                    massvalue += particle.nam[0]+'p.real()'
                else:
                    massvalue += particle.nam[0]+'.real()'
            elif (particle in self.Model.leptons):
                if particle.sym['Charge']!=0:
                    massvalue = '0.0'
                else:
                    massvalue += particle.nam[:1]+'.real()'
            elif isinstance(particle,Boson):
                if (particle.nam=='g' or particle.nam == 'A'):
                    massvalue = '0.0'
                else:
                    massvalue = 'sqrt(Proc->pc.m'+particle.nam[0]+'2.real())'
            return massvalue

        ## 
        ##  This function has waaaay too many reponsabilities.... 
        ##

        TAB3 = '    '
        TAB6 = TAB3+TAB3
        TAB9 = TAB6+TAB3

        INTDICT = {'Include Integrands' : '' ,\
                   'nRadiative'         : str(len(self.RadiProc.subproc)) ,\
                   'Integrand Catalogue': ''}

        Count = 0
        for Radiative in self.RadiProc.subproc:

            ClassName = Radiative+'_Integrands'

            INTDICT['Include Integrands'] += '#include "Integrands/'+Radiative+'/'+ClassName+'.h"\n'
            INTDICT['Integrand Catalogue'] += TAB9+'Channels['+str(Count)+'] = new '+ClassName+'(Proc);\n'
            INTDICT['Integrand Catalogue'] += TAB9+'ChannelMap.insert({"'+Radiative+'",'+str(Count)+'});\n'

            Count += 1
            ChlDir = self.SrcDir+'/Integrands/'+Radiative
            MakeDir(ChlDir)

            DICT={'SubProcHeader'    : '' ,\
                  'SubProcConst'     : '' ,\
                  'SubProcName'      : '' ,\
                  'Next'             : str(self.RadiProc.nex) ,\
                  'SubProcMat'       : '' ,\
                  'SubProcSub'       : '' ,\
                  'SubProcPlu'       : '',\
                  'SubProcEnd'       : ''}
            
            DICT['SubProcHeader'] = Radiative.upper()
            DICT['SubProcName'] = ClassName

            ## 
            ##  Integrands constructor 
            ##

            count = 0 
            for particle in self.RadiProc.subproc[Radiative]:
                massvalue = ProcessConstMassName(particle)
                DICT['SubProcConst'] += TAB3+'Masses['+str(count)+'] = '+massvalue+';\n'
                count += 1
    
            DICT['SubProcConst'] += '\n'   
            DICT['SubProcConst'] += TAB3+'nBorn = '+str(len(self.Linked[Radiative]))+';\n'   
            DICT['SubProcConst'] += TAB3+'BornMasses = new double* [nBorn];\n' 
            DICT['SubProcConst'] += TAB3+'for(int i=0;i<nBorn;i++) BornMasses[i]= new double[Next-1];\n\n' 
             
            Bcount = 0
            for Born in self.Linked[Radiative]:
                count = 0
                DICT['SubProcConst'] += TAB3+'BornMap.insert({"'+Born+'",'+str(Bcount)+'});\n'
                for particle in self.Linked[Radiative][Born]:
                    DICT['SubProcConst'] += TAB3+'BornMasses['+str(Bcount)+']['+str(count)+']='+ProcessConstMassName(particle)+';\n'
                    count += 1
                Bcount += 1
            
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

            def BORNTAG(Dipole):
                return Dipole['SBORNTAG']
            EWKDipoles.sort(key=BORNTAG)
            QCDDipoles.sort(key=BORNTAG)

            DIPDIC = {'II':'ab','IF':'ai','FI':'ia','FF':'ij'}

            ## 
            ##  The DipoleTree has in it all the information we require:
            ##  The postprocessing of the Dipole tree is as follows:
            ##   - Fisrt we split the tree in two: QCD and EWK, mainly beca 
            ##  
            ##  

            for CP in RadiCPS:

                CPHEAD = TAB3+'if (cp =="'+CP+'"){\n'
                CACHEDTAG = ''

                if len(EWKDipoles):
                    CPHEAD += TAB6+'double born[3];\n'
                if len(QCDDipoles):
                    CPHEAD += TAB6+'double borncc['+str(1+(self.RadiProc.nex-1)*(self.RadiProc.nex-2)/2)+'];\n'
                    
                DICT['SubProcSub'] += CPHEAD
                DICT['SubProcSub'] += TAB6+'i = Proc->AmpMap.at("'+Radiative+'");\n'
                DICT['SubProcSub'] += TAB6+'Proc->evaluate_alpha(i,"tree_tree","'+CP+'",p,'+str(Next)+',mu,radiative,&acc);\n'
                DICT['SubProcSub'] += TAB6+'rval[2] += radiative[2];\n\n'

                DICT['SubProcPlu'] += CPHEAD

                DICT['SubProcEnd'] += CPHEAD
                DICT['SubProcEnd'] += TAB6+'double aux[3];\n'
                        
                ##
                ##   EWK Dipoles 
                ##

                for Emitter in EWKDipoles:
                    for Spectator in EWKDipoles:
                        if Emitter == Spectator:
                            continue
                        
                        ## This throws away any CP combination that has no Born counterpart
                        if RadiCPS[CP]['ae'] <= 0:
                            continue
                        CPB = 'as'+str(RadiCPS[CP]['as'])+'ae'+str(RadiCPS[CP]['ae']-1)
                        
                        K = []

                        try:
                            K = [ i for i in Emitter['IJ'] if i in Spectator['IJ']][0]
                        except:
                            # print Emitter['NAME'],Spectator['NAME'],'-> Skipped: Different Radiation'
                            continue 
                        
                        I = [ i for i in Emitter['IJ'] if i != K ][0]
                        J = [ i for i in Spectator['IJ'] if i != K ][0]

                        SigmaI = ( 1 if (( I<self.RadiProc.lni and self.RadiProc.subproc[Radiative][I].pid>0)\
                                         or  ( I>self.RadiProc.lni and self.RadiProc.subproc[Radiative][I].pid<0)) else -1 )
                        SigmaJ = ( 1 if (( J<self.RadiProc.lni and self.RadiProc.subproc[Radiative][J].pid>0)\
                                         or  ( J>self.RadiProc.lni and self.RadiProc.subproc[Radiative][J].pid<0)) else -1 )
                        
                        QI=0
                        QJ=0

                        try:                            
                            QI = self.RadiProc.subproc[Radiative][I].sym['Charge']
                            QJ = self.RadiProc.subproc[Radiative][J].sym['Charge']
                        except:
                            continue
                        
                        if CACHEDTAG != Emitter['SBORNTAG']:
                            AMPMAPLINE = TAB6+'i = Proc->AmpMap.at("'+Emitter['SBORNTAG']+'");\n'
                            BORNMAPLINE  = TAB6+'j = BornMap.at("'+Emitter['SBORNTAG']+'");\n'
                            BORNMAPLINE += TAB6+'SetInMom(j);\n'
                            BORNMAPLINE += TAB6+'SetFiMom(j,rand,&J);\n'
                            BORNMAPLINE += TAB6+'for(int k=0;k<=2;k++) p.push_back(BornMomenta[k]);\n'
                            CACHEDTAG = Emitter['SBORNTAG']
                        else:
                            AMPMAPLINE = ''
                            BORNMAPLINE = ''

                        ##
                        ## There is a factor of 3 between each charge in StandardModel.py vs RealWorld!
                        ## Hence the (1/9) for squared charges. 
                        ## StandardModel uses Charges*3 to enforce charge conservation using integers.
                        ##

                        AMPMAPLINE  += TAB6+'DipFac = '+str(SigmaI*SigmaJ*QI*QJ)+'.0/9*EWKFac;\n'
                        
                        ## 
                        ##  Subtracted Function 
                        ## 

                        PREFIX = Emitter['SUBTYP']+Spectator['SUBTYP']
                        DIPFUN = DIPDIC[PREFIX]

                        DICT['SubProcSub'] += AMPMAPLINE
                        DICT['SubProcSub'] += TAB6+'Build_'+PREFIX+'_Momenta('+str(Next-1)+\
                                                   ',p,p_tilde,'+str(I)+','+str(J)+','+str(K)+');\n'
                        DICT['SubProcSub'] += TAB6+'Proc->evaluate_alpha(i,"tree_tree","'+CPB+'",p_tilde,'+str(Next-1)+',mu,born,&acc);\n'
                        DICT['SubProcSub'] += TAB6+'*rval -= DipFac*g_'+DIPFUN+'_fermion(p['+str(I)+'],p['+str(J)+'],p['+str(K)+'],Masses['+str(I)+'],Masses['+str(J)+'])*born[2];\n\n'
                        
                        ##
                        ##  Plus Distribution
                        ##

                        DICT['SubProcPlu'] += AMPMAPLINE
                        DICT['SubProcPlu'] += TAB6+';\n'
                        DICT['SubProcPlu'] += TAB6+';\n'
                        
                        ##
                        ##  Endpoint Function
                        ##

                        DICT['SubProcEnd'] += AMPMAPLINE
                        DICT['SubProcEnd'] += BORNMAPLINE
                        DICT['SubProcEnd'] += TAB6+'Invariant = BornMasses[j]['+str(I)+']*BornMasses[j]['+str(I)+']'\
                                                              '+ BornMasses[j]['+str(J)+']*BornMasses[j]['+str(J)+']'\
                                                              +('+' if (PREFIX=='II'or PREFIX=='FF') else '-')+\
                                                              '2*(BornMomenta['+str(I)+']*BornMomenta['+str(J)+']);\n'
                        DICT['SubProcEnd'] += TAB6+'Proc->evaluate_alpha(i,"tree_tree","'+CPB+'",p,'+str(Next-1)+',mu,born,&acc);\n'
                        DICT['SubProcEnd'] += TAB6+'G_'+DIPFUN+'_fermion(Invariant,BornMasses[j]['+str(I)+'],BornMasses[j]['+str(J)+'],mu,aux);\n'
                        DICT['SubProcEnd'] += TAB6+'for(int k=0;k<=2;k++) rval[k] += DipFac*aux[k]*born[2];\n\n'
                        
                ##
                ##  QCD Dipoles
                ##

                for Emitter in QCDDipoles:
                    for Spectator in QCDDipoles:
                        if Emitter == Spectator:
                            continue

                        ## This throws away any CP combination that has no Born counterpart
                        if RadiCPS[CP]['as'] <= 0:
                            continue
                        CPB = 'as'+str(RadiCPS[CP]['as']-1)+'ae'+str(RadiCPS[CP]['ae'])

                        K = []
                        try:
                            K = [ i for i in Emitter['IJ'] if i in Spectator['IJ']][0]
                        except:
                            continue 

                        I = [ i for i in Emitter['IJ'] if i != K ][0]
                        J = [ i for i in Spectator['IJ'] if i != K ][0]

                        FIRST = min(I,J)
                        SECON = max(I,J)

                        if CACHEDTAG != Emitter['SBORNTAG']:
                            AMPMAPLINE = TAB6+'i = Proc->AmpMap.at("'+Emitter['SBORNTAG']+'");\n'
                            BORNMAPLINE  = TAB6+'j = BornMap.at("'+Emitter['SBORNTAG']+'");\n'
                            BORNMAPLINE += TAB6+'SetInMom(j);\n'
                            BORNMAPLINE += TAB6+'SetFiMom(j,rand,&J);\n'
                            BORNMAPLINE += TAB6+'for(int k=0;k<=2;k++) p.push_back(BornMomenta[k]);\n'
                            CACHEDTAG = Emitter['SBORNTAG']
                        else:
                            AMPMAPLINE = ''
                            BORNMAPLINE = ''

                        ## This unflattens a symmetric (non-repeating) matrix, which borncc is

                        WHICHCC = (self.RadiProc.nex-1)*FIRST + SECON - FIRST - FIRST*(FIRST-1)/2

                        PREFIX = Emitter['SUBTYP']+Spectator['SUBTYP']
                        DIPFUN = DIPDIC[PREFIX]

                        ## 
                        ##  Subtracted Function 
                        ## 

                        DICT['SubProcSub'] += AMPMAPLINE
                        DICT['SubProcSub'] += TAB6+'Build_'+PREFIX+'_Momenta('+str(Next-1)+\
                                                   ',p,p_tilde,'+str(I)+','+str(J)+','+str(K)+');\n'
                        DICT['SubProcSub'] += TAB6+'Proc->evaluate_alpha_cc(i,"tree_tree","'+CPB+'",p_tilde,'+str(Next-1)+',borncc,&acc);\n'
                        DICT['SubProcSub'] += TAB6+'*rval += QCDFac*g_'+DIPFUN+'_fermion(p['+str(I)+'],p['+str(J)+'],p['+str(K)+'],Masses['+str(I)+'],Masses['+str(J)+'])*borncc['+str(WHICHCC)+'];\n\n'
                        
                        ##
                        ##  Plus Distribution
                        ##

                        DICT['SubProcPlu'] += AMPMAPLINE

                        ## 
                        ##  Endpoint Function 
                        ##

                        DICT['SubProcEnd'] += AMPMAPLINE
                        DICT['SubProcEnd'] += BORNMAPLINE
                        DICT['SubProcEnd'] += TAB6+'Invariant = BornMasses[j]['+str(I)+']*BornMasses[j]['+str(I)+']'\
                                                              '+ BornMasses[j]['+str(J)+']*BornMasses[j]['+str(J)+']'\
                                                              +('+' if (PREFIX=='II'or PREFIX=='FF') else '-')+\
                                                              '2*(BornMomenta['+str(I)+']*BornMomenta['+str(J)+']);\n'
                        DICT['SubProcEnd'] += TAB6+'Proc->evaluate_alpha_cc(i,"tree_tree","'+CPB+'",p,'+str(Next-1)+',borncc,&acc);\n'
                        DICT['SubProcEnd'] += TAB6+'G_'+DIPFUN+'_fermion(Invariant,BornMasses[j]['+str(I)+'],BornMasses[j]['+str(J)+'],mu,aux);\n'
                        DICT['SubProcEnd'] += TAB6+'for(int k=0;k<=2;k++) rval[k] += QCDFac*aux[k]*borncc['+str(WHICHCC)+'];\n\n'
                        

                DICT['SubProcSub'] += TAB3+'}\n'
                DICT['SubProcPlu'] += TAB3+'}\n'
                DICT['SubProcEnd'] += TAB3+'}\n'
            
            SubFunctionH = seek_and_destroy(self.TplDir+'/Dipole_Functions_tpl.h',DICT)
            SubFunctionC = seek_and_destroy(self.TplDir+'/Dipole_Functions_tpl.cpp',DICT)
            WriteFile(ChlDir+'/'+Radiative+'_Integrands.cpp',SubFunctionC)
            WriteFile(ChlDir+'/'+Radiative+'_Integrands.h',SubFunctionH)
        
        IntegrandClass = seek_and_destroy(self.TplDir+'/Integrands_tpl.h',INTDICT)
        WriteFile(self.SrcDir+'/Integrands.h',IntegrandClass)
        
def main(CONFIGFILE):
    
    AMPLITUDE = Amplitude(CONFIGFILE)

if __name__ == "__main__":
    if len(sys.argv)<2:
        print "\33[31mError\33[0m: No input file specified"
    elif len(sys.argv)==2:
        main(sys.argv[1])
    else:
        print "\33[31mError\33[0m: Too many arguments"  
        
    