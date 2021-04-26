#!/usr/bin/env python

from time import localtime,strftime,mktime
from Dipoles import *
from OLP import *

class Madisqe:
    def __init__(self,CONFIGFILE):
        
        StartTime = localtime()

        LogFile = open('Madisqe.log','w')
        cmd = "git log --pretty=format:'commit %H%nDate:   %ad' -n 1 > Madisqe.log"
        os.system(cmd)
        LogFile.close()
        
        self.ReadConfig(CONFIGFILE)
        
        print '     Code Generation Started at '+strftime('%d %b %Y %H:%M:%S %Z',StartTime)
        print '                                                                   '
        print '###################################################################'
        print '                                                                   '
        print '     Sorting the necesary channels                                 '
        print '                                                                   '
        
        ##
        ##  Both the Born and Radiative processes are classified by
        ##  Processes which are the same up to sorting are not kept 
        ##  this must be all summed back by the c++ interface
        ##  
        ##  The veto here is on all processes that are not possible at tree level 
        ##  The issue is that if we allow these we can only compute them @ LO
        ##

        self.LoadConfig()
        self.VetoProcs()
        self.CreateDirectories()
        self.WriteModel()

        print '                                                                   '
        print '     Sorting completed at '+strftime('%d %b %Y %H:%M:%S %Z',localtime())
        print '     Total elapsed time: ' + str(mktime(localtime())-mktime(StartTime)) + ' s'
        print '                                                                   '        
        print '###################################################################'
        print '                                                                   '
        print '     The Born and Virtuel include the channels('+str(len(self.BornProc.subproc))+'):'
        print '                                                                   '
        
        for Ch in self.BornProc.subproc:
            print '     '+Ch
        
        print '                                                                   '
        print '###################################################################'
        print '                                                                   '
        print '     The Radiative includes the channels ('+str(len(self.RadiProc.subproc))+'):'
        print '                                                                   '
        
        for Ch in self.RadiProc.subproc:
            print '     '+Ch

        print '                                                                   '
        print '###################################################################'
        print '                                                                   '
        print '     Building Real Integrands                                      '
        print '                                                                   '
        
        self.BuildDipoleTree()
        self.BuildRealIntegrands()
        self.Build_Dipole_Structures()
        self.BuildRealIntegrands_New()

        print '                                                                   '
        print '     Dipole fetching completed at '+strftime('%d %b %Y %H:%M:%S %Z',localtime())
        print '     Total elapsed time: ' + str(mktime(localtime())-mktime(StartTime)) + ' s'
        print '                                                                   '
        print '###################################################################'
        print '                                                                   '
        print '     Generating matrix elements and Virtual Integrands             '
        print '                                                                   '
        
        self.BuildOLP(False)
        self.BuildVirtualIntegrands()

        print '                                                                   '
        print '     Generation of matrix elements completed at '+strftime('%d %b %Y %H:%M:%S %Z',localtime())
        print '     Total elapsed time: ' + str(mktime(localtime())-mktime(StartTime)) + ' s'
        print '                                                                   '
        print '###################################################################'
        print '                                                                   '
        print '     Building Madisqe Interface                                    '
        
        self.BuildInterface()
        
        print '                                                                   '
        print '     Interface building completed at '+strftime('%d %b %Y %H:%M:%S %Z',localtime())
        print '     Total elapsed time: ' + str(mktime(localtime())-mktime(StartTime)) + ' s'
        print '                                                                   '
        print '###################################################################'
        print '                                                                   '
        print '     Code generation finished at '+strftime('%d %b %Y %H:%M:%S %Z',localtime())
        print '     Total elapsed time: ' + str(mktime(localtime())-mktime(StartTime)) + ' s'
        print '                                                                   '
        print '###################################################################'
        
        ## 
        ##  TO DO: At the moment the code is a bit unstable 
        ##  We need to stabilize and clean as much as possible 
        ##  The next task will be to centralize masses into the 
        ##

        # for DipoleStrucure in self.DipoleStructures:
            # print self.DipoleStructures[DipoleStrucure]
            # self.DipoleStructures[DipoleStrucure].Show()

    ##
    ##  Configfile reading and loading into the Madisqe class
    ##

    def LoadPaths(self):
        PATHS = { 'GSL PATH'    : [''], \
                  'CUBA PATH'   : [''], \
                  'LHAPDF PATH' : [''], }

        PATHSFILE=os.getcwd()+'/util/utils.path'
        COMMCHAR = '#'
        REQUIRED = set([])

        f = open(PATHSFILE,'r')
        configfile= []
        for l in f:
            if l.find(COMMCHAR)==0 or len(l)==0 or l=='\n':
                continue
            configfile.append([False,l])

        for configuration in PATHS.keys():
            warn = True
            for line in configfile:
                if configuration in line[1]:
                    thisconfig = line[1].strip().split('=')
                    PATHS[thisconfig[0].strip()] = (thisconfig[1].strip()).split(',')[0]
                    warn = False
                    line[0] = True
            if warn:
                if configuration not in REQUIRED:
                    print "\33[95mWarning\33[0m: Configuartion",configuration,"not found in configuration file, using default:",PATHS[configuration][0]
                else:
                    print "\33[31mError\33[0m: Configuration",configuration,"not found in configuration, this is a required configuration."
                    sys.exit()
        
        for line in configfile:
            if not line[0]:
                print "\33[95mWarning\33[0m: Configuartion",(((line[1].strip()).split('='))[0]).strip(),"is not a valid configuration setting, please check the syntax"
        f.close()
        self.paths = PATHS

    def ReadConfig(self,CONFIGFILE):
        CONFIG = {'INITIAL STATE': [0],\
                  'FINAL STATE'  : [0],\
                  'RADIATIVE'    : [0],\
                  'MODEL FILE'   : ['StandardModel'],\
                  'PATH'         : [os.getcwd()],\
                  'NLOX PATH'    : [0] ,\
                  'RECOLA PATH'  : [0] ,\
                  'SEED TEMP'    : ['tpl/NLOX_seed_tpl.in'] ,\
                  'VERBOSE'      : [0] ,\
                  'STAGE'        : [0] }
        
        REQUIRED = set(['INITIAL STATE','FINAL STATE','NLOX PATH','RECOLA PATH','RADIATIVE']) 

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
                
        self.Model = LoadedModel.WorkingModel
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
                self.fin.append(self.ParticleContent[particle])
            except:
                print "\33[31mError\33[0m: Particle class",particle,"not defined."
                sys.exit()
        
        try:    
            self.finrad.append(self.ParticleContent[self.config['RADIATIVE'][0]])
        except:
            print "\33[31mError\33[0m: Particle class",self.config['RADIATIVE'][0],"not defined"
            sys.exit()

        self.RadiProc = Process(self.ini,self.finrad,self.Model)
        self.BornProc = Process(self.ini,self.fin,self.Model)
        self.path = self.config['PATH'][0]+"/"+str(self.RadiProc)

        self.seed_tpl = self.config['SEED TEMP'][0]
        if not os.path.exists(self.seed_tpl):
            print '\33[31mError\33[0m: NLOX input seed',self.seed_tpl,'not found'
            sys.exit()
        
        self.nlox = self.config['NLOX PATH'][0]+'/nlox.py'
        if not os.path.exists(self.nlox):
            print '\33[31mError\33[0m: NLOX executable',self.nlox,'not found'
            sys.exit()

        self.recola = self.config['RECOLA PATH'][0]
        if not os.path.exists(self.recola+'/recola2-2.2.2/librecola.so'):
            print '\33[31mError\33[0m: RECOLA library: librecola.so not found at',self.recola+'/recola2-2.2.2'
            sys.exit()

        self.stage = int(self.config['STAGE'][0])

        if self.verbose:
            print 'Configuartion loaded succesfuly'
            for conf in self.config:
                print conf,'=',(self.config[conf] if len(self.config[conf]) > 1 else self.config[conf][0])

        self.TplDir = 'tpl'
        self.IntDir = 'src'
        self.Inputs = 'input'
        self.SrcDir = self.path

    def VetoProcs(self):
        ToDelete = []
        for subproc in self.RadiProc.subproc:
            CC = []
            self.Model.ConnectedComponents(self.RadiProc.subproc[subproc],CC)
            if len(CC) != 1:
                ToDelete.append(subproc)
        for proc in ToDelete:
            del self.RadiProc.subproc[proc]
        
    def WriteModel(self):
        MOD_DICT = {'Build Model' : '' ,\
                    'NParticles'  : '' ,\
                    }
                    
        TAB4 = '    '
        TAB8 = TAB4+'    '

        MOD_DICT['NParticles'] += str(len(self.Model.Fundamentals))
        
        for Particle in self.Model.Fundamentals:
            P = self.Model.Fundamentals[Particle]
            Charge = '0'
            if 'Charge' in P.sym:
                ##
                ## The charges in SandardModel have a relative factor of 3 with respect to the real world
                ## this is to enforce charge conservation using integers
                ##
                Charge = str(P.sym['Charge'])+'/3'
            MOD_DICT['Build Model'] += TAB8+'Particle '+P.nam+' = Particle("'+P.nam+'",0.0,0.0,'+Charge+','+str(P.pid)+');\n' 
    
            
        MODEL = Template(self.TplDir+'/Model_tpl.h',self.SrcDir+'/Code/Model.h',MOD_DICT)
        MODEL.Write()

    ##
    ##  Dipole Tree creation from the Radiative process
    ##  BuildDipoleTree, GetDipoles and Build_Dipole_Structures_OLD are deprecated...
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
                    
                    partyps = ''
                    syms = {}
                    first = True
                    for par in dipole:
                        partyps += ('b' if isinstance(par,Boson) else 'f')
                        sign = -1
                        if first:
                            sign = +1
                            first = False
                        for sym in par.sym:
                            if sym in syms:
                                syms[sym] += sign*par.sym[sym]
                            else:
                                syms[sym] = sign*par.sym[sym]
                    # if partyps == 'fbf':
                    #     partyps = 'ffb'

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
                                                           'NAME':dipole_nam,'DIPOLE':dipole,'PARTYPS':partyps})

    def Build_Dipole_Structures_OLD(self,Verbose=False):

        # We first fetch the allowed Underlying Born tags allowed by the 
        # particle content specification in the input file

        ALLOWEDTAGS={}
        for Born in self.BornProc.subproc:
            TAG=""
            FINALS=self.BornProc.subproc[Born][self.RadiProc.lni:]
            for Particle in FINALS:
                TAG+=Particle.nam
            ALLOWEDTAGS[TAG]=FINALS

        # Now we proceed to build the Dipole structures

        self.DipoleStructures = {}
        for Radiative in self.RadiProc.subproc:
            self.DipoleStructures[Radiative] = []

            if Verbose:
                print "Building Dipoles for",Radiative

            Dipoles = []
            for Ri in range(self.RadiProc.lni,len(self.RadiProc.subproc[Radiative])): # Radiation must be a final particle
                Radiation = self.RadiProc.subproc[Radiative][Ri]


                if Verbose:
                    print "  Evaluating the possibility of considering "+Radiation.nam+"("+str(Ri)+") as radiation"

                # Up to this point we have selected a particle to be our Radiation
                # now we need to loop over the particle content of the Model to 
                # and build all possible radiation channels
                
                for Ei in range(len(self.RadiProc.subproc[Radiative])):
                    Emitter = self.RadiProc.subproc[Radiative][Ei]


                    if Verbose:
                        print "    Considering emission off of "+Emitter.nam+"("+str(Ei)+")"
                    
                    if Ei >= Ri: # Only consider pairs of emitter-radiation
                        continue

                    if not self.Model.CanMakeVertex(Emitter,Radiation):
                        continue

                    # We have to scenarios either the Emitter is initial or is final.
                    # The conservation of Model Charges needs to be applied differently in 
                    # each case:
                    #           - Initial Emitter: Emitter    -> REParticle + Radiation
                    #           - Final   Emitter: REParticle -> Emitter    + Radiation

                    for NewPar in self.Model.Fundamentals: # Loop over the Mdoel's Fundamentals
                        REParticle = self.Model.Fundamentals[NewPar]
                        
                        if REParticle.mas and Emitter.mas and Radiation.mas: # Skip dipoles with 3 massive particles
                            continue
                        
                        if Ei < self.RadiProc.lni:
                            DecayChannel = [Emitter,REParticle,Radiation]
                        else:
                            DecayChannel = [REParticle,Emitter,Radiation]

                        # We need to hand-skip A->AA , A->ZA, A->ZZ, this needs to be better addressed by self.Model.CanMakeVertex 
                        if DecayChannel[0] == self.Model.ParticleContent["A"] and DecayChannel[1] == self.Model.ParticleContent["A"] and DecayChannel[1] == self.Model.ParticleContent["A"]:
                            continue
                        if self.Model.ParticleContent["Z"] in DecayChannel and self.Model.ParticleContent["A"] in DecayChannel:
                            continue
                        
                        # This next block uses particle engine to veto using charge conservation

                        Symmetries = {}
                        First = True
                        for Particle in DecayChannel:
                            Sign = +1
                            if First:
                                Sign = -1
                                First = False
                            for Symmetry in Particle.sym:
                                if Symmetry in Symmetries:
                                    Symmetries[Symmetry] += Sign*Particle.sym[Symmetry]
                                else:
                                    Symmetries[Symmetry]  = Sign*Particle.sym[Symmetry]
                        
                        Reject = False
                        for Symmetry in Symmetries:
                            if Symmetries[Symmetry] != 0:
                                Reject = True
                        if Reject:
                            continue
                        
                        if Verbose:
                            print "      Cosidering Radiation+Emitter: "+REParticle.nam

                        # Now that we have built and veto the possible channels we 
                        # proceed to build the dependency tree for this radiative

                        UnderlyingBorn = []
                        for Index in range(len(self.RadiProc.subproc[Radiative])):
                            BornParticle = self.RadiProc.subproc[Radiative][Index]
                            UnderlyingBorn.append(BornParticle)
                        UnderlyingBorn[Ei] = REParticle
                        del UnderlyingBorn[Ri]

                        # Now we check if any of the resulting trees is disconected,
                        # if so this Decay channel is disallowed

                        def IsPossible(LIST):
                            CC=[]
                            self.Model.ConnectedComponents(LIST,CC)
                            if len(CC)==1:
                                return True
                            else:
                                return False

                        if not IsPossible(DecayChannel) or not IsPossible(UnderlyingBorn):
                            continue

                        # Now we require that the Underlying born configuration
                        # exists in the expansion of the Born tags

                        def MyF(p):
                                return -p.pid
                        FIN = UnderlyingBorn[self.RadiProc.lni:]
                        FIN.sort(key=MyF)
                        TAG=""
                        for Particle in FIN:
                            TAG+=Particle.nam

                        if TAG not in ALLOWEDTAGS:
                            continue

                        if REParticle in self.Model.QCDPars:
                            Dipoles.append(QCDDipole(DecayChannel,UnderlyingBorn))

                        if REParticle in self.Model.EWKPars:
                            Dipoles.append(EWKDipole(DecayChannel,UnderlyingBorn))

            self.DipoleStructures[Radiative]=DipoleStructure(Radiative,Dipoles)       
    
    def Build_Dipole_Structures(self,Verbose=False):

        ##
        ## We first fetch the allowed Underlying Born tags allowed by the 
        ## particle content specification in the input file
        ##

        ALLOWEDTAGS={}
        for Born in self.BornProc.subproc:
            TAG=""
            FINALS=self.BornProc.subproc[Born][self.RadiProc.lni:]
            for Particle in FINALS:
                TAG+=Particle.nam
            ALLOWEDTAGS[TAG]=FINALS

        ## 
        ## Now we proceed to build the Dipole structures
        ##

        self.DipoleStructures = {}
        for Radiative in self.RadiProc.Channels:
            RadiativeChannel = self.RadiProc.Channels[Radiative]
            self.DipoleStructures[Radiative] = []

            if Verbose:
                print "Building Dipoles for",Radiative

            Dipoles = []
            for Ri in range(self.RadiProc.lni,len(RadiativeChannel.Particles)): # Radiation must be a final particle[Ei,Ri],
                Radiation = RadiativeChannel.Particles[Ri]

                if Verbose:
                    print "  Evaluating the possibility of considering "+Radiation.nam+"("+str(Ri)+") as radiation"

                ##
                ## Up to this point we have selected a particle to be our Radiation
                ## now we need to loop over the particle content of the Model
                ## and build all possible radiation channels
                ##
                
                for Ei in range(len(RadiativeChannel.Particles)):
                    Emitter = RadiativeChannel.Particles[Ei]


                    if Verbose:
                        print "    Considering emission off of "+Emitter.nam+"("+str(Ei)+")"
                    
                    if Ei >= Ri: # Only consider pairs of emitter-radiation
                        continue

                    if not self.Model.CanMakeVertex(Emitter,Radiation):
                        continue

                    ##
                    ## We have to scenarios either the Emitter is initial or is final.
                    ## The conservation of Model Charges needs to be applied differently in 
                    ## each case:
                    ##           - Initial Emitter: Emitter    -> REParticle + Radiation
                    ##           - Final   Emitter: REParticle -> Emitter    + Radiation
                    ##

                    for NewPar in self.Model.Fundamentals: # Loop over the Mdoel's Fundamentals
                        REParticle = self.Model.Fundamentals[NewPar]
                        
                        if REParticle.mas and Emitter.mas and Radiation.mas: # Skip dipoles with 3 massive particles
                            continue
                        
                        Type = ''
                        if Ei < self.RadiProc.lni: # If the radiation is off the initial state only a massless RE is allowed
                            Type = 'ia'
                            ParentIndices = [Ei,Ri]
                            DecayChannel = Channel([Emitter],[REParticle,Radiation],self.Model)
                        else:
                            Type = 'ij'
                            ParentIndices = [Ei,Ri]
                            DecayChannel = Channel([REParticle],[Emitter,Radiation],self.Model)
                        ##
                        ## We need to hand-skip A->AA , A->ZA, A->ZZ, this needs to be better addressed by self.Model.CanMakeVertex 
                        ##

                        if DecayChannel.Particles[0] == self.Model.ParticleContent["A"] and DecayChannel.Particles[1] == self.Model.ParticleContent["A"] and DecayChannel.Particles[2] == self.Model.ParticleContent["A"]:
                            continue
                        if self.Model.ParticleContent["Z"] in DecayChannel.Particles and self.Model.ParticleContent["A"] in DecayChannel.Particles:
                            continue
                        
                        if not DecayChannel.IsPossible():
                            continue
                        
                        if Verbose:
                            print "      Cosidering Radiation+Emitter: "+REParticle.nam

                        ##
                        ## Now that we have built and veto the possible channels we 
                        ## proceed to build the dependency tree for this radiative
                        ##

                        def MyF(p):
                            return -p.pid

                        UBList = []
                        UBMap = {}
                        Count = 0
                        for Index in range(len(RadiativeChannel.Particles)):
                            BornParticle = RadiativeChannel.Particles[Index]
                            UBList.append(BornParticle)
                            if Index != Ri:
                                # if Index != Ei:
                                UBMap[Index] = Count
                                Count +=1
                        UBList[Ei] = REParticle
                        del UBList[Ri]

                        Initial = UBList[:self.RadiProc.lni]
                        Final = UBList[self.RadiProc.lni:]
                        Initial.sort(key=MyF)
                        Final.sort(key=MyF)
                        
                        ##
                        ##  There is a lot going on here and it would be better to get 
                        ##  the c++ OLP to do all this sorting when called...
                        ##  I think for now this is ok, however we should rework this later
                        ##

                        def BuildMap(List,offset=0):
                            MappingList = []
                            Counter = offset
                            for Particle in List:
                                MappingList.append([Particle,Counter])
                                Counter+=1
                            def SortFun(p):
                                return -p[0].pid
                            MappingList.sort(key=SortFun)
                            Mapping = {}
                            Count = offset
                            for Entry in MappingList:
                                Mapping[Count]=Entry[1]
                                Count += 1
                            return Mapping

                        Mapping = BuildMap(UBList[:self.RadiProc.lni])
                        Mapping.update(BuildMap(UBList[self.RadiProc.lni:],self.RadiProc.lni))
                        UnderlyingBorn = Channel(Initial,Final,self.Model)

                        ## Unsorted
                        # UnderlyingBorn = Channel(UBList[:self.RadiProc.lni],UBList[self.RadiProc.lni:],self.Model)

                        if not UnderlyingBorn.IsPossibleAtTreeLevel() or not DecayChannel.IsPossibleAtTreeLevel():
                            continue

                        ## 
                        ## Now we require that the Underlying born configuration
                        ## exists in the expansion of the Born tags
                        ## 

                        TAG=""
                        for Particle in Final:
                            TAG+=Particle.nam

                        if TAG not in ALLOWEDTAGS:
                            continue

                        if self.Model.Fundamentals['g'] in DecayChannel.Particles:
                            Dipoles.append(QCDDipole(Type,ParentIndices,UBMap,Mapping,RadiativeChannel,DecayChannel,UnderlyingBorn,self.Model))

                        if self.Model.Fundamentals['A'] in DecayChannel.Particles:
                            Dipoles.append(EWKDipole(Type,ParentIndices,UBMap,Mapping,RadiativeChannel,DecayChannel,UnderlyingBorn,self.Model))

            self.DipoleStructures[Radiative]=DipoleStructure(RadiativeChannel,Dipoles)       

    ##
    ##  OLP Scheduler: Prepare Channels to schedule and Schedule
    ##

    def CollectChannels(self):

        ##
        ##  This collects the channels and prepares them 
        ##  to be passed down to the OLP Scheduler 
        ##

        self.Virts={}
        self.Reals={}
        self.Borns={}

        for Ch in self.BornProc.subproc:
            Ini = self.BornProc.subproc[Ch][:self.RadiProc.lni]
            Fin = self.BornProc.subproc[Ch][self.RadiProc.lni:]
            VirtChannel=Channel(Ini,Fin,self.Model)
            self.Virts[Ch]=VirtChannel

        for Ch in self.RadiProc.subproc:
            Ini = self.RadiProc.subproc[Ch][:self.RadiProc.lni]
            Fin = self.RadiProc.subproc[Ch][self.RadiProc.lni:]
            RealChannel=Channel(Ini,Fin,self.Model)
            self.Reals[Ch]=RealChannel

            for Dipole in self.DipoleStructures[Ch].Dipoles:
                Ini = Dipole.UnderlyingBorn.Particles[:self.RadiProc.lni]
                Fin = Dipole.UnderlyingBorn.Particles[self.RadiProc.lni:]
                def Key(p):
                    return -p.pid
                Ini.sort(key=Key)
                Fin.sort(key=Key)
                BornChannel=Channel(Ini,Fin,self.Model)
                HASH=""
                for Particle in Ini:
                    HASH+=Particle.nam
                HASH+="_"
                for Particle in Fin:
                    HASH+=Particle.nam
                self.Borns[HASH]=BornChannel

    def BuildOLP(self,GenerateCode):
        
        ## Here we will collect the Channels 
        ## and pass them to the OLP Schedulers 
        ## and Interface constructors to weave 
        ## it all up

        self.CollectChannels()

        ## 
        ##  Generate code using NLOX : Since NLOX is an analytical code 
        ##  we fisrt generate the source code and then build the OLP interface
        ##

        self.nlox_process = self.TplDir+'/nlox_process_tpl.h'
        self.nlox_olp = self.TplDir+'/NLOX_OLP_tpl.h'
        self.NLOX_OLP = NLOX_OLP(self.nlox,self.Virts,self.Reals,self.Borns,self.Model,self.seed_tpl,self.nlox_process,self.nlox_olp,self.SrcDir)
        self.NLOX_OLP.WriteSeeds()
        if GenerateCode:
            self.NLOX_OLP.GenerateCode()
            self.NLOX_OLP.GenerateInterface()
        self.NLOX_OLP.WriteOLPClass()

        ## 
        ##  Generate RECOLA inerface : Since RECOLA is numerical 
        ##  we simply need to write the OLP interface and the actual 
        ##  generation of code will be done at runtime
        ##

        self.recola_olp = self.TplDir+'/RECOLA_OLP_tpl.h'
        self.RECOLA_OLP = RECOLA_OLP(self.recola,self.Virts,self.Reals,self.Borns,self.Model,self.recola_olp,self.SrcDir)
        self.RECOLA_OLP.WriteOLPClass()

    ##
    ##  Build the Interface 
    ##
    
    def CreateDirectories(self):
        MakeDir(self.SrcDir)
        MakeDir(self.SrcDir+'/NLOX_Process')
        MakeDir(self.SrcDir+'/Code')

    def BuildInterface(self):

        ## 
        ## OLD Deprecated, this lightens the load on generated processes.
        ## They now all point to the same directory inside src.
        ## This could possibly impede the swift test at a PSP like we were 
        ## doing by simply editting Virtual_Structure's BGenerate to be 
        ## a cosntant. This also makes it more difficult to have process dependent
        ## PSP remappings. 
        ##
        ## The solution to this would be to have an implementation file per process 
        ## to define an in-house PSP Generator, then each integrand can use the 
        ## implementation it pleases, this would be controlled by Vrtual.h and Real.h 
        ##
        ## The makefiles have also been edited for this reason
        ##

        ## CODEFILES = ['Dipole_Structure.h','Dipole_Definitions.h', \
        ##              'PSP_Generator.h','Utilities.h','Virtual_Structure.h',\
        ##              'Montecarlo_Integrator.h','GSL_Integrator.h','CUBA_Integrator.h', \
        ##              'PDF_Set.h','LHA_PDF.h','Kinematics.h','Analysis.h', \
        ##              'Dummy_PDF.h' ,\
        ##              'Constants.h','XSection.h','XSection_Integrator.h',\
        ##              'Integrand.h','OLP.h','Four_Vector.h','Input.h']

        INPUTFILES = ['Run_Settings.input']

        TOPLAYERFILES  = ['Main.cpp','Analysis.cpp','Model.cpp']

        ## for file in CODEFILES:
        ##     CopyFile(self.IntDir+'/'+file,self.SrcDir+'/Code/'+file)
        
        for file in TOPLAYERFILES:
            CopyFile(self.IntDir+'/'+file,self.SrcDir+'/'+file)

        for file in INPUTFILES:
            CopyFile(self.Inputs+'/'+file,self.SrcDir+'/'+file)

        self.LoadPaths()
        self.paths['NLOX PATH']=self.config['NLOX PATH'][0]
        self.paths['RECOLA PATH']=self.config['RECOLA PATH'][0]
        self.paths['MADISQE PATH']=os.getcwd()
        
        Makefile = Template(self.TplDir+'/makefile',self.SrcDir+'/makefile',self.paths)
        Makefile.Write()
        
    ##
    ##  Build the Integrands Classes: Virtual and Real
    ##

    def BuildRealIntegrands_New(self):
        for DipoleStructure in self.DipoleStructures:
            ChlDir = self.SrcDir+'/Real/'+str(self.DipoleStructures[DipoleStructure].Radiative)
            MakeDir(ChlDir)
            self.DipoleStructures[DipoleStructure].Show(self.TplDir,ChlDir)


    def BuildRealIntegrands(self):

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

            ClassName = Radiative+'_Real'

            INTDICT['Include Integrands'] += '#include "../Real/'+Radiative+'/'+ClassName+'.h"\n'
            INTDICT['Integrand Catalogue'] += TAB9+'Channels['+str(Count)+'] = new '+ClassName+'(*Provider,*model);\n'
            INTDICT['Integrand Catalogue'] += TAB9+'ChannelMap.insert({"'+Radiative+'",'+str(Count)+'});\n'

            Count += 1
            ChlDir = self.SrcDir+'/Real_Old/'+Radiative
            MakeDir(ChlDir)

            DICT={'SubProcHeader'    : '' ,\
                  'SubProcConst'     : '' ,\
                  'SubProcName'      : '' ,\
                  'Next'             : str(self.RadiProc.nex) ,\
                  'SubProcMat'       : '' ,\
                  'SubProcSub'       : '' ,\
                  'SubProcPlu'       : '' ,\
                  'SubProcEnd'       : '' ,\
                  'Model'            : self.Model.Name}
            
            DICT['SubProcHeader'] = Radiative.upper()
            DICT['SubProcName'] = ClassName

            ## 
            ##  Integrands constructor 
            ##

            count = 0 
            for particle in self.RadiProc.subproc[Radiative]:
                DICT['SubProcConst'] += TAB3+'Masses['+str(count)+'] = model->'+particle.nam+'.Mass;\n'
                DICT['SubProcConst'] += TAB3+'PID['+str(count)+'] = model->'+particle.nam+'.PID;\n'
                count += 1
    
            DICT['SubProcConst'] += '\n'   
            DICT['SubProcConst'] += TAB3+'nBorn = '+str(len(self.Linked[Radiative]))+';\n'   
            DICT['SubProcConst'] += TAB3+'BornMasses = new double* [nBorn];\n'
            DICT['SubProcConst'] += TAB3+'BornPID = new int* [nBorn];\n' 
            DICT['SubProcConst'] += TAB3+'for(int i=0;i<nBorn;i++) BornPID[i]= new int[NextR-1];\n\n'  
            DICT['SubProcConst'] += TAB3+'for(int i=0;i<nBorn;i++) BornMasses[i]= new double[NextR-1];\n\n' 
             
            Bcount = 0
            for Born in self.Linked[Radiative]:
                count = 0
                DICT['SubProcConst'] += TAB3+'BornMap.insert({"'+Born+'",'+str(Bcount)+'});\n'
                for particle in self.Linked[Radiative][Born]:
                    DICT['SubProcConst'] += TAB3+'BornMasses['+str(Bcount)+']['+str(count)+']= model->'+particle.nam+'.Mass;\n'
                    DICT['SubProcConst'] += TAB3+'BornPID['+str(Bcount)+']['+str(count)+']= model->'+particle.nam+'.PID;\n'
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

            CPHEAD = ''

            for CP in RadiCPS:
                
                if len(CPHEAD):
                    CPHEAD = TAB3+'else if (cp =="'+CP+'"){\n'
                else:
                    CPHEAD = TAB3+'if (cp =="'+CP+'"){\n'
                
                CACHEDTAG = ''

                if len(EWKDipoles):
                    CPHEAD += TAB6+'double born[3];\n'
                if len(QCDDipoles):
                    CPHEAD += TAB6+'double borncc['+str(1+(self.RadiProc.nex-1)*(self.RadiProc.nex-2)/2)+'];\n'
                    CPHEAD += TAB6+'ColorAndSpinMatrix BornCCSC('+str(1+(self.RadiProc.nex-1)*(self.RadiProc.nex-2)/2)+');\n'
                    
                    
                DICT['SubProcSub'] += CPHEAD
                DICT['SubProcSub'] += TAB6+'i = Proc->AmpMap.at("'+Radiative+'");\n'
                DICT['SubProcSub'] += TAB6+'for(int j=0;j<NextR;j++){pp[5*j]=p.at(j).p0;pp[5*j+1]=p.at(j).p1;pp[5*j+2]=p.at(j).p2;pp[5*j+3]=p.at(j).p3;pp[5*j+4]=0.0;}\n'
                DICT['SubProcSub'] += TAB6+'Proc->evaluate_alpha(i,"tree_tree","'+CP+'",pp,'+str(Next)+',mu,radiative,&acc);\n'
                DICT['SubProcSub'] += TAB6+'*rval = radiative[2];\n\n'

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

                        QI=0
                        QJ=0

                        SigmaIJ = ( -1 if Emitter['SUBTYP'] != Spectator['SUBTYP'] else 1 )

                        try:                            
                            QI = self.RadiProc.subproc[Radiative][I].sym['Charge']
                            QJ = self.RadiProc.subproc[Radiative][J].sym['Charge']
                        except:
                            continue
                        
                        if CACHEDTAG != Emitter['SBORNTAG']:
                            AMPMAPLINE   = TAB6+'i = Proc->AmpMap.at("'+Emitter['SBORNTAG']+'");\n'
                            BORNMAPLINE  = TAB6+'j = BornMap.at("'+Emitter['SBORNTAG']+'");\n'
                            BORNMAPLINE += TAB6+'BGenerate("'+Emitter['SBORNTAG']+'",sqrts,rand,&J);\n'
                            BORNMAPLINE += TAB6+'for(int j=0;j<NextR-1;j++){pp[5*j+0]=BornMomenta[j].p0;pp[5*j+1]=BornMomenta[j].p1;pp[5*j+2]=BornMomenta[j].p2;pp[5*j+3]=BornMomenta[j].p3;pp[5*j+4]=0.0;}\n'
                            CACHEDTAG = Emitter['SBORNTAG']
                        else:
                            AMPMAPLINE = ''
                            BORNMAPLINE = ''

                        ##
                        ## There is a factor of 3 between each charge in StandardModel.py vs RealWorld!
                        ## Hence the (1/9) for squared charges. 
                        ## StandardModel uses Charges*3 to enforce charge conservation using integers.
                        ##

                        AMPMAPLINE  += TAB6+'DipFac = '+str(SigmaIJ*QI*QJ)+'.0/9*EWKFac;\n'
                        
                        ## 
                        ##  Subtracted Function 
                        ## 

                        PREFIX = Emitter['SUBTYP']+Spectator['SUBTYP']
                        DIPFUN = DIPDIC[PREFIX]
                        TYP = Emitter['PARTYPS']

                        DICT['SubProcSub'] += AMPMAPLINE
                        DICT['SubProcSub'] += TAB6+'Build_'+PREFIX+'_Momenta(p,&p_tilde,'+str(I)+','+str(J)+','+str(K)+');\n'
                        DICT['SubProcSub'] += TAB6+'for(int j=0;j<NextR-1;j++){pp_tilde[5*j+0]=p_tilde.at(j).p0;pp_tilde[5*j+1]=p_tilde.at(j).p1;pp_tilde[5*j+2]=p_tilde.at(j).p2;pp_tilde[5*j+3]=p_tilde.at(j).p3;pp_tilde[5*j+4]=0.0;}\n'
                        DICT['SubProcSub'] += TAB6+'Proc->evaluate_alpha(i,"tree_tree","'+CPB+'",pp_tilde,'+str(Next-1)+',mu,born,&acc);\n'
                        DICT['SubProcSub'] += TAB6+'*rval += DipFac*g_'+DIPFUN+'_'+TYP+'(p['+str(I)+'],p['+str(J)+'],p['+str(K)+'],Masses['+str(I)+'],Masses['+str(J)+'])*born[2];\n\n'
                        
                        ##
                        ##  Plus Distribution
                        ##

                        DICT['SubProcPlu'] += AMPMAPLINE
                        DICT['SubProcPlu'] += BORNMAPLINE
                        DICT['SubProcPlu'] += TAB6+'Invariant = BornMasses[j]['+str(I)+']*BornMasses[j]['+str(I)+']'\
                                                              '+ BornMasses[j]['+str(J)+']*BornMasses[j]['+str(J)+']'\
                                                              +('+' if (PREFIX=='II'or PREFIX=='FF') else '-')+\
                                                              '2*(BornMomenta['+str(I)+']*BornMomenta['+str(J)+']);\n'
                        DICT['SubProcPlu'] += TAB6+'Proc->evaluate_alpha(i,"tree_tree","'+CPB+'",pp,'+str(Next-1)+',mu,born,&acc);\n'
                        DICT['SubProcPlu'] += TAB6+'CurlyG_'+DIPFUN+'_'+TYP+'(Invariant,BornMasses[j]['+str(I)+'],BornMasses[j]['+str(J)+'],mu,aux);\n'
                        DICT['SubProcPlu'] += TAB6+'for(int k=0;k<=2;k++) rval[k] += DipFac*aux[k]*born[2];\n\n'
                        DICT['SubProcPlu'] += TAB6+';\n'
                        
                        ##Channel
                        ##  Endpoint Function
                        ##

                        DICT['SubProcEnd'] += AMPMAPLINE
                        DICT['SubProcEnd'] += BORNMAPLINE
                        DICT['SubProcEnd'] += TAB6+'Invariant = BornMasses[j]['+str(I)+']*BornMasses[j]['+str(I)+']'\
                                                              '+ BornMasses[j]['+str(J)+']*BornMasses[j]['+str(J)+']'\
                                                              +('+' if (PREFIX=='II'or PREFIX=='FF') else '-')+\
                                                              '2*(BornMomenta['+str(I)+']*BornMomenta['+str(J)+']);\n'
                        DICT['SubProcEnd'] += TAB6+'Proc->evaluate_alpha(i,"tree_tree","'+CPB+'",pp,'+str(Next-1)+',mu,born,&acc);\n'
                        DICT['SubProcEnd'] += TAB6+'G_'+DIPFUN+'_'+TYP+'(Invariant,BornMasses[j]['+str(I)+'],BornMasses[j]['+str(J)+'],mu,aux);\n'
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

                        ## Offset to shift the color-correlator counter to the correct value

                        if K < FIRST:
                            FIRST -= 1
                        if K < SECON:
                            SECON -= 1

                        if CACHEDTAG != Emitter['SBORNTAG']:
                            AMPMAPLINE = TAB6+'i = Proc->AmpMap.at("'+Emitter['SBORNTAG']+'");\n'
                            BORNMAPLINE  = TAB6+'j = BornMap.at("'+Emitter['SBORNTAG']+'");\n'
                            BORNMAPLINE += TAB6+'BGenerate("'+Emitter['SBORNTAG']+'",sqrts,rand,&J);\n'
                            BORNMAPLINE += TAB6+'for(int j=0;j<NextR-1;j++){pp[5*j+0]=BornMomenta[j].p0;pp[5*j+1]=BornMomenta[j].p1;pp[5*j+2]=BornMomenta[j].p2;pp[5*j+3]=BornMomenta[j].p3;pp[5*j+4]=0.0;}\n'
                            CACHEDTAG = Emitter['SBORNTAG']
                        else:
                            AMPMAPLINE = ''
                            BORNMAPLINE = ''

                        ## This unflattens a symmetric (non-repeating) matrix, which borncc is

                        WHICHCC = (self.RadiProc.nex-2)*FIRST + SECON - FIRST - FIRST*(FIRST-1)/2

                        PREFIX = Emitter['SUBTYP']+Spectator['SUBTYP']
                        DIPFUN = DIPDIC[PREFIX]
                        # TYP = ('fermion' if isinstance(Emitter['DIPOLE'][0],Fermion) else 'boson')
                        TYP     = Emitter['PARTYPS'] 
                        BORNTYP = ('cc'     if TYP == 'ffb' else 'cc_and_sc')
                        BORNPTR = ('borncc' if TYP == 'ffb' else str(I)+',BornCCSC.ccsc')
                        TR      = (''       if TYP == 'ffb' else 'Trace')
                        INI     = (''       if TYP == 'ffb' else 'FMatrixT<double>')
                        BORN    = ('borncc' if TYP == 'ffb' else 'BornCCSC.ccsc')

                        ## 
                        ##  Subtracted Function 
                        ## 

                        DICT['SubProcSub'] += AMPMAPLINE
                        DICT['SubProcSub'] += TAB6+'Build_'+PREFIX+'_Momenta(p,&p_tilde,'+str(I)+','+str(J)+','+str(K)+');\n'
                        DICT['SubProcSub'] += TAB6+'for(int j=0;j<NextR-1;j++){pp_tilde[5*j+0]=p_tilde.at(j).p0;pp_tilde[5*j+1]=p_tilde.at(j).p1;pp_tilde[5*j+2]=p_tilde.at(j).p2;pp_tilde[5*j+3]=p_tilde.at(j).p3;pp_tilde[5*j+4]=0.0;}\n'
                        DICT['SubProcSub'] += TAB6+'Proc->evaluate_alpha_'+BORNTYP+'(i,"tree_tree","'+CPB+'",pp_tilde,'+str(Next-1)+','+BORNPTR+');\n'
                        DICT['SubProcSub'] += TAB6+'*rval += QCDFac*'+TR+'(g_'+DIPFUN+'_'+TYP+'(p['+str(I)+'],p['+str(J)+'],p['+str(K)+'],Masses['+str(I)+'],Masses['+str(J)+'])*'+INI+'('+BORN+'['+str(WHICHCC)+']));\n\n'
                        
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
                        DICT['SubProcEnd'] += TAB6+'Proc->evaluate_alpha_cc(i,"tree_tree","'+CPB+'",pp,'+str(Next-1)+',borncc);\n'
                        DICT['SubProcEnd'] += TAB6+'G_'+DIPFUN+'_'+TYP+'(Invariant,BornMasses[j]['+str(I)+'],BornMasses[j]['+str(J)+'],mu,aux);\n'
                        DICT['SubProcEnd'] += TAB6+'for(int k=0;k<=2;k++) rval[k] += QCDFac*aux[k]*borncc['+str(WHICHCC)+'];\n\n'
                        

                DICT['SubProcSub'] += TAB3+'}\n'
                DICT['SubProcPlu'] += TAB3+'}\n'
                DICT['SubProcEnd'] += TAB3+'}\n'
            

            SubFunctionH = Template(self.TplDir+'/Dipole_Functions_tpl.h',ChlDir+'/'+Radiative+'_Real.h',DICT)
            SubFunctionC = Template(self.TplDir+'/Dipole_Functions_tpl.cpp',ChlDir+'/'+Radiative+'_Real.cpp',DICT)

            SubFunctionH.Write()
            SubFunctionC.Write()

        RealIntegrand = Template(self.TplDir+'/Real_tpl.h',self.SrcDir+'/Code/Real.h',INTDICT)
        RealIntegrand.Write()

    def BuildVirtualIntegrands(self):

        TAB3 = '    '
        TAB6 = TAB3+TAB3
        TAB9 = TAB6+TAB3

        VIRTDICT = {'Include Virtuals' : '' ,\
                    'nChannels'         : str(len(self.Borns)) ,\
                    'Virtual Catalogue': ''}

        Count = 0
        for Born in self.Borns:

            ClassName = Born+'_Virtual'

            VIRTDICT['Include Virtuals'] += '#include "../Virtual/'+Born+'/'+ClassName+'.h"\n'
            VIRTDICT['Virtual Catalogue'] += TAB9+'Channels['+str(Count)+'] = new '+ClassName+'(*Provider,*model);\n'
            VIRTDICT['Virtual Catalogue'] += TAB9+'ChannelMap.insert({"'+Born+'",'+str(Count)+'});\n'

            Count += 1
            ChlDir = self.SrcDir+'/Virtual/'+Born
            MakeDir(ChlDir)

            DICT={'SubProcHeader'    : '' ,\
                  'SubProcConst'     : '' ,\
                  'SubProcName'      : '' ,\
                  'Next'             : str(self.RadiProc.nex-1) ,\
                  'SubProcBorn'      : '',\
                  'SubProcVirt'      : ''}

            DICT['SubProcHeader'] = Born.upper()
            DICT['SubProcName'] = ClassName
 
            count = 0
            for particle in self.Borns[Born].Particles:
                DICT['SubProcConst'] += TAB3+'Particles['+str(count)+']= &(model->'+particle.nam+');\n'
                count += 1
            
            BornCPS = self.Model.GetCPS(self.Borns[Born].Particles)
            VirtCPS = {}
            for CP in BornCPS:
                asp = BornCPS[CP]['as']
                aep = BornCPS[CP]['ae']
                asc = 'as'+str(asp+1)+'ae'+str(aep)
                aec = 'as'+str(asp)+'ae'+str(aep+1)
                VirtCPS[asc] = {'as':asp+1,'ae':aep}
                VirtCPS[aec] = {'as':asp,'ae':aep+1}

            CPHEAD = ''
            for CP in BornCPS:    
                if len(CPHEAD):
                    CPHEAD = TAB3+'else if (cp =="'+CP+'"){\n'
                else:
                    CPHEAD = TAB3+'if (cp =="'+CP+'"){\n'

                DICT['SubProcBorn'] += CPHEAD
                DICT['SubProcBorn'] += TAB6+'Args.SubProc = "'+Born+'";\n'
                DICT['SubProcBorn'] += TAB6+'Args.CouplingPower[0] = '+str(BornCPS[CP]['as'])+';\n'
                DICT['SubProcBorn'] += TAB6+'Args.CouplingPower[1] = '+str(BornCPS[CP]['ae'])+';\n'
                DICT['SubProcBorn'] += TAB6+'Provider->Evaluate(&Args);\n'
                DICT['SubProcBorn'] += TAB6+'}\n'
            
            CPHEAD = ''  
            for CP in VirtCPS:
                if len(CPHEAD):
                    CPHEAD = TAB3+'else if (cp =="'+CP+'"){\n'
                else:
                    CPHEAD = TAB3+'if (cp =="'+CP+'"){\n'

                DICT['SubProcVirt'] += CPHEAD
                DICT['SubProcVirt'] += TAB6+'Args.SubProc = "'+Born+'";\n'
                DICT['SubProcVirt'] += TAB6+'Args.CouplingPower[0] = '+str(VirtCPS[CP]['as'])+';\n'
                DICT['SubProcVirt'] += TAB6+'Args.CouplingPower[1] = '+str(VirtCPS[CP]['ae'])+';\n'
                DICT['SubProcVirt'] += TAB6+'Provider->Evaluate(&Args);\n'
                DICT['SubProcVirt'] += TAB6+'}\n'
    
            SubFunctionH = Template(self.TplDir+'/Virtual_Functions_tpl.h',ChlDir+'/'+Born+'_Virtual.h',DICT)
            SubFunctionC = Template(self.TplDir+'/Virtual_Functions_tpl.cpp',ChlDir+'/'+Born+'_Virtual.cpp',DICT)
            
            SubFunctionH.Write()
            SubFunctionC.Write()

        VirtualIntegrand = Template(self.TplDir+'/Virtual_tpl.h',self.SrcDir+'/Code/Virtual.h',VIRTDICT)
        VirtualIntegrand.Write()

def main(CONFIGFILE):

    print'''
####################################################################
#                                                                  #        
#      __  __               _   _                                  #
#     |  \/  |   __ _    __| | (_)  ___    __ _    ___             #
#     | |\/| |  / _` |  / _` | | | / __|  / _` |  / _ \            #
#     | |  | | | (_| | | (_| | | | \__ \ | (_| | |  __/            #
#     |_|  |_|  \__,_|  \__,_| |_| |___/  \__, |  \___|            #
#                                            |_|                   #
#                                                                  #        
#                           Version 1.0.0                          #                   
#                          by D. Figueroa                          #                   
#                                                                  #        
#                                                                  #        
####################################################################
'''
    MADISQE = Madisqe(CONFIGFILE)

if __name__ == "__main__":
    if len(sys.argv)<2:
        print "\33[31mError\33[0m: No input file specified"
    elif len(sys.argv)==2:
        main(sys.argv[1])
    else:
        print "\33[31mError\33[0m: Too many arguments"  
        
    