from Physics import *
from Template import *

class OLP:
    
    template = ''

    def Write():
        ChannelDict = {'BuildChannels':''}
        for Ch in self.ChannelMap:
            ChannelDict += 'ChannelMap.insert({"'+Ch+'",'+self.ChannelMap[Ch]+'})\n'

class NLOX_OLP(OLP):
    def __init__(self,NLOX,Virts,Reals,Borns,Model,Seed,Process,Path):
        
        self.nlox=NLOX
        self.Virts=Virts
        self.Reals=Reals
        self.Borns=Borns
        self.Model=Model
        self.Path =Path
        self.Seed=Seed
        self.Process=Process
        self.ConfigureChannels()

    def ChannelToDict(self,Channel,Type):
        BASEDICT = {'Initial State' : '' ,\
                    'Final State'   : '' ,\
                    'Tree CP'       : 'bornAlphaQCD = ' ,\
                    'Loop CP'       : '## virtAlphaQCD = ' ,\
                    'enableTerms'   : '4' ,\
                    'ColorCorr'     : 'false' ,\
                    'SpinCorr'      : 'false' ,\
                    'Path'          : '' }

        count = 1
        for Particle in Channel.Initial:
            BASEDICT['Initial State'] += Particle.nam
            if count != len(Channel.Initial):
                BASEDICT['Initial State'] += ','
            count += 1
        count = 1
        for Particle in Channel.Final:
            BASEDICT['Final State'] += Particle.nam
            if count != len(Channel.Final):
                BASEDICT['Final State'] += ','
            count += 1

        BASEDICT['Path'] = self.Path
        BASEDICT['Tree CP'] += CSS(Channel.MaxQCD,Channel.MinQCD)
        
        if(Type=='Virt'):
            BASEDICT['ColorCorr'] = 'true'
            BASEDICT['SpinCorr'] = 'true'
            BASEDICT['enableTerms'] = '7'
            BASEDICT['Loop CP']  = 'virtAlphaQCD = '
            BASEDICT['Loop CP'] += CSS(Channel.MaxQCD+1,Channel.MinQCD)

        if(Type=='Born'):
            BASEDICT['ColorCorr'] = 'true'
            BASEDICT['SpinCorr'] = 'true'

        return BASEDICT

    def ConfigureChannels(self,Fixed=False):
        
        ## We have three kind of Channels to Schedule to the OLP
        ## Virts -> 2 to n @NLO Channels
        ## Reals -> 2 to n+1 @LO Channles
        ## Borns -> 2 to n @LO Channels -> Both Spin and Color Correlated 

        self.Channels={}

        for Ch in self.Borns:
            BornDict = self.ChannelToDict(self.Borns[Ch],'Born')
            self.Channels[Ch]={'Channel':self.Borns[Ch],'Dictionary':BornDict,'Type':'Born'}
        
        for Ch in self.Reals:
            RealDict = self.ChannelToDict(self.Reals[Ch],'Real')
            self.Channels[Ch]={'Channel':self.Reals[Ch],'Dictionary':RealDict,'Type':'Real'}

        for Ch in self.Virts:
            VirtDict = self.ChannelToDict(self.Virts[Ch],'Virt')
            self.Channels[Ch]={'Channel':self.Virts[Ch],'Dictionary':VirtDict,'Type':'Virt'}

        self.ReducedChannels={}
        self.ChannelMap={}

        for Ch in self.Channels:
            self.ReducedChannels[Ch]=self.Channels[Ch]
            self.ChannelMap[Ch]=self.Channels[Ch]

        ## We need to add some padding to avoid sending
        ## unecessary resquests to the OLP. The padding 
        ## will work as follows:
        ##        - The c++ OLP object has all the Channels listed
        ##        - It also has a map from the Channels to the subset requested to the OLP
        
    def WriteSeeds(self):
        for Ch in self.ReducedChannels:
            Dictionary = self.ReducedChannels[Ch]['Dictionary']
            FileName=self.Path+'/'+Ch+'.in'
            ChannelSeed = Template(self.Seed,FileName,Dictionary)
            ChannelSeed.Write()

    def Generate(self):
        Here = os.getcwd()
        os.chdir(self.Path)
        
        count = 1
        tot = str(len(self.ReducedChannels))
        for Ch in self.ReducedChannels:
            print 'Generating',self.ReducedChannels[Ch]['Type'],'channel:',Ch,'('+str(count)+'/'+tot+')'
            cmd = self.nlox+' '+Ch +'.in > '+Ch+'.log'
            os.system(cmd)
            cmd = 'mv '+Ch+'.log '+self.Path+'/'+Ch
            os.system(cmd)
            cmd = 'mv '+Ch+'.in '+self.Path+'/'+Ch
            os.system(cmd)
            count += 1

        os.chdir(Here)

    def GenerateInterface(self):
        INTERFACEDICT = {'Include Channels'   : '' ,\
                         'NSubProcesses'      : str(len(self.ReducedChannels)) ,\
                         'Construct Channels' : '\n' ,\
                         'Amplitude Map'      : '\n'}
        TAB3 = '    '
        TAB6 = TAB3+TAB3
        
        count = 0
        for Ch in self.ReducedChannels:
            INTERFACEDICT['Include Channels']   += '#include "../'+Ch+'/code/Sub_'+Ch+'.h"\n'
            INTERFACEDICT['Construct Channels'] += TAB6+'subproc['+str(count)+'] = '+'new Sub_'+Ch+'(pc);\n'
            INTERFACEDICT['Amplitude Map']      += TAB6+'AmpMap.insert({"'+Ch+'",'+str(count)+'});\n'
            count += 1

        FILENAME = self.Path+'/code/nlox_process.h'
        cmd = 'rm '+FILENAME
        os.system(cmd)
        NLOXPROCESS = Template(self.Process,FILENAME,INTERFACEDICT)
        NLOXPROCESS.Write()

class Recola_OLP(OLP):
    pass