
from Physics import *
import gc
        
class StandardModel(Model):
    def __init__(self):

        self.Name = 'StandardModel'

        ##
        ## Bosons
        ## 

        g  = Boson('g', {},0,False)
        A  = Boson("A", {},22,False)
        Z  = Boson("Z", {},23,True)
        Wp = Boson("Wp",{'Charge': 3},24,True)
        Wm = Boson("Wm",{'Charge':-3},-24,True)
        h  = Boson("h", {},25,True)

        ##
        ## Fermions: Quarks
        ## 

        d = Fermion("d",{'udness':1,'Charge':-1},1,False)
        u = Fermion("u",{'udness':1,'Charge': 2},2,False)
        s = Fermion("s",{'scness':1,'Charge':-1},3,False)
        c = Fermion("c",{'scness':1,'Charge': 2},4,False)
        b = Fermion("b",{'btness':1,'Charge':-1},5,True)
        t = Fermion("t",{'btness':1,'Charge': 2},6,True)

        dbar = Fermion("dbar",{'udness':-1,'Charge': 1},-1,False)
        ubar = Fermion("ubar",{'udness':-1,'Charge':-2},-2,False)
        sbar = Fermion("sbar",{'scness':-1,'Charge': 1},-3,False)
        cbar = Fermion("cbar",{'scness':-1,'Charge':-2},-4,False)
        bbar = Fermion("bbar",{'btness':-1,'Charge': 1},-5,True)
        tbar = Fermion("tbar",{'btness':-1,'Charge':-2},-6,True)

        self.quarks = set([u,d,c,s,b,t,ubar,dbar,cbar,sbar,bbar,tbar])

        ##
        ## Fermoions: Leptons
        ##
        
        em   = Fermion("em",  {'eness':1,'Charge':-3},11,False)
        ne   = Fermion("ne",  {'eness':1,'Charge': 0},12,False)
        mum  = Fermion("mum", {'mness':1,'Charge':-3},13,False)
        nm   = Fermion("nm",  {'mness':1,'Charge': 0},14,False)
        taum = Fermion("taum",{'tness':1,'Charge':-3},15,False)
        nt   = Fermion("nt",  {'tness':1,'Charge': 0},16,False)

        ep    = Fermion("ep",   {'eness':-1,'Charge':3},-11,False)
        nebar = Fermion("nebar",{'eness':-1,'Charge':0},-12,False)
        mup   = Fermion("mup",  {'mness':-1,'Charge':3},-13,False)
        nmbar = Fermion("nmbar",{'mness':-1,'Charge':0},-14,False)
        taup  = Fermion("taup", {'tness':-1,'Charge':3},-15,False)
        ntbar = Fermion("ntbar",{'tness':-1,'Charge':0},-16,False)

        self.leptons = set([ep,em,mup,mum,taup,taum,ne,nebar,nm,nmbar,nt,ntbar])        

        self.QCDPars = set([g])
        self.QCDPars = self.QCDPars.union(self.quarks)  
        
        self.EWKPars = set([A,Z,Wp,Wm])
        self.EWKPars = self.EWKPars.union(self.quarks)
        self.EWKPars = self.EWKPars.union(self.leptons)

        self.YUKPars = set([Wp,Wm,Z,t,tbar,b,bbar,h])

        self.Couplings = [self.QCDPars,self.EWKPars,self.YUKPars]      

        ########################################################################
        ## Declare composite particles:

        ## Composite particles are declared with a name and a list of fundmental 
        ## particles

        p = CompositeParticle("p",[d,dbar,u,ubar,s,sbar,c,cbar,b,bbar,g,A])

        j = CompositeParticle("j",[d,dbar,u,ubar,s,sbar,c,cbar,b,bbar,g,A])
        
        lp = CompositeParticle("lp",[ep,mup,taup])
        
        lm = CompositeParticle("lm",[em,mum,taum])
        
        Zj = CompositeParticle("Zj",[d,dbar,u,ubar,s,sbar,c,cbar,b,bbar,ep,em,mup,mum,taup,taum,ne,nm,nt,nebar,nmbar,ntbar])
        
        Wj = CompositeParticle("Wj",[d,dbar,u,ubar,s,sbar,c,cbar,b,bbar,ep,em,mup,mum,taup,taum,ne,nm,nt,nebar,nmbar,ntbar])
        
        jb = CompositeParticle("jb",[b,bbar])

        X = CompositeParticle("X",[g,A])

        ##
        ##
        ########################################################################

        ParticleContent = {}
        for obj in gc.get_objects():
            if isinstance(obj,Particle) or isinstance(obj,CompositeParticle):
                ParticleContent[obj.nam] = obj
        
        self.ParticleContent = ParticleContent
        self.BuildModel()

WorkingModel = StandardModel()
