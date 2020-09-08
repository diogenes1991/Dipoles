
from Particle_Engine import *
import gc
        
class StandardModel(Model):
    def __init__(self):

        self.name = 'Standard Model'

        ##
        ## Bosons
        ## 

        g  = Particle("g", {},0,'Boson')
        A  = Particle("A", {},22,'Boson')
        Z  = Particle("Z", {},23,'Boson')
        Wp = Particle("Wp",{'Charge': 3},24,'Boson')
        Wm = Particle("Wm",{'Charge':-3},-24,'Boson')
        h  = Particle("h", {},25,'Boson')
        
        ##
        ## Fermions: Quarks
        ## 

        d = Particle("d",{'udness':1,'Charge':-1},1,'Fermion')
        u = Particle("u",{'udness':1,'Charge': 2},2,'Fermion')
        s = Particle("s",{'scness':1,'Charge':-1},3,'Fermion')
        c = Particle("c",{'scness':1,'Charge': 2},4,'Fermion')
        b = Particle("b",{'btness':1,'Charge':-1},5,'Fermion')
        t = Particle("t",{'btness':1,'Charge': 2},6,'Fermion')

        dbar = Particle("dbar",{'udness':-1,'Charge': 1},-1,'Fermion')
        ubar = Particle("ubar",{'udness':-1,'Charge':-2},-2,'Fermion')
        sbar = Particle("sbar",{'scness':-1,'Charge': 1},-3,'Fermion')
        cbar = Particle("cbar",{'scness':-1,'Charge':-2},-4,'Fermion')
        bbar = Particle("bbar",{'btness':-1,'Charge': 1},-5,'Fermion')
        tbar = Particle("tbar",{'btness':-1,'Charge':-2},-6,'Fermion')

        self.quarks = set([u,d,c,s,b,t,ubar,dbar,cbar,sbar,bbar,tbar])

        ##
        ## Fermoions: Leptons
        ##
        
        em   = Particle("em",  {'eness':1,'Charge':-3},11,'Fermion')
        ne   = Particle("ne",  {'eness':1,'Charge': 0},12,'Fermion')
        mum  = Particle("mum", {'mness':1,'Charge':-3},13,'Fermion')
        nm   = Particle("nm",  {'mness':1,'Charge': 0},14,'Fermion')
        taum = Particle("taum",{'tness':1,'Charge':-3},15,'Fermion')
        nt   = Particle("nt",  {'tness':1,'Charge': 0},16,'Fermion')

        ep    = Particle("ep",   {'eness':-1,'Charge':3},-11,'Fermion')
        nebar = Particle("nebar",{'eness':-1,'Charge':0},-12,'Fermion')
        mup   = Particle("mup",  {'mness':-1,'Charge':3},-13,'Fermion')
        nmbar = Particle("nmbar",{'mness':-1,'Charge':0},-14,'Fermion')
        taup  = Particle("taup", {'tness':-1,'Charge':3},-15,'Fermion')
        ntbar = Particle("ntbar",{'tness':-1,'Charge':0},-16,'Fermion')

        self.leptons = set([ep,ep,mup,mum,taup,taum,ne,nebar,nm,nmbar,nt,ntbar])

        self.QCDPars = set([g])
        self.QCDPars = self.QCDPars.union(self.quarks)  
        
        self.EWKPars = set([A,Z,Wp,Wm,h])
        self.EWKPars = self.EWKPars.union(self.quarks)
        self.EWKPars = self.EWKPars.union(self.leptons)      

        ########################################################################
        ## Declare composite particles:

        ## Composite particles are declared with a name and a list of fundmental 
        ## particles

        p = CompositeParticle("p",[d,dbar,u,ubar,s,sbar,c,cbar,b,bbar,g,A])

        j = CompositeParticle("j",[d,dbar,u,ubar,s,sbar,c,cbar,b,bbar,g,A])
        
        lp = CompositeParticle("lp",[ep,mup,taup])
        
        lm = CompositeParticle("lm",[em,mum,taum])
        
        Zj = CompositeParticle("Zj",[d,dbar,u,ubar,s,sbar,c,cbar,b,bbar,ep,em,mup,mum,taup,taum,ne,nm,nt,nebar,nmbar,ntbar])
        # Zj.DaughterOf(Z)
        
        Wj = CompositeParticle("Wj",[d,dbar,u,ubar,s,sbar,c,cbar,b,bbar,ep,em,mup,mum,taup,taum,ne,nm,nt,nebar,nmbar,ntbar])
        # Wj.DaughterOf(Wp)
        
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
