import sys,os,copy

class Particle:    
    def __init__(self,Name,Sym,PID,Type):
        self.nam = Name
        self.pid = PID
        self.typ = Type
        self.sym = Sym
        
    def __str__(self):
        return self.nam

class Boson(Particle):
    def Carries(Symmetry):
        self.carrier = Symmetry

class Fermion(Particle):
    def Dummy(self):
        return 1

class CompositeParticle:
    def __init__(self,Name,Subparticles,Mother=None,Degree=2):
        self.nam = Name
        for particle in Subparticles:
            if not isinstance(particle,Particle):
                print "\33[31mError\33[0m: In composite particle",self.nam,"a non-particle constituent has been passed:",particle
                sys.exit()
        self.sub = Subparticles
        
        if Mother!=None:
            if not isinstance(Mother,Particle):
                print "\33[31mError\33[0m: In composite particle",self.nam,"a non-particle parent has been passed:",Mother
                sys.exit()
            self.parent = Mother
            self.degree = 2
            self.InheritSymmetries()
            
    def DaughterOf(self,parent,degree=2):
        if not isinstance(parent,Particle):
            print "\33[31mError\33[0m: In composite particle",self.nam,"a non-particle parent has been passed:",parent
            sys.exit()
        self.parent = parent
        self.degree = degree
        self.InheritSymmetries()
            
    def __str__(self):
        return self.nam
    
    def InheritSymmetries(self):
        self.sym_weights = {}
        for index in range(len(self.sub)):
            #
            # If one wants to make use of the abstract symmetry engine to enforce 
            # parent -> daughters relations one must create a new set abstract symmetries that all the 
            # daughters of the parent have a non-zero weight under. The values of the 
            # weights are all the non-zero values of abstract symmetries that the daughters
            # have but the parent doesn't. For example "Zj" is daughter of "Z", "Z" has zero 
            # weight under "Charge" so for all particles in "Zj" we need to add a new abstract
            # symmetry "Zness" that takes the value of the charge of the daughter, this way 
            # only subsets of daughters that have zero total charge would be considered valid 
            # subsets of the decay "Z" -> "Zj" "Zj". The generation symmetry is also present, but 
            # it is a multiple of the charge symmetry and thus redundant, of course we don't know this
            # and a priori we include both.
            #
            new_particle = copy.deepcopy(self.sub[index])
            for symmetry in self.sub[index].sym.keys():
                new_sym_nam = self.parent.nam+"_"+self.nam+"_"+symmetry
                new_particle.sym[new_sym_nam] = self.sub[index].sym[symmetry]
                try:
                    self.sym_weights[new_sym_nam] = -float(self.parent.sym[symmetry])/self.degree
                except:
                    self.sym_weights[new_sym_nam] = 0
            self.sub[index] = copy.deepcopy(new_particle)

class Process:           
    def __init__(self,Initial_State,Final_State):
        self.ini = Initial_State
        self.fin = Final_State
        self.nex = len(Initial_State)+len(Final_State)
        self.lni = len(Initial_State)
        self.lnf = len(Final_State)
        
        self.proc = []
        for i in self.ini:
            self.proc.append(i)
        for i in self.fin:
            self.proc.append(i)
           
        self.InheritSymmetries()
        self.BuildSubProcesses()
        
    def __str__(self):
        aux = ""
        for i in self.ini:
            aux += i.nam +""
        aux += "_"
        for i in self.fin:
            aux += i.nam + ""
        # aux += " ("+str(len(self.scatters))+")"
        return aux
    
    def Print(self):
        line = ""
        for subprocess in self.scatters:
            count = 0
            for particle in subprocess:
                if count == len(self.ini):
                    line += "-> "
                line += particle.nam+" "
                count += 1
            line += "\n"
        print line
     
    def write(self):
        name = str(self)
        f = open(name+".subprocesses","w+")
        
        for subprocess in self.scatters:
            w = ''
            count = 0
            l = ''
            for particle in subprocess:
                if count == len(self.ini):
                    w += '-> '
                w += particle.nam+' '
                count += 1
            w += '\n'
            f.write(w)
        print len(self.scatters),"subprocesses written"
        
    def InheritSymmetries(self): 
        self.sym_weights = {}
        for particle in self.proc:
            try: 
                particle.sym_weights
            except:
                continue
            for symmetry in particle.sym_weights.keys():
                try:
                    self.sym_weights[symmetry] += particle.sym_weights[symmetry]
                except:
                    self.sym_weights[symmetry] = particle.sym_weights[symmetry]

    def BuildScatters(self):
        EVEN_FER = True
        SORT_SYM = True
        ABST_SYM = True
        
        self.subproc = []
        
        for i in self.proc:
            if isinstance(i,CompositeParticle):
                self.subproc.append(i.sub)
            else:
                self.subproc.append([i])
       
        self.scatters = []
        l = 1
        bounds = []
        for particle in self.subproc:
            l *= len(particle)
            bounds.append(len(particle))
        
        for i in range(l):
            addres = []
            ac = i
            for j in bounds:
                addres.append(ac%j)
                ac = ac/j
                
            pars = []
            for j in range(len(bounds)):
                pars.append(self.subproc[j][addres[j]])
            
            self.scatters.append(pars)
        
        if EVEN_FER:
            aux = []
            for tent in self.scatters:
                par = 0
                for p in tent:
                    if(p.typ=='Fermion'):
                        par += 1
                par = par%2
                if(par==0):
                    aux.append(tent)
            self.scatters = aux
            
        if ABST_SYM:
            aux = []
            for tentative in self.scatters:
                
                ## Symmetries Initialization
                proc_syms = {}
                for particle in tentative:
                    for symmetry in particle.sym:
                        proc_syms[symmetry] = 0
                
                ## Computes the symmetry weight
                count = 0
                for particle in tentative:
                    for symmetry in proc_syms.keys():
                        if symmetry in particle.sym:
                            proc_syms[symmetry] += particle.sym[symmetry] if count < len(self.ini) else -particle.sym[symmetry]
                    count += 1
                
                ## Checks if all weights are zero or what the value they inherited from their parents
                respects_symmetries = True
                for symmetry in proc_syms.keys():
                    w = 0
                    try:
                        w = self.sym_weights[symmetry]
                    except:
                        w = 0
                    if proc_syms[symmetry] != w:
                        respects_symmetries = False
                
                ## If the tentative process respects all symmetries is stored 
                if respects_symmetries:
                        aux.append(tentative)
                
            self.scatters = aux       
            
        if SORT_SYM:
            def myf(e):
                return -e.pid
            
            aux = []
            for tent in self.scatters:
                sortedins = []
                sortedouts = []
                count = 0
                for par in tent:
                    if count < len(self.ini):
                        sortedins.append(par)
                    else:
                        sortedouts.append(par)
                    count += 1
                
                sortedins.sort(key=myf)        
                sortedouts.sort(key=myf)
                
                sortedins.extend(sortedouts)
                
                if len(aux):
                    shouldadd = True
                    for subproc in aux:
                        areequal = True
                        for n in range(len(subproc)):
                            if subproc[n].nam != sortedins[n].nam:
                                areequal = False
                        if areequal:
                           shouldadd = False
                    if shouldadd:
                        aux.append(sortedins)
                else:
                    aux.append(sortedins)
                    
            self.scatters = aux
            
            def SortFunction(Proc):
                n = 0
                for i in Proc:
                    if i.typ == "Boson":
                        n += len(Proc)
                    if i.pid > 0:
                        n += -1
                    else:
                        n += 1
                return n
            self.scatters.sort(key=SortFunction)

    def BuildSubProcesses(self):
        self.BuildScatters()
        self.subproc = {}
        for scatter in self.scatters:
            count = 0
            key = ""
            for particle in scatter:
                count += 1
                key += particle.nam
                if count == len(self.ini):
                    key += "_"
            self.subproc[key] = scatter

    def BuildString(self,LIST):
        key = ""
        count = 0
        for particle in LIST:
            count += 1
            key += particle.nam
            if count == len(self.ini):
                key += "_"
        return key

class Model:
    def __str__(self):
        retval = '{'
        for particle in self.ParticleContent:
            retval += str(particle) + ' '
        retval += '}'
        return retval

    def BuildModel(self):
        self.sym = {}
        # print self.ParticleContent
        # for particle in self.ParticleContent.keys():
            # if isinstance(self.ParticleContent[particle],Particle):
                # self.sym[]

        self.CrossDictionary = {}
        for particle in self.ParticleContent:
            if not isinstance(self.ParticleContent[particle],Particle):
                continue
            self.CrossDictionary[self.ParticleContent[particle].nam] = self.ParticleContent[particle].nam
            for otherparticle in self.ParticleContent:
                if not isinstance(self.ParticleContent[otherparticle],Particle):
                    continue
                if self.ParticleContent[particle].pid == -self.ParticleContent[otherparticle].pid:
                    self.CrossDictionary[self.ParticleContent[particle].nam] = self.ParticleContent[otherparticle].nam

        self.Composites = {n:p for (n,p) in self.ParticleContent.items() if isinstance(p,CompositeParticle)}
        self.Fundamentals = {n:p for (n,p) in self.ParticleContent.items() if isinstance(p,Particle)}
