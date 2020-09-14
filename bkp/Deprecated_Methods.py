def LinkProcs(self):
        self.link = {}
        for subproc in self.RadiProc.subproc:
            self.link[subproc] = {}
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

                # Otherwise we eliminate the Boson from the process and build a newlist 

                newlist = []
                for index in range(len(self.RadiProc.subproc[subproc])):
                    if index == counter:
                        continue
                    newlist.append(self.RadiProc.subproc[subproc][index])
                
                if counter >= len(self.RadiProc.ini):
                    tentborn = self.RadiProc.BuildString(newlist)
                    self.link[subproc][tentborn] = newlist
                else:
                    for index in range(len(self.RadiProc.ini),len(self.RadiProc.ini)+len(self.RadiProc.fin)-1):
                        # Copy 
                        swaped = [particle for particle in self.RadiProc.subproc[subproc]]
                        
                        # Swap 
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

def LinkProcs2(self):
        for radi in self.RadiProc.subproc:
            print radi
            for id1 in range(len(self.RadiProc.subproc[radi])):
                for id2 in range(len(self.RadiProc.subproc[radi])):
                    if id2 <= id1:
                        continue
                    # we want to pair states with total charge zero
                    sign = 1
                    if (id1 < self.RadiProc.lni and id2 < self.RadiProc.lni)\
                    or (id1 > self.RadiProc.lni and id2 > self.RadiProc.lni):
                        sign = +1
                    if (id1 < self.RadiProc.lni and id2 > self.RadiProc.lni)\
                    or (id1 > self.RadiProc.lni and id2 < self.RadiProc.lni):
                        sign = -1
                    syms = {}
                    for sym in self.RadiProc.subproc[radi][id1].sym:
                        syms[sym] = self.RadiProc.subproc[radi][id1].sym[sym]
                    for sym in self.RadiProc.subproc[radi][id2].sym:
                        if sym in syms:
                            syms[sym] += sign*self.RadiProc.subproc[radi][id2].sym[sym]
                        else:
                            syms[sym] = self.RadiProc.subproc[radi][id2].sym[sym]
                    exclude = False
                    for sym in syms:
                        if syms[sym]!=0:
                            exclude = True
                    if not exclude:
                        print self.RadiProc.subproc[radi][id1].nam,self.RadiProc.subproc[radi][id2].nam