# Construction of the EWK Subtracted Function
# This uses auto coupling-power detection

EWKc = 0
QCDc = 0
for particle in self.RadiProc.subproc[subproc]:
    if particle in self.Model.EWKPars:
        EWKc += 1
    if particle in self.Model.QCDPars:
        QCDc += 1

MaxQCD = QCDc - 2
MinQCD = self.RadiProc.nex - EWKc   
MaxEWK = EWKc - 2
MinEWK = self.RadiProc.nex - QCDc
CPSQCD = []
for val in range(MinQCD,MaxQCD+1):
    CPSQCD.append('as'+str(val)+'ae'+str(self.RadiProc.nex-val-2))
CPSEWK = []
for val in range(MinEWK,MaxEWK+1):
    CPSEWK.append('as'+str(self.RadiProc.nex-val-2)+'ae'+str(val))
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


def WriteTests(self):
        DICT = {    'nExtRadi'         : str(self.RadiProc.nex) ,\
                    'Evaluate BornsC'  :  '' ,\
                    'Evaluate BornsF'  :  '' }
        
        ## 
        ##  Looping over Borns:
        ##
        ##  - Set masses 
        ##  - Call for generation of a PSP
        ##  - Print the PSP
        ##  - Call NLOX
        ##  - Print the result 
        ##

        for born in self.Borns:
            DICT['Evaluate BornsC'] += '/////////////////////////////////////////\n'
            DICT['Evaluate BornsC'] += '// Evaluation of block for:'+born+'\n'

            count = 0
            for particle in self.Borns[born]:
                PARTYP = ('In' if count<self.RadiProc.lni else 'FB')
                PARMAS = 'm'+(particle.nam[0] if (particle in self.Model.quarks or particle in self.Model.EWKPars) else particle.nam)
                if count==self.RadiProc.lni-1:
                    count = 0
                NUM = str(count)
                DICT['Evaluate BornsC'] += 'M'+PARTYP+'['+str(count)+'] = proc->pc.'+PARMAS+';\n'
                count += 1
            
            BornCPS = self.Model.GetCPS(self.Borns[born])
            for cp in BornCPS:
            
                DICT['Evaluate BornsC'] += 'std::cout << Now evaluating: '+born+'<< at <<'+cp+'<<std::endl;\n'
            DICT['Evaluate BornsC'] += '/////////////////////////////////////////\n'
        
        TestProcessC = seek_and_destroy(self.TplDir+'/Test_ProcessC.tpl',DICT)
        WriteFile(self.SrcDir+'/Test_Process.cc',TestProcessC)
    