from Physics import *
from Template import *

##
##  Dipoles overhead classes to build the Dipole Structures 
##

class Dipole:

    prefix = ''

    def __init__(self,DecayChannel,UnderlyingBorn):
        self.DecayChannel=DecayChannel
        self.UnderlyingBorn=UnderlyingBorn

    def __str__(self):
        Name =self.prefix+' Dipole\n'
        Name+='Decay Channel: '+str(self.DecayChannel)
        Name+='\nUnderlying Born: '+str(self.UnderlyingBorn)
        Name+='\n'
        return Name

class EWKDipole(Dipole):

    prefix = 'EWK'
    
class QCDDipole(Dipole):

    prefix = 'QCD'

class DipoleStructure:
    def __init__(self,Radiative,Dipoles):
        self.Radiative=Radiative
        self.Dipoles=Dipoles

    def __str__(self):
        Name="Dipole Structures for "+self.Radiative+"\n"
        for Dipole in self.Dipoles:
            Name+=str(Dipole)
        return Name

    def Write(self):
        pass
