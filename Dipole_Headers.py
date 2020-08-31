def Charge(PID):
    aux = 0
    Qu = 2/3.
    Ql = -1.
    Qw = 1
    if  abs(PID) <= 6 and int(abs(PID))%2==1: aux = Qu-Qw
    if abs(PID) <= 6 and int(abs(PID))%2==0:  aux = Qu
    if abs(PID) == 11 or abs(PID) == 13 or abs(PID) == 15 : aux = Ql
    if abs(PID) == 12 or abs(PID) == 14 or abs(PID) == 16 : aux = Ql+Qw
    if abs(PID) == 24 : aux = Qw
    if  PID < 0 : aux *= -1.0
    return aux

def Symb_Charge(PID):
    aux = ""
    if  abs(PID) <= 6 and int(abs(PID))%2==1: aux = "Qu-Qw"
    if abs(PID) <= 6 and int(abs(PID))%2==0:  aux = "Qu"
    if abs(PID) == 11 or abs(PID) == 13 or abs(PID) == 15 : aux = "Ql"
    if abs(PID) == 12 or abs(PID) == 14 or abs(PID) == 16 : aux = "Ql+Qw"
    if abs(PID) == 24 : aux = "Qw"
    if  PID < 0 : aux = ["-"]+aux 
    return aux
  
def signo( j ):
    if  j >= 0 : return 1.
    else: return -1.

def PID_to_Name(PID):
    out = ""
    aux = abs(PID)
    # Quarks #
    if aux<=6:
        if  aux==0: out = out + "g"
        if  aux==1: out = out + "d"
        if  aux==2: out = out + "u"
        if  aux==3: out = out + "s"
        if  aux==4: out = out + "c"
        if  aux==5: out = out + "b"
        if  aux==6: out = out + "t"
        if  PID<0 : out = out + "bar"
    
    # Leptons #
    elif aux>10 and aux<17:
        if aux%2:
            if aux==11:  out = out + "e"
            if aux==13:  out = out + "mu"
            if aux==15:  out = out + "tau"
            if PID>0  :  out = out + "m"
            else: out = out + "p"
    
        else:
            out = out + "nu"
            if aux==12: out = out + "e"
            if aux==14: out = out + "mu"
            if aux==16: out = out + "tau"
            if PID<0:   out = out + "bar"


    # Gauge Bosons and Higgs #

    elif (aux>20 and  aux<26):
        if aux==22:  out = out + "A"
        if aux==23:  out = out + "Z"
        if PID==24:  out = out + "Wp"
        if PID==-24: out = out + "Wp"
        if aux==25:  out = out + "H"
        
    return out
    

    