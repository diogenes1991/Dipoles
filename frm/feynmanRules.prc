#procedure feynmanRules()
********************************************************************************
*
* feynman Rules not covered by Python (most standard model vertices/propagators)
*
********************************************************************************

cf DVd;

*#ifdef MASSRENORM
*id SF(n?,K?,m?) = SF(n,K,m);
*#endif

* Label d dimensional indices for fermions in loop for later (by Seth)
* Only needed for certain IR HV schemes
#if `LOOPS' == 1 && `dimScheme' == HVIR
splitarg DV,3;
id DV(v1?,v2?,?args1,-k,?args2) = DVd(v1,v2,?args1,-k,?args2);
id DV(v1?,v2?,?args1,k,?args2) = DVd(v1,v2,?args1,k,?args2);
repeat id DV(v1?,v2?,k1?,k2?,?args,m?) = DV(v1,v2,k1+k2,?args,m);
repeat id DVd(v1?,v2?,k1?,k2?,?args,m?) = DVd(v1,v2,k1+k2,?args,m);
#else
id DV(?args) = DVd(?args); * just treat all propagators as normal
#endif * `dimScheme' == HVIR

id SF(n?,k?,m?) = (G(n,k) + m*G(n))*DS(k,m,1);
id DG(k?,m?) = DS(k,m,1);

b DVd;
.sort
Keep Brackets;

#do i=1,1
    id once DVd(v1?,v2?,k?,0) = 
        -DS(k,0,1)*d_(v1,v2)+`GAUGE'*k(v1)*k(v2)*DS(k,0,2);
    id once DVd(v1?,v2?,k?,m?) =
       -DS(k,m,1)*d_(v1,v2)+`GAUGE'*k(v1)*k(v2)*DS(k,(1-`GAUGE')*m,1)*DS(k,m,1);
    if (count(DVd,1) != 0) redefine i "0";
    b DVd;
    .sort:DVd;
    Keep Brackets;
#enddo

#if `dimScheme' == HVIR
b DV;
.sort
Keep Brackets;

#do i=1,1
    id once DV(v1?,v2?,k?,0) = 
        -DS(k,0,1)*d4d(v1,v2)+`GAUGE'*k(v1)*k(v2)*DS(k,0,2);
    id once DV(v1?,v2?,k?,m?) = 
        -DS(k,m,1)*d4d(v1,v2)+`GAUGE'*k(v1)*k(v2)*DS(k,(1-`GAUGE')*m,1)*DS(k,m,1);
    if (count(DV,1) != 0) redefine i "0";
    b DV;
    .sort:DV;
    Keep Brackets;
#enddo

#endif * `dimScheme' == HVIR

*#if `dimScheme' == CDR || `dimScheme' == HV || `LOOPS' == 0
*    b d4d;
*    .sort
*    Keep Brackets;
*    id d4d(?args) = d_(?args);
*#endif

b VVV;
.sort
Keep Brackets;
*#do i=1,1
*    id once VVV(k1?,v1?,k2?,v2?,k3?,v3?) =
*        (k1(v3)-k2(v3))*d_(v1,v2)+(k2(v1)-k3(v1))*d_(v2,v3)+(k3(v2)-k1(v2))*d_(v1,v3);
*    id once VVV(0,v1?,k2?,v2?,k3?,v3?) =
*        (-k2(v3))*d_(v1,v2)+(k2(v1)-k3(v1))*d_(v2,v3)+(k3(v2))*d_(v1,v3);
*    id once VVV(k1?,v1?,0,v2?,k3?,v3?) =
*        (k1(v3))*d_(v1,v2)+(-k3(v1))*d_(v2,v3)+(k3(v2)-k1(v2))*d_(v1,v3);
*    id once VVV(k1?,v1?,k2?,v2?,0,v3?) =
*        (k1(v3)-k2(v3))*d_(v1,v2)+(k2(v1))*d_(v2,v3)+(-k1(v2))*d_(v1,v3);
*    if (count(VVV,1) != 0) redefine i "0";
*    b VVV;
*    .sort:VVV;
*    Keep Brackets;
*#enddo
#do i=1,1
    id once VVV(k1?,v1?,k2?,v2?,k3?,v3?) =
        -1*((k1(v3)-k2(v3))*d_(v1,v2)+(k2(v1)-k3(v1))*d_(v2,v3)+(k3(v2)-k1(v2))*d_(v1,v3));
    id once VVV(0,v1?,k2?,v2?,k3?,v3?) =
        -1*((-k2(v3))*d_(v1,v2)+(k2(v1)-k3(v1))*d_(v2,v3)+(k3(v2))*d_(v1,v3));
    id once VVV(k1?,v1?,0,v2?,k3?,v3?) =
        -1*((k1(v3))*d_(v1,v2)+(-k3(v1))*d_(v2,v3)+(k3(v2)-k1(v2))*d_(v1,v3));
    id once VVV(k1?,v1?,k2?,v2?,0,v3?) =
        -1*((k1(v3)-k2(v3))*d_(v1,v2)+(k2(v1))*d_(v2,v3)+(-k1(v2))*d_(v1,v3));
    if (count(VVV,1) != 0) redefine i "0";
    b VVV;
    .sort:VVV;
    Keep Brackets;
#enddo

*b d4d;
*.sort
*Keep Brackets;
*repeat;
*id d4d(v1?,v2?)*d4d(v1?,v3?) = d4d(v3,v2);
*id d4d(v1?,v2?)*d4d(v3?,v1?) = d4d(v3,v2);
*id d4d(v1?,v2?)*d4d(v3?,v2?) = d4d(v1,v3);
*endrepeat;
.sort

#endprocedure
