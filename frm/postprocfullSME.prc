#procedure postprocfullSME

b SimpSMEl, SimpSMEr;
.sort
Keep Brackets;
id SimpSMEl(n1?)*SimpSMEr(n2?) = SME(n1, n2);
id SimpSMEl(n1?) = SME(n1, 0);

b DS;
.sort
Keep Brackets;
id DS(p1?,m1?,n1?) = DS(p1,m1,1)^n1;
normalize, (0), DS;
*splitarg DS;
*b DS, p1, p2, p3, p4, p5, p6;
*.sort
*Keep Brackets;
*
*repeat;
*    id DS(p1?, p2?, m?, n?)*p1?.p2? = 1/2*DS(p1,p2,m,n-1) - 1/2*(p1.p1 + p2.p2 + m^2)*DS(p1,p2,m,n);
*    id DS(p1?, ?args, m?, 0) = 1;
*endrepeat;
*
*repeat id DS(p1?,p2?,?args,m?, n?) = DS(p1 + p2, ?args,m, n);
b DS;
.sort
Keep Brackets;
repeat id DS(p1?, m?,n1?)*DS(p1?,m?, n2?) = DS(p1,m,n1+n2);

#include processSpecific.inc;

b p1,p2,p3,p4,p5,p6;
.sort
format mathematica;
Keep Brackets;

b p1,p2,p3,p4,p5,p6;
.sort
Keep Brackets;

id p1?.p2? = SS(p1,p2);

* Expand SMEs in epsilon if they can have any non-const pieces
* The form of fullSME will need to be detected by C code generator to
* determine SME type (Poly3 or complex)
#if (`dimScheme' != HV)
b d;
.sort
Keep Brackets;
id d = 4 - 2*ep;
b ep;
.sort
Keep Brackets;
#do i = {`SMElist'}
    #do j = 0, {`nLOparts'-1}
        l [fullSME({`i'*`nLOparts' + `j'},0)] = [fullSME({`i'*`nLOparts' + `j'})][1];
        l [fullSME({`i'*`nLOparts' + `j'},1)] = [fullSME({`i'*`nLOparts' + `j'})][ep];
        l [fullSME({`i'*`nLOparts' + `j'},2)] = [fullSME({`i'*`nLOparts' + `j'})][ep^2];
    #enddo
#enddo
.sort
#do i = {`SMElist'}
    #do j = 0, {`nLOparts' -1}
        hide [fullSME({`i'*`nLOparts' + `j'})];
    #enddo
#enddo
#endif

*b DS, ep, i_,`bracketconstants',`bracketmasses', e_;
b DS, ep, i_,`bracketconstants', e_;
print;
.end

#endprocedure
