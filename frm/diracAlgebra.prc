#procedure diracAlgebra
#define G5Scheme "AG5"

.sort:diracAlgebra presort;
dimension d;
b QQ;
.sort:diracAlgebra QQ bracket;
Keep Brackets;
#message >> diracAlgebra

* Tag quark lines with spin index so we can restore a unit gamma later if
* needed
* Mark the left spinor with the spinline index so we don't forget
* once the identity matrices are removed.
id QQ(n1?,?args)*G5(n2?) = QQ(n1,n2,?args)*G5(n2);
id QQ(n1?,?args)*G(n2?,?args2) = QQ(n1,n2,?args)*G(n2,?args2);
id QQ(n1?,?args)*PL(n2?) = QQ(n1,n2,?args)*PL(n2);
id QQ(n1?,?args)*PR(n2?) = QQ(n1,n2,?args)*PR(n2);

b G, G5;
.sort:diracAlgebra G/G5 bracket;
on properorder;
Keep Brackets;

id G(n1?) = 1;

repeat;
    id G(n1?,v2?)*G(n1?,v2?) = d;
    repeat id G(n1?,p1?)*G(n1?,p1?) = p1.p1;;

* only move G5 past known 4-dim gammas in HV (otherwise unsafe); assume 
* anticommuting in other schemes
    #if `G5Scheme' == HVG5
        repeat id G5(n1?)*G(n1?,v1?ve) = -G(n1,v1)*G5(n1);
    #else
        repeat id G5(n1?)*G(n1?,v1?) = -G(n1,v1)*G5(n1);
    #endif
    repeat id G5(n1?)*G5(n1?) = 1;

*    id once disorder G(n1?,p1?)*G(n1?,p2?) = 2*p1.p2*G(n1) - G(n1,p2)*G(n1,p1);
    id once disorder G(n1?,v1?)*G(n1?,v2?) = 2*d_(v1,v2) - G(n1,v2)*G(n1,v1);

endrepeat;

* Restore unit gamma
b QQ;
.sort:QQ bracket;
Keep Brackets;

id QQ(n1?,n2?,k?,?args)*QQ(?args2) = QQ(n1,k,?args)*G(n2)*QQ(?args2);
id QQ(n1?,n2?,k?,?args) = QQ(n1,k,?args);

#include processSpecific.inc

.sort:diracAlgebra cleanup;

#endprocedure
