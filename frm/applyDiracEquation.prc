#procedure applyDiracEquation()

***********************************************************************
*
* Application of dirac equation on slashed external momenta sandwiched
* between spinors
*
* Method: a tag, starting from the spinor is commuted between gamma
* matrices until the correct slashed momentum is found. If the tag
* hits a spinor, it is removed. Else, the tag and the gamma matrix are
* commuted back to the spinor and the Dirac equation is applied.
*
***********************************************************************

#define SPINLINES "3"

.sort:applyDirac presort;

cf mass;
s l,r;
f ps;
cf a;
b QQ;
.sort:applyDiracEquation QQ;
on properorder;
Keep Brackets;

* Note: all fields below used to have a required name match to the set "quarks"
* This has been removed for clarity (Dirac eq should apply to all fermions)

* mark fields
#do i=1,{`SPINLINES'*2}
  id once QQ(?args, p1?) = QQ(?args, p1, (-1)^(`i'),0);
#enddo

#if `applyDiracEquation' != 1
    id QQ(?args,p1?, -1,0) = QQ(?args,p1, l);
    id QQ(?args,p1?, 1,0) = QQ(?args,p1, r);
#else
    .sort:Dirac equation;
    #message >> Dirac equation
    id QQ(?args,p1?, -1,0) = QQ(?args,p1, l,0);
    id QQ(?args,p1?, 1,0) = QQ(?args,p1, r,0);

    id QQ(Bar(U),k?, ?args)*G(n1?,k?) = mass(k)*QQ(Bar(U),k,?args)*G(n1);
    id QQ(Bar(V),k?, ?args)*G(n1?,k?) = -mass(k)*QQ(Bar(V),k,?args)*G(n1);
    id G(n1?,k?)*QQ(U,k?, ?args) = mass(k)*G(n1)*QQ(U,k,?args);
    id G(n1?,k?)*QQ(V,k?, ?args) = -mass(k)*G(n1)*QQ(V,k,?args);

* Mark the left spinor with the spinline index so we don't forget in case
* once the identity matrices are removed.
    id QQ(n1?,k?,l,0)*G5(n2?) = QQ(n1,k,n2,l,0)*G5(n2);
    id QQ(n1?,k?,l,0)*G(n2?,?args) = QQ(n1,k,n2,l,0)*G(n2,?args);
    id QQ(n1?,k?,l,0)*PL(n2?) = QQ(n1,k,n2,l,0)*PL(n2);
    id QQ(n1?,k?,l,0)*PR(n2?) = QQ(n1,k,n2,l,0)*PR(n2);

* Remove identity matrices for efficiency, restore later.
* Warning: this will not work if it removes the last gamma matrix in a fermion
* loop. Do not use this procedure until after loop traces have been done.
    id G(n1?) = 1;

*   apply dirac equation to left spinors
    #do j=1,`SPINLINES'
*     #do j=1,4
*       bring all unit matrices in spinor space to the left, get rid
*       of powers of them and absorb them in adjacent gammas
*       *Note*: in this version, we have removed the unit matrices.
*	b G, G5;
*	.sort:Dirac G/G5 bracket;
*	Keep Brackets;
*	repeat;
*            id G(n1?,v1?)*G(n1?) = G(n1)*G(n1,v1);
*            id G5(n1?)*G(n1?) = G(n1)*G5(n1);
*        endrepeat;
*	repeat id G(n1?)*G(n1?) = G(n1);
*	id G(n1?)*G(n1?,v1?) = G(n1,v1);
*	id G(n1?)*G5(n1?) = G5(n1);
*	id G(n1?)*PL(n1?) = PL(n1);
*	id G(n1?)*PR(n1?) = PR(n1);

	b QQ;
	.sort:Dirac QQ bracket;
	Keep Brackets;

        id once QQ(?args,k?,?args2,l,0) = QQ(?args,k,?args2,l,1)*ps(`j', k);
*       if there is a G(n?), spinor chain has no slashed momenta, continue.
*		id ps(`j',K?)*G(n1?) = G(n1);
* If there is nothing between ps and the right QQ, spinor chain has no slashed momenta, continue.
        id ps(`j',k?)*QQ(?args) = QQ(?args);

        #do i=1,1
            .sort:Finished with unit mat;
            if(match(ps(`j', k1?)*G(n1?,k1?)) == 0);
                id once ps(`j',k1?)*G(n1?,k2?) = G(n1,k2)*ps(`j',k1);
            endif;
            id once
                ps(`j',k1?)*G(n1?,v2?!{p1,p2,p3,p4,p5,p6,q1,q2,q3,q4,q5,q6,k}) = G(n1,v2)*ps(`j',k1);
            id once
                ps(`j',k1?)*G5(n1?) = G5(n1)*ps(`j',k1);
            id once
                ps(`j',k1?)*PL(n1?) = PL(n1)*ps(`j',k1);
            id once
                ps(`j',k1?)*PR(n1?) = PR(n1)*ps(`j',k1);
*           arrived to the right, absorb ps
            id once ps(`j', k?)*QQ(?args,r,n2?) = QQ(?args,r,n2);
*           correct p_slash found
            id once ps(`j', k?)*G(n1?,k?) = ps(`j', k, 1)*G(n1,k);
            if(match(ps(`j', k?)) != 0) redefine i "0";
            .sort;some ps matching;
        #enddo

*       move p_slash to the left
        repeat;
*            id once G(n1?,v2?)*ps(`j',p1?,1)*G(n1?,p1?) =  2*d_(p1,v2)*G(n1) - ps(`j', p1,1)*G(n1,p1)*G(n1,v2);
            id once G(n1?,v2?)*ps(`j',p1?,1)*G(n1?,p1?) =  2*d_(p1,v2) - ps(`j', p1,1)*G(n1,p1)*G(n1,v2);
            id once G5(n1?)*ps(`j',p1?,1)*G(n1?,p1?) = - ps(`j',p1,1)*G(n1,p1)*G5(n1);
            id once PL(n1?)*ps(`j',p1?,1)*G(n1?,p1?) = ps(`j',p1,1)*G(n1,p1)*PR(n1);
            id once PR(n1?)*ps(`j',p1?,1)*G(n1?,p1?) = ps(`j',p1,1)*G(n1,p1)*PL(n1);
        endrepeat;

*        id QQ(Bar(U), K?,?args)*ps(`j',K?,1)*G(n1?,K?) = mass(K)*QQ(Bar(U),K, ?args)*G(n1);
*        id QQ(Bar(V), K?,?args)*ps(`j',K?,1)*G(n1?,K?) = -mass(K)*QQ(Bar(V),K, ?args)*G(n1);
        id QQ(Bar(U), k?,?args)*ps(`j',k?,1)*G(n1?,k?) = mass(k)*QQ(Bar(U),k, ?args);
        id QQ(Bar(V), k?,?args)*ps(`j',k?,1)*G(n1?,k?) = -mass(k)*QQ(Bar(V),k, ?args);
        .sort:mass stuff;
    #enddo

*   apply dirac equation to right spinors
*    #do j=1,4
    #do j=1,`SPINLINES'
*       bring all unit matrices in spinor space to the left, get rid
*       of powers of them and absorb them in adjacent gammas
*	b G, G5;
*	.sort:Dirac right spin G/G5;
*        #message right spinors
*	Keep Brackets;
*	repeat;
*            id G(n1?,v1?)*G(n1?) = G(n1)*G(n1,v1);
*            id G5(n1?)*G(n1?) = G(n1)*G5(n1);
*            id PL(n1?)*G(n1?) = G(n1)*PL(n1);
*            id PR(n1?)*G(n1?) = G(n1)*PR(n1);
*        endrepeat;
*	repeat id G(n1?)*G(n1?) = G(n1);
*	id G(n1?)*G(n1?,v1?) = G(n1,v1);
*	id G(n1?)*G5(n1?) = G5(n1);
*	id G(n1?)*PL(n1?) = PL(n1);
*	id G(n1?)*PR(n1?) = PR(n1);

	b QQ;
	.sort:Dirac right QQ;
	Keep Brackets;
        id once QQ(?args2,k?,?args,r,0) = ps(`j', k)*QQ(?args2,k,?args,r,1);
*		id G(n1?)*ps(`j',k?) = G(n1);
        id QQ(?args)*ps(`j',k?) = QQ(?args);
        #do i=1,1
            if(match(G(n1?,k1?)*ps(`j', k1?)) == 0);
                id once G(n1?,k2?)*ps(`j',k1?)= ps(`j',k1)*G(n1,k2);
            endif;
            id once G(n1?,v2?!{p1,p2,p3,p4,p5,p6,q1,q2,q3,q4,q5,q6,k})*ps(`j',k1?) = ps(`j',k1)*G(n1,v2);
            id once G5(n1?)*ps(`j',k1?) = ps(`j',k1)*G5(n1);
            id once PL(n1?)*ps(`j',k1?) = ps(`j',k1)*PL(n1);
            id once PR(n1?)*ps(`j',k1?) = ps(`j',k1)*PR(n1);

            id once QQ(?args, l, n2?)*ps(`j', k?) = QQ(?args,l,n2);

            id once G(n1?,k?)*ps(`j', k?) = G(n1,k)*ps(`j', k, 1);
            if(match(ps(`j', k?)) != 0) redefine i "0";
            .sort
        #enddo

        repeat;
*            id once G(n1?,p1?)*ps(`j',p1?,1)*G(n1?,v2?) = 2*d_(p1,v2)*G(n1) - G(n1,v2)*G(n1,p1)*ps(`j', p1,1);
            id once G(n1?,p1?)*ps(`j',p1?,1)*G(n1?,v2?) = 2*d_(p1,v2) - G(n1,v2)*G(n1,p1)*ps(`j', p1,1);
            id once G(n1?,p1?)*ps(`j',p1?,1)*G5(n1?)= - G5(n1)*G(n1,p1)*ps(`j',p1,1);
            id once G(n1?,p1?)*ps(`j',p1?,1)*PL(n1?)=  PR(n1)*G(n1,p1)*ps(`j',p1,1);
            id once G(n1?,p1?)*ps(`j',p1?,1)*PR(n1?)=  PL(n1)*G(n1,p1)*ps(`j',p1,1);
        endrepeat;

*        id G(n1?,K?)*ps(`j', K?,1)*QQ(U,K?, ?args) = +mass(K)*G(n1)*QQ(U,K,?args);
*        id G(n1?,K?)*ps(`j', K?,1)*QQ(V,K?, ?args) = -mass(K)*G(n1)*QQ(V,K,?args);
        id G(n1?,k?)*ps(`j', k?,1)*QQ(U,k?, ?args) = +mass(k)*QQ(U,k,?args);
        id G(n1?,k?)*ps(`j', k?,1)*QQ(V,k?, ?args) = -mass(k)*QQ(V,k,?args);
        .sort:mass stuff;
    #enddo
    b QQ;
    .sort
    Keep Brackets;
    id QQ(?args, n1?) = QQ(?args);
* Reinsert unit matrix (in G notation) if missing between QQ.
    id QQ(?args, n1?, l)*QQ(?args2, r) = QQ(?args, n1, l)*G(n1)*QQ(?args2, r);
    id QQ(?args, n1?, l) = QQ(?args, l);
    .sort
#endif
#endprocedure
