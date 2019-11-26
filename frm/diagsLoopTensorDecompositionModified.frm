
** CR
** Taken from diagsLoopTensorDecomposition.frm in
** the amptools2 development branch (Christian2),
** on 20170213 around 2:30pm.

**************************************************
** Loop Tensor Decomposition
**************************************************
** TI stands for Tensor Integral.
** TC stands for Tensor Coefficient.
** Our N-point loop integrals of rank P are represented in the form
**   TI(N,P,v1,...,v{P},q1,...,q{N-1},m0,m1,...,m{N-1}) .
** Our tensor coefficients are represented in the form
**   TC(N,P,i1,i2,...,q1,q2,...,m0,m1,m2,...) .

b KK, TI;
.sort

#message >> tensor decomposition
repeat id KK(k(v1?))*TI(n1?, ?args) = TI(n1, k(v1),?args);
repeat id TI(n?, ?args, k(v1?)) = TI(n,k(v1),?args);
.sort
** CR
*#message CRPRINT loop tensors identified 1
*print +s;
*.sort
** CR

** Declare some needed objects.
auto cf TC;
cf Q;
cf D1,DD,DDD;
cf DD,DDD;
cf DQ,DQQ,DQQQ,DQQQQ;
cf DDQ,DDQQ;
s si,sj,sk,sl,sm,sn;
auto s m;

** Translate the TI tensor integral notation from master branch notation
** to amptools2 development branch notation.
** Add in rank IDs.
id TI(n1?,p1?!{k},?args) = TI(n1,0,p1,?args);
.sort
id TI(n1?,k(v1?),k(v2?),k(v3?),k(v4?),k(v5?),k(v6?),?args) = TI(n1,6,v1,v2,v3,v4,v5,v6,?args);
.sort
id TI(n1?,k(v1?),k(v2?),k(v3?),k(v4?),k(v5?),?args) = TI(n1,5,v1,v2,v3,v4,v5,?args);
.sort
id TI(n1?,k(v1?),k(v2?),k(v3?),k(v4?),?args) = TI(n1,4,v1,v2,v3,v4,?args);
.sort
id TI(n1?,k(v1?),k(v2?),k(v3?),?args) = TI(n1,3,v1,v2,v3,?args);
.sort
id TI(n1?,k(v1?),k(v2?),?args) = TI(n1,2,v1,v2,?args);
.sort
id TI(n1?,k(v1?),?args) = TI(n1,1,v1,?args);
.sort
** Shift momenta (not needed for 2-point and lower).
id TI(6,n1?,?args1,q1?,q2?,q3?,q4?,q5?,m0?,m1?,m2?,m3?,m4?,m5?)
  = TI(6,n1,?args1,q1,q1+q2,q1+q2+q3,q1+q2+q3+q4,q1+q2+q3+q4+q5,m0,m1,m2,m3,m4,m5);
id TI(5,n1?,?args1,q1?,q2?,q3?,q4?,m0?,m1?,m2?,m3?,m4?)
  = TI(5,n1,?args1,q1,q1+q2,q1+q2+q3,q1+q2+q3+q4,m0,m1,m2,m3,m4);
id TI(4,n1?,?args1,q1?,q2?,q3?,m0?,m1?,m2?,m3?)
  = TI(4,n1,?args1,q1,q1+q2,q1+q2+q3,m0,m1,m2,m3);
id TI(3,n1?,?args1,q1?,q2?,m0?,m1?,m2?)
  = TI(3,n1,?args1,q1,q1+q2,m0,m1,m2);
.sort

** CR
*#message CRPRINT loop tensors identified 2
*print +s;
*.sort
** CR

*b TI;
.sort
*keep brackets;

** 1-point is special, and only those of even rank are non-zero.
id TI(1,0,?args) = TC(1,0,0,?args);
id TI(1,2,v1?,v2?,?args) = d_(v1,v2)*TC(1,2,0,0,?args);
id TI(1,4,v1?,v2?,v3?,v4?,?args) = DD(v1,v2,v3,v4)*TC(1,4,0,0,0,0,?args);
id TI(1,6,v1?,v2?,v3?,v4?,v5?,v6?,?args) = DDD(v1,v2,v3,v4,v5,v6)*TC(1,6,0,0,0,0,0,0,?args);
id TI(1,5,v1?,v2?,v3?,v4?,v5?,?args) = 0;
id TI(1,3,v1?,v2?,v3?,?args) = 0;
id TI(1,1,v1?,?args) = 0;

** Now from 2-point on
** Rank 0
id TI(n1?,0,?args) = TC(n1,0,0,?args);
** Rank 1
id TI(n1?,1,v1?,?args) =
+ sum_(si,1,n1-1,
    Q(n1,si,v1,?args)*TC(n1,1,si,?args)
  )
;
** Rank 2
id TI(n1?,2,v1?,v2?,?args) =
+ d_(v1,v2)*TC(n1,2,0,0,?args)
+ sum_(si,1,n1-1,
    sum_(sj,1,n1-1,
      Q(n1,si,v1,?args)*
      Q(n1,sj,v2,?args)*
      TC(n1,2,si,sj,?args)
    )
  )
;
** Rank 3
id TI(n1?,3,v1?,v2?,v3?,?args) = 
+ sum_(si,1,n1-1,
    DQ(n1,si,v1,v2,v3,?args)*
    TC(n1,3,0,0,si,?args)
  )
+ sum_(si,1,n1-1,
    sum_(sj,1,n1-1,
      sum_(sk,1,n1-1,
        Q(n1,si,v1,?args)*
        Q(n1,sj,v2,?args)*
        Q(n1,sk,v3,?args)*
        TC(n1,3,si,sj,sk,?args)
      )
    )
  )
;
** Rank 4
id TI(n1?,4,v1?,v2?,v3?,v4?,?args) = 
+ DD(v1,v2,v3,v4)*TC(n1,4,0,0,0,0,?args)
+ sum_(si,1,n1-1,
    sum_(sj,1,n1-1,
      DQQ(n1,si,sj,v1,v2,v3,v4,?args)*
      TC(n1,4,0,0,si,sj,?args)
    )
  )
+ sum_(si,1,n1-1,
    sum_(sj,1,n1-1,
      sum_(sk,1,n1-1,
        sum_(sl,1,n1-1,
          Q(n1,si,v1,?args)*
          Q(n1,sj,v2,?args)*
          Q(n1,sk,v3,?args)*
          Q(n1,sl,v4,?args)*
          TC(n1,4,si,sj,sk,sl,?args)
        )
      )
    )
  )
;
** Rank 5
id TI(n1?,5,v1?,v2?,v3?,v4?,v5?,?args) = 
+ sum_(si,1,n1-1,
    DDQ(n1,si,v1,v2,v3,v4,v5,?args)*
    TC(n1,5,0,0,0,0,si,?args)
  )
+ sum_(si,1,n1-1,
    sum_(sj,1,n1-1,
      sum_(sk,1,n1-1,
        DQQQ(n1,si,sj,sk,v1,v2,v3,v4,v5,?args)*
        TC(n1,5,0,0,si,sj,sk,?args)
      )
    )
  )
+ sum_(si,1,n1-1,
    sum_(sj,1,n1-1,
      sum_(sk,1,n1-1,
        sum_(sl,1,n1-1,
          sum_(sm,1,n1-1,
            Q(n1,si,v1,?args)*
            Q(n1,sj,v2,?args)*
            Q(n1,sk,v3,?args)*
            Q(n1,sl,v4,?args)*
            Q(n1,sm,v5,?args)*
            TC(n1,5,si,sj,sk,sl,sm,?args)
          )
        )
      )
    )
  )
;
** Rank 6
id TI(n1?,6,v1?,v2?,v3?,v4?,v5?,v6?,?args) = 
+ DDD(v1,v2,v3,v4,v5,v6)*TC(n1,6,0,0,0,0,0,0,?args)
+ sum_(si,1,n1-1,
    sum_(sj,1,n1-1,
      DDQQ(n1,si,sj,v1,v2,v3,v4,v5,v6,?args)*
      TC(n1,6,0,0,0,0,si,sj,?args)
    )
  )
+ sum_(si,1,n1-1,
    sum_(sj,1,n1-1,
      sum_(sk,1,n1-1,
        sum_(sl,1,n1-1,
          DQQQQ(n1,si,sj,sk,sl,v1,v2,v3,v4,v5,v6,?args)*
          TC(n1,6,0,0,si,sj,sk,sl,?args)
        )
      )
    )
  )
+ sum_(si,1,n1-1,
    sum_(sj,1,n1-1,
      sum_(sk,1,n1-1,
        sum_(sl,1,n1-1,
          sum_(sm,1,n1-1,
            sum_(sn,1,n1-1,
              Q(n1,si,v1,?args)*
              Q(n1,sj,v2,?args)*
              Q(n1,sk,v3,?args)*
              Q(n1,sl,v4,?args)*
              Q(n1,sm,v5,?args)*
              Q(n1,sn,v6,?args)*
              TC(n1,6,si,sj,sk,sl,sm,sn,?args)
            )
          )
        )
      )
    )
  )
;
*.sort

*b Q,DQQQQ,DQQQ,DDQQ,DQQ,DDQ,DDD,DD;
.sort
*keep brackets;

id DQQQQ(n1?,si?,sj?,sk?,sl?,v1?,v2?,v3?,v4?,v5?,v6?,?args) =
  (d_(v1,v2)*Q(n1,si,v3,?args)*Q(n1,sj,v4,?args)*Q(n1,sk,v5,?args)*Q(n1,sl,v6,?args)
  +d_(v1,v3)*Q(n1,si,v2,?args)*Q(n1,sj,v4,?args)*Q(n1,sk,v5,?args)*Q(n1,sl,v6,?args) 
  +d_(v1,v4)*Q(n1,si,v2,?args)*Q(n1,sj,v3,?args)*Q(n1,sk,v5,?args)*Q(n1,sl,v6,?args) 
  +d_(v1,v5)*Q(n1,si,v2,?args)*Q(n1,sj,v3,?args)*Q(n1,sk,v4,?args)*Q(n1,sl,v6,?args) 
  +d_(v1,v6)*Q(n1,si,v2,?args)*Q(n1,sj,v3,?args)*Q(n1,sk,v4,?args)*Q(n1,sl,v5,?args) 
  +d_(v2,v3)*Q(n1,si,v1,?args)*Q(n1,sj,v4,?args)*Q(n1,sk,v5,?args)*Q(n1,sl,v6,?args) 
  +d_(v2,v4)*Q(n1,si,v1,?args)*Q(n1,sj,v3,?args)*Q(n1,sk,v5,?args)*Q(n1,sl,v6,?args) 
  +d_(v2,v5)*Q(n1,si,v1,?args)*Q(n1,sj,v3,?args)*Q(n1,sk,v4,?args)*Q(n1,sl,v6,?args) 
  +d_(v2,v6)*Q(n1,si,v1,?args)*Q(n1,sj,v3,?args)*Q(n1,sk,v4,?args)*Q(n1,sl,v5,?args) 
  +d_(v3,v4)*Q(n1,si,v1,?args)*Q(n1,sj,v2,?args)*Q(n1,sk,v5,?args)*Q(n1,sl,v6,?args) 
  +d_(v3,v5)*Q(n1,si,v1,?args)*Q(n1,sj,v2,?args)*Q(n1,sk,v4,?args)*Q(n1,sl,v6,?args) 
  +d_(v3,v6)*Q(n1,si,v1,?args)*Q(n1,sj,v2,?args)*Q(n1,sk,v4,?args)*Q(n1,sl,v5,?args) 
  +d_(v4,v5)*Q(n1,si,v1,?args)*Q(n1,sj,v2,?args)*Q(n1,sk,v3,?args)*Q(n1,sl,v6,?args) 
  +d_(v4,v6)*Q(n1,si,v1,?args)*Q(n1,sj,v2,?args)*Q(n1,sk,v3,?args)*Q(n1,sl,v5,?args) 
  +d_(v5,v6)*Q(n1,si,v1,?args)*Q(n1,sj,v2,?args)*Q(n1,sk,v3,?args)*Q(n1,sl,v4,?args));
id DDQQ(n1?,si?,sj?,v1?,v2?,v3?,v4?,v5?,v6?,?args) =
  (DD(v3,v4,v5,v6)*Q(n1,si,v1,?args)*Q(n1,sj,v2,?args)
  +DD(v2,v4,v5,v6)*Q(n1,si,v1,?args)*Q(n1,sj,v3,?args)
  +DD(v2,v3,v5,v6)*Q(n1,si,v1,?args)*Q(n1,sj,v4,?args)
  +DD(v2,v3,v4,v6)*Q(n1,si,v1,?args)*Q(n1,sj,v5,?args)
  +DD(v2,v3,v4,v5)*Q(n1,si,v1,?args)*Q(n1,sj,v6,?args)
  +DD(v1,v4,v5,v6)*Q(n1,si,v2,?args)*Q(n1,sj,v3,?args)
  +DD(v1,v3,v5,v6)*Q(n1,si,v2,?args)*Q(n1,sj,v4,?args)
  +DD(v1,v3,v4,v6)*Q(n1,si,v2,?args)*Q(n1,sj,v5,?args)
  +DD(v1,v3,v4,v5)*Q(n1,si,v2,?args)*Q(n1,sj,v6,?args)
  +DD(v1,v2,v5,v6)*Q(n1,si,v3,?args)*Q(n1,sj,v4,?args)
  +DD(v1,v2,v4,v6)*Q(n1,si,v3,?args)*Q(n1,sj,v5,?args)
  +DD(v1,v2,v4,v5)*Q(n1,si,v3,?args)*Q(n1,sj,v6,?args)
  +DD(v1,v2,v3,v6)*Q(n1,si,v4,?args)*Q(n1,sj,v5,?args)
  +DD(v1,v2,v3,v5)*Q(n1,si,v4,?args)*Q(n1,sj,v6,?args)
  +DD(v1,v2,v3,v4)*Q(n1,si,v5,?args)*Q(n1,sj,v6,?args));
id DQQQ(n1?,si?,sj?,sk?,v1?,v2?,v3?,v4?,v5?,?args) =
  (d_(v1,v2)*Q(n1,si,v3,?args)*Q(n1,sj,v4,?args)*Q(n1,sk,v5,?args) 
  +d_(v1,v3)*Q(n1,si,v2,?args)*Q(n1,sj,v4,?args)*Q(n1,sk,v5,?args) 
  +d_(v1,v4)*Q(n1,si,v2,?args)*Q(n1,sj,v3,?args)*Q(n1,sk,v5,?args) 
  +d_(v1,v5)*Q(n1,si,v2,?args)*Q(n1,sj,v3,?args)*Q(n1,sk,v4,?args) 
  +d_(v2,v3)*Q(n1,si,v1,?args)*Q(n1,sj,v4,?args)*Q(n1,sk,v5,?args) 
  +d_(v2,v4)*Q(n1,si,v1,?args)*Q(n1,sj,v3,?args)*Q(n1,sk,v5,?args) 
  +d_(v2,v5)*Q(n1,si,v1,?args)*Q(n1,sj,v3,?args)*Q(n1,sk,v4,?args) 
  +d_(v3,v4)*Q(n1,si,v1,?args)*Q(n1,sj,v2,?args)*Q(n1,sk,v5,?args) 
  +d_(v3,v5)*Q(n1,si,v1,?args)*Q(n1,sj,v2,?args)*Q(n1,sk,v4,?args) 
  +d_(v4,v5)*Q(n1,si,v1,?args)*Q(n1,sj,v2,?args)*Q(n1,sk,v3,?args));
id DDQ(n1?,si?,v1?,v2?,v3?,v4?,v5?,?args) =
  (DD(v1,v2,v3,v4)*Q(n1,si,v5,?args)
  +DD(v1,v2,v3,v5)*Q(n1,si,v4,?args)
  +DD(v1,v2,v4,v5)*Q(n1,si,v3,?args)
  +DD(v1,v5,v3,v4)*Q(n1,si,v2,?args)
  +DD(v5,v2,v3,v4)*Q(n1,si,v1,?args));
id DQQ(n1?,si?,sj?,v1?,v2?,v3?,v4?,?args) =
  (d_(v1,v2)*Q(n1,si,v3,?args)*Q(n1,sj,v4,?args)
  +d_(v1,v3)*Q(n1,si,v2,?args)*Q(n1,sj,v4,?args)
  +d_(v1,v4)*Q(n1,si,v2,?args)*Q(n1,sj,v3,?args)
  +d_(v2,v3)*Q(n1,si,v1,?args)*Q(n1,sj,v4,?args)
  +d_(v2,v4)*Q(n1,si,v1,?args)*Q(n1,sj,v3,?args)
  +d_(v3,v4)*Q(n1,si,v1,?args)*Q(n1,sj,v2,?args));
id DQ(n1?,si?,v1?,v2?,v3?,?args) =
  (d_(v1,v2)*Q(n1,si,v3,?args) 
  +d_(v1,v3)*Q(n1,si,v2,?args) 
  +d_(v2,v3)*Q(n1,si,v1,?args));
.sort

** i runs from 1 to pointedness-1, and is the summation index from above.
** The first argument of Q is pointedness.
** For example Q(4,2,vLGY2,-p1,p2,p2-p3,0,0,0,0) should be p2(vLGY2).
#do i = 1, 5
*id Q(6,`i',v1?,q1?,q2?,q3?,q4?,q5?,?args) = q`i'(v1);
*id Q(5,`i',v1?,q1?,q2?,q3?,q4?,?args) = q`i'(v1);
*id Q(4,`i',v1?,q1?,q2?,q3?,?args) = q`i'(v1);
*id Q(3,`i',v1?,q1?,q2?,?args) = q`i'(v1);
*id Q(2,`i',v1?,q1?,?args) = q`i'(v1);
#do t = 2, 6
id Q(`t',`i',v1?,q1?,...,q{`t'-1}?,?args) = q`i'(v1);
#enddo
#enddo
id DDD(?args) = 1/3*distrib_(+1,2,D1,DD,?args);
.sort
id DD(?args) = 1/2*distrib_(+1,2,D1,D1,?args);
.sort
id D1(v1?,v2?) = d_(v1,v2);

*b TC;
.sort
*keep brackets;

** Resolve the permutations of n-point coefficients.
** Rank 0 is trivial (just rename it).
id TC(n1?,0,?args) = TC0(n1,0,?args);
** Rank 1 does not need to be symmetrized (just rename it).
id TC(n1?,1,?args) = TC1(n1,1,?args);
** Rank 2
id TC(n1?,2,?args) = TC2(n1,2,?args);
symmetrize TC2 3,4;
** Rank 3
id TC(n1?,3,?args) = TC3(n1,3,?args);
symmetrize TC3 3,4,5;
** Rank 4
id TC(n1?,4,?args) = TC4(n1,4,?args);
symmetrize TC4 3,4,5,6;
** Rank 5
id TC(n1?,5,?args) = TC5(n1,5,?args);
symmetrize TC5 3,4,5,6,7;
** Rank 6
id TC(n1?,6,?args) = TC6(n1,6,?args);
symmetrize TC6 3,4,5,6,7,8;

*b TC;
.sort
*keep brackets;

** Convert TCs back to TI notation (drop the rank as second argument again)
id TC6(n1?,6,?args) = TI(n1,?args);
id TC5(n1?,5,?args) = TI(n1,?args);
id TC4(n1?,4,?args) = TI(n1,?args);
id TC3(n1?,3,?args) = TI(n1,?args);
id TC2(n1?,2,?args) = TI(n1,?args);
id TC1(n1?,1,?args) = TI(n1,?args);
.sort
** Rank zero tensor coefficients do not carry an index 0 in NLOX notation.
** So, drop also the 0 index for all rank zero coefficients.
id TC0(n1?,0,?args) = TC0(n1,?args);
.sort
id TC0(n1?,0,?args) = TI(n1,?args);
.sort


** Shift momenta back (again not needed for 2-point and lower).
id TI(6,?args1,q1?,q2?,q3?,q4?,q5?,m0?,m1?,m2?,m3?,m4?,m5?)
  = TI(6,?args1,q1,-q1+q2,-q2+q3,-q3+q4,-q4+q5,m0,m1,m2,m3,m4,m5);
id TI(5,?args1,q1?,q2?,q3?,q4?,m0?,m1?,m2?,m3?,m4?)
  = TI(5,?args1,q1,-q1+q2,-q2+q3,-q3+q4,m0,m1,m2,m3,m4);
id TI(4,?args1,q1?,q2?,q3?,m0?,m1?,m2?,m3?)
  = TI(4,?args1,q1,-q1+q2,-q2+q3,m0,m1,m2,m3);
id TI(3,?args1,q1?,q2?,m0?,m1?,m2?)
  = TI(3,?args1,q1,-q1+q2,m0,m1,m2);
.sort

** CR
*#message CRPRINT loop tensors decomposed
*print +s;
*.sort
** CR

