#procedure treatDiagram(CONJ)
.sort
dimension d;

* Start by making sure there is something with color in the diagram, else
* color code will fail later.
    mul cOlOne;
    .sort:COlOne mult;

    #call feynmanRules

    $kcount = 0;
    .sort:kk start;
    cf KKt;
    s KKflip;

    b PL;
    .sort:PL bracket;
    Keep Brackets;

    id PL = 1/2*(G - G5);

    b G, G5;
    print;
    .sort:Treat diagram Gs;
    Keep Brackets;

    id G(n1?,v1?)*G5 = G(n1,v1)*G5(n1);
    id G5*G(n1?,v1?) = G5(n1)*G(n1,v1);
    id G5*G(n1?) = G5(n1)*G(n1);
    id G(n1?)*G5 = G(n1)*G5(n1);
    id G(n1?,v1?)*G = G(n1,v1)*G(n1);
    id G*G(n1?,v1?) = G(n1)*G(n1,v1);
    id G*G(n1?) = G(n1)*G(n1);
    id G(n1?)*G = G(n1)*G(n1);

*    id once G(n1?,v1?)*G5(n1?) = - G5(n1)*G(n1,v1);

*    id G(n1?,v1?)*G5(n1?) = 1/2*(G(n1,v1)*G5(n1)-G5(n1)*G(n1,v1));

    b G, G5;
    print;
    .sort:Treat diagram print;

    b G;
    .sort:Treat diagram G bracket;
    Keep Brackets;
    splitarg G;

    repeat id once G(n1?,?args,k1?, k2?) = G(n1,?args,k1) + G(n1,k2);
    id G(n1?,-k1?{k,p1,p2,p3,p4,p5,p6}) = - G(n1,k1);
    .sort:Gamma momentum fixing;


#if (`LOOPS' > 0)
*   extract loop momentum, so we can shift it into tensor integral definition
* Development note: there is an extreme amount of repetition here.
* Try to unify several of these blocks with a pattern for the functions.
    #do i=1,1
        $kcount = `$kcount' + 1;
        b G;
        .sort:Loop momentum extraction;
        Keep Brackets;
        id once G(n?,k) = G(n,vk`$kcount')*KKt(k(vk`$kcount'));
        if(match(G(n1?,k))) redefine i "0";
        .sort:More loop extraction;
    #enddo

    id k1?.k2? = SS(k1,k2);
    #do i=1,1
        $kcount = `$kcount' + 1;
        b SS;
        .sort:Bracket SS;
        Keep Brackets;
        id once SS(k,k1?) = KKt(k(vk`$kcount'))*k1(vk`$kcount');
        if(match(SS(k,k1?))) redefine i "0";
        .sort:Finding SS;
    #enddo

    #do i=1,1
         $kcount = `$kcount' + 1;
         b EpsAStar,EpsA,EpsWStar,EpsW,EpsZStar,EpsZ;
         .sort:Eps Bracket;
         Keep Brackets;
         id once EpsA?{EpsAStar,EpsA,EpsWStar,EpsW,EpsZStar,EpsZ}(k,k1?) = KKt(k(vk`$kcount'))*EpsA(vk`$kcount', k1);
         if(match(EpsA?{EpsAStar,EpsA,EpsWStar,EpsW,EpsZStar,EpsZ}(k,k1?))) redefine i "0";
         .sort:Eps loop mom search;
    #enddo

*    #do i=1,1
*        $kcount = `$kcount' + 1;
*        b EpsGStar;
*        .sort:Bracketing glue;
*        Keep Brackets;
*        id once EpsGStar(k,k1?, cOli1?) = KKt(k(vk`$kcount'))*EpsGStar(vk`$kcount', k1, cOli1);
*        if(match(EpsGStar(k,k1?, cOli1?))) redefine i "0";
*        .sort:Finding glue;
*    #enddo

*    #do i=1,1
*        $kcount = `$kcount' + 1;
*        b EpsG;
*        .sort:Bracketing glue;
*        Keep Brackets;
*        id once EpsG(k,k1?, cOli1?) = KKt(k(vk`$kcount'))*EpsG(vk`$kcount', k1, cOli1);
*        if(match(EpsG(k,k1?, cOli1?))) redefine i "0";
*        .sort:Finding glue;
*    #enddo

    #do i=1,1
        $kcount = `$kcount' + 1;
        b EpsG,EpsGStar;
        .sort:Bracketing glue;
        Keep Brackets;
        id once EpsG?{EpsG,EpsGStar}(k,k1?, cOli1?) = KKt(k(vk`$kcount'))*EpsG(vk`$kcount', k1, cOli1);
        if(match(EpsG?{EpsG,EpsGStar}(k,k1?, cOli1?))) redefine i "0";
        .sort:Finding glue;
    #enddo

    #do i=1,1
        $kcount = `$kcount' + 1;
        b d4d;
        .sort
        Keep Brackets;
*        id d4d(k,v1?) = d4d(v1,k);
        id once d4d(v1?,k) = KKt(k(vk`$kcount'))*d4d(v1,vk`$kcount');
        if(match(d4d(v1?,k))) redefine i "0";
        .sort
    #enddo

    b k;
    .sort:Bracket loop momentum;
    Keep Brackets;

    id k(v1?) = KKt(k(v1));
*   standard form of propagators with q = +k +- ...
    b DS;
    .sort:Bracket DS;
    Keep Brackets;
    splitarg ((k)), DS;
    id DS(-k, n1?, n2?) = DS(k,n1,n2)*KKflip;
    id DS(p1?, -k, n1?, n2?) = DS(-p1, k, n1, n2)*KKflip;
    id DS(p1?!{k},m?,n?) = sDS(p1,m,n);
    repeat id DS(?args, k1?, k2?, n1?, n2?) = DS(?args, k1+k2, n1, n2);
    b KKt, KKflip;
    .sort:Propagator std form;
    Keep Brackets;
    id KKflip^(n1?) = KKflip;
* repeat id KKt(k(v1?))*KKflip = -KK(k(v1))*KKflip;

    id KKflip = 1;
    id KKt(k(v1?)) = KK(k(v1));
*   replacement of non-loop propagators
    b DS, SS;
    .sort:Propagator replace prep;
    Keep Brackets;

*   replace propagators with ints
    id `DSREP' = `TIREP';
    b SS;
    .sort:Int replacement;
    Keep Brackets;

    id SS(k1?,k2?) = k1.k2;

    b KK, TI;
    .sort:TI bracket;

    repeat id KK(k(v1?))*TI(n1?, ?args) = TI(n1, k(v1),?args);
    repeat id TI(n?, ?args, k(v1?)) = TI(n,k(v1),?args);
*    #call tensorDecompose
** CR 20170330: Alternatively there is another Lorentz decomposition routine, to
** be switched with the input-file switch useOldTensorDecomposition (true/false)
** corresponding to oldTensorDec (1/0).
    #if `oldTensorDec' == 1
    #call tensorDecompose
    #elseif `oldTensorDec' == 0
    #include diagsLoopTensorDecompositionModified.frm
    #endif
#endif
.sort:Loop treat diagram;

b G;
.sort:G bracket;
Keep Brackets;

splitarg G;
repeat id G(n1?,?args,k1?, k2?) = G(n1,?args,k1) + G(n1,k2);
repeat id G(n1?,?args,v2?, v3?) = G(n1,?args,v2) + G(n1,v3);
id G(n1?, -k?) = -G(n1, k);
.sort:G momentum fixing;

.sort:Redundant sort?;
dimension d;
#call redefineFields(`CONJ',`LOOPS')

#endprocedure
