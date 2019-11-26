#procedure Traces(N)

.sort:trace presort;

* dimension d;

#message >> traces
#do i = 0, `N'
* deal with free unsummed momenta
    b g_,p1,p2,p3,p4,p5,p6;
    .sort;
    Keep Brackets;

    #do j=1,1
        if ( match(p1?(v1?vhvext$temp)*g_(v2?,?args,v1?,?args2)) );
            sum $temp;
            redefine j "0";
        endif;
        b g_,p1,p2,p3,p4,p5,p6;
        .sort;
        Keep Brackets;
    #enddo

    b g_,d_;
    .sort;
    Keep Brackets;
* simplify 4-d tensors as much as possible, specify HV scheme later
**    repeat; 
**        id d_(v1?vhvext,v2?vhvext)*d_(v1?,v3?) = d_(v2,v3);
**        id d_(v1?vhvext,v2?vhvext)*g_(v3?,?args,v1?,?args2) = g_(v3,?args,v2,?args2);
**        id d_(v1?vhvext,v1?) = 4;
**    endrepeat;
**    b g_,d_;
**    .sort:trace g bracket;
**    Keep Brackets;
    tracen, `i';
#enddo

#endprocedure
