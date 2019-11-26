#procedure Traces4(N)

.sort
#message >> traces
dimension 4;

#do i = 0, `N'
    b g_;
    .sort
    Keep Brackets;
*    #do j = 10,19
*        id g_(?a, v1?, ?b, v1?, ?c) = g_(?a, vhv`j', ?b, vhv`j', ?c);
*    #enddo
    trace4, `i';
#enddo

#endprocedure
