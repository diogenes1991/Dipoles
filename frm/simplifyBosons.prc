#procedure simplifyBosons
    b VB;
    .sort:simplifyBosons VB;
    cf VBc;
    Keep Brackets;

    #if `transversePhotons' == 1
        id VB(fEpsA, p1?,p1?) = 0;
        id VB(fEpsAStar, p1?,p1?) = 0;
    #endif
    #if `transverseGluons' == 1
        id VB(fEpsG, p1?, p1?, ?args) = 0;
        id VB(fEpsGStar, p1?, p1?, ?args) = 0;
    #endif
    #if `transverseZW' == 1
        id VB(fEpsW, p1?, p1?) = 0;
        id VB(fEpsWStar, p1?, p1?) = 0;
        id VB(fEpsZ, p1?, p1?) = 0;
        id VB(fEpsZStar, p1?, p1?) = 0;
    #endif
*   gauge conditions
    #if `gluGaugeConditions' == 1
        id VB(fEpsG,?args) = VBc(fEpsG,?args);
        id VB(fEpsGStar, ?args) = VBc(fEpsGStar, ?args);
        id VB(fEpsA,?args) = VBc(fEpsA,?args);
        id VB(fEpsAStar, ?args) = VBc(fEpsAStar, ?args);

* initial pair
        id VBc(fEpsG, p3?{p1,p2}, p4?{p1,p2}, ?args1)*VBc(fEpsG, v1?, p3?, ?args2)  = 0;
        id VBc(fEpsG, p3?{p1,p2}, p4?{p1,p2}, ?args1)*VBc(fEpsG, p5?, p3?, ?args2)  = 0;
        id VBc(fEpsA, p3?{p1,p2}, p4?{p1,p2}, ?args1)*VBc(fEpsA, v1?, p3?, ?args2)  = 0;
        id VBc(fEpsA, p3?{p1,p2}, p4?{p1,p2}, ?args1)*VBc(fEpsA, p5?, p3?, ?args2)  = 0;
* this pair will be matched in axial scheme by computeSCContractions; remove from consideration
        id VBc(fEpsG, v1?, p3?{p1,p2}, ?args1)*VBc(fEpsG, v2?, p4?{p1,p2}, ?args2) 
          = VB(fEpsG, v1, p3, ?args1)*VB(fEpsG, v2, p4, ?args2);
        id VBc(fEpsA, v1?, p3?{p1,p2}, ?args1)*VBc(fEpsA, v2?, p4?{p1,p2}, ?args2) 
          = VB(fEpsA, v1, p3, ?args1)*VB(fEpsA, v2, p4, ?args2);
        id VBc(fEpsGStar, p3?{p1,p2}, p4?{p1,p2}, ?args1)*VBc(fEpsGStar, v1?, p3?, ?args2)  = 0;
        id VBc(fEpsGStar, p3?{p1,p2}, p4?{p1,p2}, ?args1)*VBc(fEpsGStar, p5?, p3?, ?args2)  = 0;
        id VBc(fEpsAStar, p3?{p1,p2}, p4?{p1,p2}, ?args1)*VBc(fEpsAStar, v1?, p3?, ?args2)  = 0;
        id VBc(fEpsAStar, p3?{p1,p2}, p4?{p1,p2}, ?args1)*VBc(fEpsAStar, p5?, p3?, ?args2)  = 0;
* this pair will be matched in axial scheme by computeSCContractions; remove from consideration
        id VBc(fEpsGStar, v1?, p3?{p1,p2}, ?args1)*VBc(fEpsGStar, v2?, p4?{p1,p2}, ?args2) 
          = VB(fEpsGStar, v1, p3, ?args1)*VB(fEpsGStar, v2, p4, ?args2);
        id VBc(fEpsAStar, v1?, p3?{p1,p2}, ?args1)*VBc(fEpsAStar, v2?, p4?{p1,p2}, ?args2) 
          = VB(fEpsAStar, v1, p3, ?args1)*VB(fEpsAStar, v2, p4, ?args2);

* final pair
        #do i = 3, 5
        #do j = 4, 6
            id VBc(fEpsG, p3?{p`i',p`j'}, p4?{p`i',p`j'}, ?args1)*VBc(fEpsG, v1?, p3?, ?args2)  = 0;
            id VBc(fEpsG, p3?{p`i',p`j'}, p4?{p`i',p`j'}, ?args1)*VBc(fEpsG, p5?, p3?, ?args2)  = 0;
            id VBc(fEpsA, p3?{p`i',p`j'}, p4?{p`i',p`j'}, ?args1)*VBc(fEpsA, v1?, p3?, ?args2)  = 0;
            id VBc(fEpsA, p3?{p`i',p`j'}, p4?{p`i',p`j'}, ?args1)*VBc(fEpsA, p5?, p3?, ?args2)  = 0;
* this pair will be matched in axial scheme by computeSCContractions; remove from consideration
            id VBc(fEpsG, v1?, p3?{p`i',p`j'}, ?args1)*VBc(fEpsG, v2?, p4?{p`i',p`j'}, ?args2) 
              = VB(fEpsG, v1, p3, ?args1)*VB(fEpsG, v2, p4, ?args2);
            id VBc(fEpsA, v1?, p3?{p`i',p`j'}, ?args1)*VBc(fEpsA, v2?, p4?{p`i',p`j'}, ?args2) 
              = VB(fEpsA, v1, p3, ?args1)*VB(fEpsA, v2, p4, ?args2);
            id VBc(fEpsGStar, p3?{p`i',p`j'}, p4?{p`i',p`j'}, ?args1)*VBc(fEpsGStar, v1?, p3?, ?args2)  = 0;
            id VBc(fEpsGStar, p3?{p`i',p`j'}, p4?{p`i',p`j'}, ?args1)*VBc(fEpsGStar, p5?, p3?, ?args2)  = 0;
            id VBc(fEpsAStar, p3?{p`i',p`j'}, p4?{p`i',p`j'}, ?args1)*VBc(fEpsAStar, v1?, p3?, ?args2)  = 0;
            id VBc(fEpsAStar, p3?{p`i',p`j'}, p4?{p`i',p`j'}, ?args1)*VBc(fEpsAStar, p5?, p3?, ?args2)  = 0;
* this pair will be matched in axial scheme by computeSCContractions; remove from consideration
            id VBc(fEpsGStar, v1?, p3?{p`i',p`j'}, ?args1)*VBc(fEpsGStar, v2?, p4?{p`i',p`j'}, ?args2) 
              = VB(fEpsGStar, v1, p3, ?args1)*VB(fEpsGStar, v2, p4, ?args2);
            id VBc(fEpsAStar, v1?, p3?{p`i',p`j'}, ?args1)*VBc(fEpsAStar, v2?, p4?{p`i',p`j'}, ?args2) 
              = VB(fEpsAStar, v1, p3, ?args1)*VB(fEpsAStar, v2, p4, ?args2);
        #enddo
        #enddo

* mixed pair
        id VBc(fEpsG, p3?, p4?, ?args1)*VBc(fEpsGStar, v1?, p3?, ?args2)  = 0;
        id VBc(fEpsG, p3?, p4?, ?args1)*VBc(fEpsGStar, p5?, p3?, ?args2)  = 0;
        id VBc(fEpsGStar, p3?, p4?, ?args1)*VBc(fEpsG, v1?, p3?, ?args2)  = 0;
        id VBc(fEpsGStar, p3?, p4?, ?args1)*VBc(fEpsG, p5?, p3?, ?args2)  = 0;
        id VBc(fEpsA, p3?, p4?, ?args1)*VBc(fEpsAStar, v1?, p3?, ?args2)  = 0;
        id VBc(fEpsA, p3?, p4?, ?args1)*VBc(fEpsAStar, p5?, p3?, ?args2)  = 0;
        id VBc(fEpsAStar, p3?, p4?, ?args1)*VBc(fEpsA, v1?, p3?, ?args2)  = 0;
        id VBc(fEpsAStar, p3?, p4?, ?args1)*VBc(fEpsA, p5?, p3?, ?args2)  = 0;

        id VBc(fEpsG, ?args) = VB(fEpsG, ?args);
        id VBc(fEpsGStar, ?args) = VB(fEpsGStar, ?args);
        id VBc(fEpsA, ?args) = VB(fEpsA, ?args);
        id VBc(fEpsAStar, ?args) = VB(fEpsAStar, ?args);

    #endif
#endprocedure
