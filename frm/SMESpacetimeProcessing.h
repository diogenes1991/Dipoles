

off statistics;
#include config.inc
#include `leftDiagDefFile' # global



#procedure ProcessGluons
   
#message >> processing gluon lines

.sort
#define NS "2"

auto v PVEC;
s X;


#do i=0,1


    id VB(fEpsG,v1?,p{`NS'*`i'+1})*VB(fEpsGStar,v2?,p{`NS'*`i'+1})
      =(-d_(v1,v2)+(1/(p{`NS'*`i'+1}.p{`NS'*`i'+2}))*(p{`NS'*`i'+1}(v1)*p{`NS'*`i'+2}(v2)+p{`NS'*`i'+1}(v2)*p{`NS'*`i'+2}(v1)));
    id VB(fEpsG,PVEC?,p{`NS'*`i'+1})*VB(fEpsGStar,v2?,p{`NS'*`i'+1})
      =(-PVEC(v2)+(1/(p{`NS'*`i'+1}.p{`NS'*`i'+2}))*((p{`NS'*`i'+1}.PVEC)*p{`NS'*`i'+2}(v2)+p{`NS'*`i'+1}(v2)*(p{`NS'*`i'+2}.PVEC)));
    id VB(fEpsGStar,PVEC?,p{`NS'*`i'+1})*VB(fEpsG,v2?,p{`NS'*`i'+1})
      =(-PVEC(v2)+(1/(p{`NS'*`i'+1}.p{`NS'*`i'+2}))*((p{`NS'*`i'+1}.PVEC)*p{`NS'*`i'+2}(v2)+p{`NS'*`i'+1}(v2)*(p{`NS'*`i'+2}.PVEC)));
    id VB(fEpsG,PVEC1?,p{`NS'*`i'+1})*VB(fEpsGStar,PVEC2?,p{`NS'*`i'+1})
      =(-(PVEC1.PVEC2)+(1/(p{`NS'*`i'+1}.p{`NS'*`i'+2}))*((p{`NS'*`i'+1}.PVEC1)*(p{`NS'*`i'+2}.PVEC2)+(p{`NS'*`i'+1}.PVEC2)*(p{`NS'*`i'+2}.PVEC1)));
    .sort
    
    id VB(fEpsG,PVEC?,PVEC?) = 0;
    id VB(fEpsGStar,PVEC?,PVEC?) = 0;
    id PVEC?.PVEC? = 0;
    .sort

    id VB(fEpsG,v1?,p{`NS'*`i'+2})*VB(fEpsGStar,v2?,p{`NS'*`i'+2})
      =(-d_(v1,v2)+(1/(p{`NS'*`i'+1}.p{`NS'*`i'+2}))*(p{`NS'*`i'+1}(v1)*p{`NS'*`i'+2}(v2)+p{`NS'*`i'+1}(v2)*p{`NS'*`i'+2}(v1)));
    id VB(fEpsG,PVEC?,p{`NS'*`i'+2})*VB(fEpsGStar,v2?,p{`NS'*`i'+2})
      =(-PVEC(v2)+(1/(p{`NS'*`i'+1}.p{`NS'*`i'+2}))*((p{`NS'*`i'+1}.PVEC)*p{`NS'*`i'+2}(v2)+p{`NS'*`i'+1}(v2)*(p{`NS'*`i'+2}.PVEC)));
    id VB(fEpsGStar,PVEC?,p{`NS'*`i'+2})*VB(fEpsG,v2?,p{`NS'*`i'+2})
      =(-PVEC(v2)+(1/(p{`NS'*`i'+1}.p{`NS'*`i'+2}))*((p{`NS'*`i'+1}.PVEC)*p{`NS'*`i'+2}(v2)+p{`NS'*`i'+1}(v2)*(p{`NS'*`i'+2}.PVEC)));
    id VB(fEpsG,PVEC1?,p{`NS'*`i'+2})*VB(fEpsGStar,PVEC2?,p{`NS'*`i'+2})
      =(-(PVEC1.PVEC2)+(1/(p{`NS'*`i'+1}.p{`NS'*`i'+2}))*((p{`NS'*`i'+1}.PVEC1)*(p{`NS'*`i'+2}.PVEC2)+(p{`NS'*`i'+1}.PVEC2)*(p{`NS'*`i'+2}.PVEC1)));
    .sort 
    
    id VB(fEpsG,PVEC?,PVEC?) = 0;
    id VB(fEpsGStar,PVEC?,PVEC?) = 0;
    id PVEC?.PVEC? = 0;
    .sort
    
    id VB(fEpsG,PVEC1?,p5)*VB(fEpsGStar,PVEC2?,p5)
      =(-(PVEC1.PVEC2));
    .sort
    id VB(fEpsG,PVEC?,PVEC?) = 0;
    id VB(fEpsGStar,PVEC?,PVEC?) = 0;
    id PVEC?.PVEC? = 0;
    .sort
    
    #if (`leftNLoops' == 0) 
    id d = 4;
    #endif

    
#enddo
      
      id PVEC1?.PVEC2?^-4 = 2*DS(PVEC1+PVEC2,0,4);
      id PVEC1?.PVEC2?^-3 = 2*DS(PVEC1+PVEC2,0,3);  
      id PVEC1?.PVEC2?^-2 = 2*DS(PVEC1+PVEC2,0,2);
      id PVEC1?.PVEC2?^-1 = 2*DS(PVEC1+PVEC2,0,1);     
    
#endprocedure