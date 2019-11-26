
*************************
** NLOXSimplifyColor.h **
*************************

** We work strictly with a normalization of T_R=1/2
** The following relations and properties are used for the color reduction in the procedure NLOXSimplifyColorReduction
** (1)
**   Replace adjoint color matrices by traces over fundamental color matrices
**   f(a,b,c)=-i_*2*(Tr(T^a*[T^b,T^c])) 
** (2)
**   Fierz ID in explicit fundamental representation (normalization T_R=1/2)
**   T_{ij}^aT_{kl}^a=1/2*(d_(i,l)*d_(kj)-1/NC*d_(i,j)*d_(k,l))
** (3)
**   Fierz ID in trace form 1 (normalization T_R=1/2)
**   Tr(T^aX)*Tr(T^aY)=1/2*(Tr(XY)-1/NC*Tr(X)*Tr(Y))
** (4)
**   Fierz ID in trace form 2 (normalization T_R=1/2)
**   Tr(T^aXT^aY)=1/2*(Tr(X)*Tr(Y)-1/NC*Tr(XY))
** (5)
**   Fundamental color matrices have normalization T_R=1/2 and are traceless
** (6)
**   The traces of fundamental color matrices are defined over fundamental (N_C-dimensional) color space

** The following procedures are implemented in NLOXSimplifyColor.h
** NLOXSimplifyColorDefineObjects, which
**   defines all the objects (symbols, functions, indices, etc.) that are needed in the procedures NLOXSimplifyColorTranslateIndices, NLOXSimplifyColorApplyIdentities and NLOXSimplifyColorOutFormat
** NLOXSimplifyColorTranslateIndices, which
**   translates the color indices that are initially written to ColorProducts.id.frm to color indices that are more clearly defined in terms of their dimensionality,
**   which enables us to use the Form built-in Kronecker delta and still get the right dimensions during saturated index contraction
** NLOXSimplifyColorApplyIdentities, which
**   carries out the color redcution on the basis of the above given identities
** NLOXSimplifyColorOutFormat, which
**   may be used to define the final output in terms of various color invariants and dimensions
** NLOXSimplifyColorPrintEnd, which
**   prints and ends
** NLOXSimplifyColorHeader, which
**  calls NLOXSimplifyColorDefineObjects
** NLOXSimplifyColorFooter, which
**  calls NLOXSimplifyColorTranslateIndices, NLOXSimplifyColorApplyIdentities, NLOXSimplifyColorOutFormat and NLOXSimplifyColorPrintEnd

** Printing of the output is done outside of the procedures at the end of NLOXSimplifyColor.h, while
** off statistics and #- are assumed to be set outside, before NLOXSimplifyColor.h is included

** Note that asterisks for comments have to start always in the first column


*************************
*************************
*************************

#procedure NLOXSimplifyColorDefineObjects

Symbols NA,NC,CA;

CF cOlT,cOlf;
CF cOlTr(cyclesymmetric);
*CF delta(symmetric);
*Form cannot handle yet argument field wildcards in (anti)symmetric functions
CF delta;
s cOlOne; * Used for colorless fields

CF color;

auto i cOla=NA;
auto i cOli=NC;
auto i cOlj=NC;

* New Dummy adjoint index *
auto i cOlij=NC;

auto i dummyc;

*****

auto i cOlajl=NA,cOlajr=NA,cOlaje=NA;
auto i cOlail=NA,cOlair=NA,cOlaie=NA;

auto i cOljl=NC,cOljr=NC,cOlje=NC;
auto i cOlil=NC,cOlir=NC,cOlie=NC;

** Stolen from declarations.h
** AutoDeclare not possible here, due to the fact that the ids `i' are inserted in between
** However, needed for the translation of these indices into something more sensible
#do i = 0, 30
  i cOli`i'e;
  i cOli`i'l;
  i cOli`i'r;
  i cOlj`i'l;
  i cOlj`i'r;
  i cOlj`i'e;
#enddo

#endprocedure


*************************
*************************
*************************

#procedure NLOXSimplifyColorTranslateIndices

** cOlf() contains adjoint color indices only
** Using argument field wildcards here seems to work
#do i = 0, 30
  id cOlf(?X,cOli`i'e,?Y) = cOlf(?X,cOlaie`i',?Y);
  id cOlf(?X,cOli`i'l,?Y) = cOlf(?X,cOlail`i',?Y);
  id cOlf(?X,cOli`i'r,?Y) = cOlf(?X,cOlair`i',?Y);
  id cOlf(?X,cOlj`i'l,?Y) = cOlf(?X,cOlajl`i',?Y);
  id cOlf(?X,cOlj`i'r,?Y) = cOlf(?X,cOlajr`i',?Y);
#enddo
.sort

** The argument convention for cOlT() is cOlT(i,j,a),
** where i and j are the fundamental (antifundamental) color indices,
** while a is the adjoint color index
** Translate the adjoint color indices
** Out-comment the use of argument field wildcards here
#do i = 0, 30
*  id cOlT(?X,cOli`i'e) = cOlT(?X,cOlaie`i');
*  id cOlT(?X,cOli`i'l) = cOlT(?X,cOlail`i');
*  id cOlT(?X,cOli`i'r) = cOlT(?X,cOlair`i');
*  id cOlT(?X,cOlj`i'l) = cOlT(?X,cOlajl`i');
*  id cOlT(?X,cOlj`i'r) = cOlT(?X,cOlajr`i');
*  id cOlT(?X,cOlj`i'e) = cOlT(?X,cOlaje`i');
  id cOlT(dummyc1?,dummyc2?,cOli`i'e) = cOlT(dummyc1,dummyc2,cOlaie`i');
  id cOlT(dummyc1?,dummyc2?,cOli`i'l) = cOlT(dummyc1,dummyc2,cOlail`i');
  id cOlT(dummyc1?,dummyc2?,cOli`i'r) = cOlT(dummyc1,dummyc2,cOlair`i');
  id cOlT(dummyc1?,dummyc2?,cOlj`i'l) = cOlT(dummyc1,dummyc2,cOlajl`i');
  id cOlT(dummyc1?,dummyc2?,cOlj`i'r) = cOlT(dummyc1,dummyc2,cOlajr`i');
  id cOlT(dummyc1?,dummyc2?,cOlj`i'e) = cOlT(dummyc1,dummyc2,cOlaje`i');
#enddo
.sort
** Translate the fundamental color indices on the right, after the bracket
** Out-comment the use of argument field wildcards here, since they don't make
** sense here
#do i = 0, 30
*  id cOlT(cOli`i'e,?Y) = cOlT(cOlie`i',?Y);
*  id cOlT(cOli`i'l,?Y) = cOlT(cOlil`i',?Y);
*  id cOlT(cOli`i'r,?Y) = cOlT(cOlir`i',?Y);
*  id cOlT(cOlj`i'l,?Y) = cOlT(cOljl`i',?Y);
*  id cOlT(cOlj`i'r,?Y) = cOlT(cOljr`i',?Y);
*  id cOlT(cOlj`i'e,?Y) = cOlT(cOlje`i',?Y);
  id cOlT(cOli`i'e,dummyc1?,dummyc2?) = cOlT(cOlie`i',dummyc1,dummyc2);
  id cOlT(cOli`i'l,dummyc1?,dummyc2?) = cOlT(cOlil`i',dummyc1,dummyc2);
  id cOlT(cOli`i'r,dummyc1?,dummyc2?) = cOlT(cOlir`i',dummyc1,dummyc2);
  id cOlT(cOlj`i'l,dummyc1?,dummyc2?) = cOlT(cOljl`i',dummyc1,dummyc2);
  id cOlT(cOlj`i'r,dummyc1?,dummyc2?) = cOlT(cOljr`i',dummyc1,dummyc2);
  id cOlT(cOlj`i'e,dummyc1?,dummyc2?) = cOlT(cOlje`i',dummyc1,dummyc2);
#enddo
.sort
** Translate the fundamental color indices in the middle
** (needs to be done after translating the fundamental color indices on the
** right, after the bracket)
** Out-comment the use of argument field wildcards here, since they don't make
** sense here (in that case also won't have to adhere to do this step after
** translating the fundamental color indices)
#do i = 0, 30
*  id cOlT(?X,cOli`i'e,?Y) = cOlT(?X,cOlie`i',?Y);
*  id cOlT(?X,cOli`i'l,?Y) = cOlT(?X,cOlil`i',?Y);
*  id cOlT(?X,cOli`i'r,?Y) = cOlT(?X,cOlir`i',?Y);
*  id cOlT(?X,cOlj`i'l,?Y) = cOlT(?X,cOljl`i',?Y);
*  id cOlT(?X,cOlj`i'r,?Y) = cOlT(?X,cOljr`i',?Y);
*  id cOlT(?X,cOlj`i'e,?Y) = cOlT(?X,cOlje`i',?Y);
  id cOlT(dummyc1?,cOli`i'e,dummyc2?) = cOlT(dummyc1,cOlie`i',dummyc2);
  id cOlT(dummyc1?,cOli`i'l,dummyc2?) = cOlT(dummyc1,cOlil`i',dummyc2);
  id cOlT(dummyc1?,cOli`i'r,dummyc2?) = cOlT(dummyc1,cOlir`i',dummyc2);
  id cOlT(dummyc1?,cOlj`i'l,dummyc2?) = cOlT(dummyc1,cOljl`i',dummyc2);
  id cOlT(dummyc1?,cOlj`i'r,dummyc2?) = cOlT(dummyc1,cOljr`i',dummyc2);
  id cOlT(dummyc1?,cOlj`i'e,dummyc2?) = cOlT(dummyc1,cOlje`i',dummyc2);
#enddo
.sort

** The arguments in the delta() function should be fundamental indices only,
** stemming from qqV vertices, with V a non-colored vector boson
** Since it's not unambiguous that the delta functions carry only fundamental
** color indices, we out-comment this step completely and take care of it in
** amptools/scripts/amplitudetools.py:runColorContractions() directly
*#do i = 0, 30
*  id delta(?X,cOli`i'e) = d_(?X,cOlie`i');
*  id delta(?X,cOli`i'l) = d_(?X,cOlil`i');
*  id delta(?X,cOli`i'r) = d_(?X,cOlir`i');
*  id delta(?X,cOlj`i'l) = d_(?X,cOljl`i');
*  id delta(?X,cOlj`i'r) = d_(?X,cOljr`i');
*#enddo
*.sort
*#do i = 0, 30
*  id delta(cOli`i'e,?Y) = d_(cOlie`i',?Y);
*  id delta(cOli`i'l,?Y) = d_(cOlil`i',?Y);
*  id delta(cOli`i'r,?Y) = d_(cOlir`i',?Y);
*  id delta(cOlj`i'l,?Y) = d_(cOljl`i',?Y);
*  id delta(cOlj`i'r,?Y) = d_(cOljr`i',?Y);
*#enddo
*.sort

#endprocedure


*************************
*************************
*************************

#procedure NLOXSimplifyColorApplyIdentities


******************************************************************************************************
** Dio 11192019: These lines are meant to help digest mixed products of Tr(Ta1.Ta2....Tan)*T(i,j,aj)
**               by means of SU(3)'s Fierz Identity, actually I would have guessed that Tr is 
**               unecessary and we could get away with just using colT(<string of adjoint indeces>,i,j).
**               This procedure is inspired by that. There is a hard coded bound of 30 on the maximum
**               number of f's on a single term.

id delta(?X) = d_(?X);

#do i=0,30

    id once cOlf(cOla1?,cOla2?,cOla3?) = -2*i_*
    (cOlT(cOlij{3*'i'},cOlij{3*'i'+1},cOla1)*cOlT(cOlij{3*'i'+1},cOlij{3*'i'+2},cOla2)*cOlT(cOlij{3*'i'+2},cOlij{3*'i'},cOla3)-
    cOlT(cOlij{3*'i'},cOlij{3*'i'+1},cOla1)*cOlT(cOlij{3*'i'+1},cOlij{3*'i'+2},cOla3)*cOlT(cOlij{3*'i'+2},cOlij{3*'i'},cOla2));

repeat;
* Fierz Identity, deletes repeated adjoint indices 
    id cOlT(cOli1?,cOli2?,cOla1?)*cOlT(cOli3?,cOli4?,cOla1?) = (1/2)*(d_(cOli1,cOli4)*d_(cOli2,cOli3) 
                                                                     - (1/NC)*d_(cOli1,cOli2)*d_(cOli3,cOli4));
* Trace Identities * 
    id cOlT(cOli1?,cOli1?,cOla1?) = 0;
    id cOlT(cOli1?,cOli2?,cOla1?)*cOlT(cOli2?,cOli1?,cOla2?) = (1/2)*d_(cOla1,cOla2);
endrepeat;
.sort

#enddo

**
**
** End of Dio's addition.
******************************************************************************************************************

** Operations on color matrices in the adjoint representation first
** (1)
*  id cOlf(cOla1?,cOla2?,cOla3?) = -i_*2*(cOlTr(cOla1,cOla2,cOla3)-cOlTr(cOla1,cOla3,cOla2));
*.sort


*****
*
*** Operations on traces that originate from the adjoint color matrices
*repeat;
*** (3)
*  id cOlTr(cOla1?,?X)*cOlTr(cOla1?,?Y)=1/2*(cOlTr(?X,?Y)-1/NC*cOlTr(?X)*cOlTr(?Y));
*** (4)
*  id cOlTr(cOla1?,?X,cOla1?,?Y)=1/2*(cOlTr(?X)*cOlTr(?Y)-1/NC*cOlTr(?X,?Y));
*** (5) and (6)
*  id cOlTr(cOla1?,cOla2?)=1/2*d_(cOla1,cOla2);
*  id cOlTr(cOla1?)=0;
*  id cOlTr()=NC;
*endrepeat;
*.sort
*
******
*
*** Operations on objects with explicit fundamental color indices
*repeat;
*** (2)
*  id cOlT(cOli1?,cOli2?,cOla1?)*cOlT(cOli3?,cOli4?,cOla1?) = 1/2*(d_(cOli1,cOli4)*d_(cOli3,cOli2)-1/NC*d_(cOli1,cOli2)*d_(cOli3,cOli4));
*** (5)
*  id cOlT(cOli1?,cOli2?,cOla1?)*cOlT(cOli2?,cOli1?,cOla2?) = 1/2*d_(cOla1,cOla2);
*  id cOlT(cOli1?,cOli1?,cOla1?) = 0;
*endrepeat;
*.sort
*
******
*
*** Identify left over traces
*repeat;
*  id cOlT(cOli1?,cOli2?,cOla1?)*cOlT(cOli2?,cOli1?,cOla2?) = cOlTr(cOla1,cOla2);
*  id cOlT(cOli1?,cOli2?,cOla1?)*cOlT(cOli2?,cOli3?,cOla2?)*cOlT(cOli3?,cOli1?,cOla3?) = cOlTr(cOla1,cOla2,cOla3);
*  id cOlT(cOli1?,cOli2?,cOla1?)*cOlT(cOli2?,cOli3?,cOla2?)*cOlT(cOli3?,cOli4?,cOla3?)*cOlT(cOli4?,cOli1?,cOla4?) = cOlTr(cOla1,cOla2,cOla3,cOla4);
*  id cOlT(cOli1?,cOli2?,cOla1?)*cOlT(cOli2?,cOli3?,cOla2?)*cOlT(cOli3?,cOli4?,cOla3?)*cOlT(cOli4?,cOli5?,cOla4?)*cOlT(cOli5?,cOli1?,cOla5?) = cOlTr(cOla1,cOla2,cOla3,cOla4,cOla5);
**  #do i = 0, 30
*  #do i = 5, 30
*    id cOlT(cOli1?,cOli2?,cOla1?)*<cOlT(cOli2?,cOli3?,cOla2?)>*...*<cOlT(cOli`i'?,cOli{`i'+1}?,cOla`i'?)>*cOlT(cOli{`i'+1}?,cOli1?,cOla{`i'+1}?) = cOlTr(cOla1,cOla2,...,cOla`i',cOla{`i'+1});
*  #enddo
*endrepeat;
*.sort
*
******
*
*** Operations on left over traces
*repeat;
*** (3)
*  id cOlTr(cOla1?,?X)*cOlTr(cOla1?,?Y)=1/2*(cOlTr(?X,?Y)-1/NC*cOlTr(?X)*cOlTr(?Y));
*** (4)
*  id cOlTr(cOla1?,?X,cOla1?,?Y)=1/2*(cOlTr(?X)*cOlTr(?Y)-1/NC*cOlTr(?X,?Y));
*** (5) and (6)
*  id cOlTr(cOla1?,cOla2?)=1/2*d_(cOla1,cOla2);
*  id cOlTr(cOla1?)=0;
*  id cOlTr()=NC;
*endrepeat;
*.sort
***
***

#endprocedure


***************
***************
***************

#procedure NLOXSimplifyColorOutFormat

** Express everything in terms of NC
** Put the outcome into a enclosing color() bracket
*id NA=NC^2-1;
*.sort
*AntiBracket NC;
*.sort
*collect color;
*normalize color;
*.sort

*****

** Translate to something that NLOX will definitely understand, i.e. to CA
** Put the outcome into a enclosing color() bracket
*id NA=NC^2-1;
*.sort
*id NC=CA;
*id 1/NC=1/CA;
*.sort
*AntiBracket CA;
*.sort
*collect color;
*normalize color;
*.sort

*****

** Output values directly
id NA=NC^2-1;
.sort
id NC=3;
id 1/NC=1/3;
.sort

#endprocedure


*************************
*************************
*************************

#procedure NLOXSimplifyColorPrintEnd

print;
.end

#endprocedure


*************************
*************************
*************************

#procedure NLOXSimplifyColorHeader

#call NLOXSimplifyColorDefineObjects

#endprocedure


*************************
*************************
*************************

#procedure NLOXSimplifyColorFooter

*#call NLOXSimplifyColorTranslateIndices
** Out-commented the above call, due to complications in it wrt an ambiguity
** in regards to adjoint vs fundamental color indices in the delta functions
** (see further comments in NLOXSimplifyColorTranslateIndices and
** amptools/scripts/amplitudetools.py:runColorContractions())
        
** Now, since the adjoint indices in the delta functions have been taken
** care of, we can take care of the rest of the indices, which should be
** only fundamental ones.
argument;
#do i = 0, 30
    id cOli`i'e = cOlie`i';
    id cOlj`i'e = cOlje`i';
#enddo
endargument;
.sort
id delta(?X) = d_(?X);
id cOlOne = 1;
.sort

** Also, since the new diagram-level color code leaves cOlTr() objects
** with cOli`i'l and possibly cOli`i'r indices, etc. we need to trans-
** late them first: In cOlTr() all indices are adjoint.
argument cOlTr;
#do i = 0, 30
    id cOli`i'l = cOlail`i';
    id cOlj`i'l = cOlajl`i';
    id cOli`i'r = cOlair`i';
    id cOlj`i'r = cOlajr`i';
#enddo
endargument;
.sort

#call NLOXSimplifyColorApplyIdentities
#call NLOXSimplifyColorOutFormat
#call NLOXSimplifyColorPrintEnd

#endprocedure


