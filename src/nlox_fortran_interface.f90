      interface
         subroutine NLOX_OLP_Start(fname, ierr)
           character(len=*) :: fname
           integer :: ierr
         end subroutine NLOX_OLP_Start
      end interface
      
      interface
         subroutine NLOX_OLP_EvalSubProcess_All(i, pp, next, mu, rval, acc)
           integer :: i, next
           double precision, dimension(next) :: pp
           double precision :: mu
           double precision, dimension(4) :: rval
           double precision :: acc
         end subroutine NLOX_OLP_EvalSubProcess_All
      end interface
      
      interface
         subroutine NLOX_OLP_EvalSubProcess(i, typ, ltyp, cp, lcp, pp, next, mu, rval2, acc)
           integer :: i, next, ltyp, lcp
           character(len=*) :: typ, cp
           double precision, dimension(next) :: pp
           double precision :: mu
           double precision, dimension(3) :: rval2
           double precision :: acc
         end subroutine NLOX_OLP_EvalSubProcess
      end interface
      
      interface
         subroutine NLOX_OLP_EvalSubProcess_CC(i, typ, ltyp, cp, lcp, pp, next, mu, rval2, acc)
           integer :: i, next, ltyp, lcp
           character(len=*) :: typ
           character(len=*) :: cp
           double precision, dimension(next) :: pp
           double precision :: mu
           double precision, dimension(3) :: rval2
           double precision :: acc
         end subroutine NLOX_OLP_EvalSubProcess_CC
      end interface
      
      interface
         subroutine NLOX_OLP_SetMass(para, lpar, mass, width, ierr)
           integer :: lpar, ierr
           character(len=*) :: para
           double precision :: mass, width
         end subroutine NLOX_OLP_SetMass
      end interface
      
      interface
         subroutine NLOX_OLP_SetParameter(para, lpar, re, im, ierr)
           integer :: lpar, ierr
           character(len=*) :: para
           double precision :: re, im
         end subroutine NLOX_OLP_SetParameter
      end interface
