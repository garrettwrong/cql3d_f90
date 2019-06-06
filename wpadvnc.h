!     wpadvnc.h
!***********************************************************************
!***********************************************************************

!..................................................................
!     wpadvnc contains the statement functions utilized for the parallel
!     transport implicit time advancement routines and its diagnostics
!..................................................................

!     coefficient u_parallel/ds
      wpweigt(i,j,l)=-vnorm*x(j)*coss(i,l)/dsz(l)
      wpweigp(i,j,l)=-vnorm*x(j)*eszp5(l)*coss(i,l)
!     f_s+1/2 (at cst i)
      fsmid(i,j,k,l)=fnp1(i,j,k,l)*dls(i,j,k,l) + &
      fnp1(ilpm1ef(i,l,+1),j,k,lpm1eff(l,+1))*(1-dls(i,j,k,l))

!     contribution to matrix for f(s-1), f(s) and f(s+1)
      wprhsmd(i,j,k,l)=fnhalf(i,j,k,l)/dtreff+velsou(i,j,k,l)
      wprhs(i,j,k,l)=fnhalf(i,j,k,l)/dtreff + &
	0.5*(velsou(i,j,k,l)+velsou(ilpm1ef(i,l,+1),j,k,lpm1eff(l,+1)))

!..................................................................
!     End of statement functions used for parallel transport
!..................................................................
