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
	  
! YuP: Not used, commented out, but still keep it, in case they are needed in future:
!     A(f)=-u_par*df/ds (at cst i)
!      wpa(i,j,k,l)=wpweigt(i,j,l)*(fsmid(i,j,k,l) - &
!       fsmid(ilpm1ef(i,l,-1),j,k,lpm1eff(l,-1)))
!      wpap(i,j,k,l)=cvmgt(wpweigp(i,j,l) * &
!       (fnp1(ilpm1ef(i,l,+1),j,k,lpm1eff(l,+1))-fnp1(i,j,k,l)), &
!       zero,ilpm1ef(i,l,+1).ne.-999)
!      wpacen(i,j,k,l)=cvmgt(wpweigt(i,j,l) * &
!       (fnp1(ilpm1ef(i,l,+1),j,k,lpm1eff(l,+1)) - &
!       fnp1(ilpm1ef(i,l,-1),j,k,lpm1eff(l,-1))) , &
!       zero,ilpm1ef(i,l,+1).ne.-999 .and. ilpm1ef(i,l,-1).ne.-999)

!     contribution to matrix for f(s-1), f(s) and f(s+1)
      wprhsmd(i,j,k,l)=fnhalf(i,j,k,l)/dtreff+velsou(i,j,k,l)
      wprhs(i,j,k,l)=fnhalf(i,j,k,l)/dtreff + &
	0.5*(velsou(i,j,k,l)+velsou(ilpm1ef(i,l,+1),j,k,lpm1eff(l,+1)))

!..................................................................
!     End of statement functions used for parallel transport
!..................................................................
