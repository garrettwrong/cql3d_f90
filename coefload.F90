module coefload_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use losstor_mod, only : losstor

  !---END USE

!
!

contains

      subroutine coefload(k)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     This routine adds in the krook operator contribution to the
!     coefficients employed in time advancement..
!..................................................................



!..................................................................
!     toroidal loss term
!..................................................................

      call losstor(k)

!..................................................................
!     orbit loss term - see subroutine losscone.
!..................................................................

      do 30 i=1,iy
        do 31 j=1,jx
          if(gone(i,j,k,indxlr_).lt.zero) then ! lost
            gon(i,j)=vnorm*x(j)*gone(i,j,k,indxlr_)/tau(i,lr_)
          elseif(gone(i,j,k,indxlr_).eq.zero)then ! confined
            gon(i,j)=zero
          else
            gon(i,j)=(1.-exp(dtreff*gone(i,j,k,indxlr_)))/dtreff
          endif

!..................................................................
!     Add in the contributions from subroutine losstor (taulos)
!     cah when multiplied by vptb(i,lr_) and by F becomes the Krook opera
!..................................................................

          cah(i,j)=gon(i,j)-1./taulos(i,j,indxlr_)
 31     continue
 30   continue

!      if(k.eq.2 .and. sum(gon).ne.zero)then
!      write(*,*)'coefload: k,lr_,sum(gon)',k,lr_,sum(gon)
!      endif

      return
      end subroutine coefload

end module coefload_mod
