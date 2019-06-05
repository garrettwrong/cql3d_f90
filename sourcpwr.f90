module sourcpwr_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast

  !---END USE

!
!

contains

      subroutine sourcpwr(k)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
!MPIINSERT_INCLUDE


!..................................................................
!     This routine computes the Neutral Beam power in Watts/cc
!..................................................................

      call bcast(tam1,zero,jx)
      do 501 i=1,iy
        do 502 j=1,jx
         if(gone(i,j,k,lr_).eq.zero)then !YuP[2017-11-21] added.
          !YuP: added gone() check: ONLY COUNT not-lost particles.
          !YuP: This is similar to the CQL3D-FOW version (added on 03/23/2015)
          tam1(j)=tam1(j)+source(i,j,k,indxlr_)*cynt2(i,l_)*vptb(i,lr_)
          !Note: tam1 is needed for sorpw_nbi below, which is only used
          !for plots and NetCDF file.
          !It has no effect on the source operator or solution.
          !The consequence of adding if(gone..) condition is
          !the drop of NBI power profile at plasma edge
          !(if lossmode is not disabled) because of banana+Larm.rad losses;
          !see tdpltjop, plot 'FSA SOURCE POWER DEN: ...';
          !Also, 'Power from NBI', and 'Power integrated over radius'
          !become smaller (values are printed in *.ps file).
         endif
 502    continue
 501  continue
      s=0.
      do 503 j=1,jx
        s=s+tam1(j)*tcsgm1(j)*cint2(j)
 503  continue

!..................................................................
!     Source power of species k in Watts/cc averaged over the flux surfa
!..................................................................

      sorpw_nbi(k,lr_)=s*fions(k)*one_*1.6022e-16*zmaxpsii(lr_) !-YuP->NBI source

!MPIINSERT_IF_RANK_EQ_0
      WRITE(*,'(a,2i5,e12.4)') &
            'sourcpwr: k,lr, sorpw_nbi', &
                       k,lr_,sorpw_nbi(k,lr_)
      WRITE(*,*)'----------------------------------------------------'
!MPIINSERT_ENDIF_RANK

      return
      end subroutine sourcpwr
      
      
end module sourcpwr_mod
