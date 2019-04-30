module diagdens_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast

  !---END USE

!
!

contains

      subroutine diagdens(xline,xmid,eline)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     Compute, for distribution stored in temp4:
!     xline=(line density for complete poloidal turn)/2
!           (particles/cm**2)
!     xmid=density of midplane distribution (particles/cm**3).
!     eline=(complete) line energy density (ergs/cm**2)/(.5*m*vnorm**2)
!..................................................................

      save

      xline=0.
      xmid=0.
      eline=0.
      call bcast(tam2,zero,jx)
      call bcast(tam1,zero,jx)
      do 10 j=1,jx
        do 11 i=1,iy
          tam2(j)=tam2(j)+temp4(i,j)*cynt2(i,l_)*vptb(i,lr_)
          tam1(j)=tam1(j)+temp4(i,j)*cynt2(i,l_)
 11     continue
 10   continue
      do 40 j=1,jx
        xline=xline+tam2(j)*cint2(j)
        xmid=xmid+tam1(j)*cint2(j)
        eline=eline+tam2(j)*cint2(j)*tcsgm1(j)
 40   continue
      return
      end
end module diagdens_mod
