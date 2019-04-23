module coefmidt_mod

  !---BEGIN USE

  use bcast_mod, only : bcast

  !---END USE

!
!

contains

      subroutine coefmidt(c,nn)
      use param_mod
      use comm_mod
      implicit integer (i-n), real*8 (a-h,o-z)

!..................................................................
!     This routine redefines the theta coefficients at the half
!     mesh points (i+1/2). nn=1 means dd; nn=2 means de;
!     nn=3 means df.
!..................................................................


      dimension c(0:iy,jx)

!.......................................................................

      do 2 j=1,jx
        do 21 i=1,iy-1
          temp1(i,j)=(c(i+1,j)+c(i,j))*.5
 21     continue
 2    continue

!..................................................................
!     set coefficients to zero near pi and 0 to force zero flux there.
!..................................................................

      do 3 j=1,jx
        temp1(iy,j)=0.
!        c(0,j)=0.    Changed 8/19/94
        temp1(0,j)=0.
!$$$c***************************************TRY 090826
!$$$c Object is to ensure no flux at theta=0 by
!$$$c setting flux at both -dy/2 and +dy/2 =0, and similarly
!$$$c around theta=pi.
!$$$        temp1(iy-1,j)=0.
!$$$        temp1(1,j)=0.
!$$$c ==> Unexpected result, of irregularity at theta=0,pi
!$$$c***************************************TRY 090826
 3    continue

!..................................................................
!     redefine fluxes near pass/trapped boundary
!..................................................................

      if (cqlpmod .ne. "enabled") then
        xx=-1.
        if (nn .eq. 3) xx=1
        do 4 j=1,jx
          temp1(itl-1,j)=c(itl-1,j)
          temp1(itl,j)=c(itl+1,j)
          temp1(itu-1,j)=xx*temp1(itl,j)
          temp1(itu,j)=c(itu+1,j)
 4      continue
      endif

!..................................................................
!     set flux=0. at v=0.
!..................................................................

!**Changed 8/19/94, for consistency (bh):      do 6 i=1,iy
      do 6 i=0,iy
        temp1(i,1)=0.
        do 61 j=1,jx
          c(i,j)=temp1(i,j)
 61     continue
 6    continue

!..................................................................
!     force 0 theta flux at pi/2.
!..................................................................

      if (symtrap .eq. "enabled") then
        do 8 j=1,jx
          c(iyh,j)=0.
 8      continue
      endif

      if (nn .ne. 3) return

!     minimum of |F| is set to 1.e-40
      do 7 j=1,jx
        do 71 i=0,iy
          if(abs(c(i,j)) .le. em40) then
            c(i,j)= sign(em40,c(i,j))
          endif
 71     continue
 7    continue

        call bcast(temp1(0,0),zero,iyjx2)
      return
      end
end module coefmidt_mod
