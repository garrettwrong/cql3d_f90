module urfmidv_mod

  !---BEGIN USE

  use bcast_mod, only : bcast
  use r8subs_mod, only : dcopy

  !---END USE

!
!

contains

      subroutine urfmidv_db(jmin,jmax)
      use param_mod
      use comm_mod
      use r8subs_mod, only : dcopy
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save
!..................................................................
!     This routine redefines the velocity mesh points at the
!     half mesh points (j+1/2) so that the fluxes are defined
!     conveniently for differencing. (nn=2 means db,
!     nn=3 means dc).
!     It follows the subroutine coefmidv.
!..................................................................
!..................................................................
!     Special averaging for db; accounts for the
!     fact that near v=0  db is proportional to x**2.
!     Also see loops at the end of this routine.
!..................................................................

      jmax0=min0(jmax,jx-1)
      jminm=max0(jmin-1,1)

        do 110 j=jmin,jmax
          do 120 i=ilim2(j),ilim1(j)
            db(i,j)=db(i,j)/xsq(j)
 120      continue
 110    continue
        call dcopy(iy,db(1:iy,2),1,db(1:iy,1),1)

!     This averaging extends diffusion coeff range to jmin-1:
      do 2 j=jminm,jmax0
        do 3 i=ilim2d(j),ilim1d(j)
          temp1(i,j)=(db(i,j)+db(i,j+1))*.5
 3      continue
 2    continue

      if(jmax.eq.jx) call bcast(temp1(1:iy,jx),zero,iy)

!..................................................................
!     Define db at the pass-trapped boundaries
!     Possibly extends i-range of diffusion coeff. by +/- 1.
!..................................................................
        do 5 j=jminm,jmax
          temp1(itl,j)=.5*( temp1(itl-1,j)/vptb(itl-1,lr_)+ &
            temp1(itl+1,j)/vptb(itl+1,lr_) )*vptb(itl,lr_)
          temp1(itu,j)=.5*( temp1(itu-1,j)/vptb(itu-1,lr_)+ &
            temp1(itu+1,j)/vptb(itu+1,lr_) )*vptb(itu,lr_)
 5      continue

!     ilim2d, ilim1d are extended by 1 in i-direction.
!..................................................................
!     special averaging for db
!..................................................................
      do 6 j=jminm,jmax
        do 7 i=ilim2d(j),ilim1d(j)
          db(i,j)=temp1(i,j)*xcensq(j)
 7      continue
 6    continue

      return
      end








!====================================================================
!====================================================================
      subroutine urfmidv_dc(jmin,jmax)
      use param_mod
      use comm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save
!..................................................................
!     This routine redefines the velocity mesh points at the
!     half mesh points (j+1/2) so that the fluxes are defined
!     conveniently for differencing.
!     It follows the subroutine coefmidv.
!..................................................................
      jmax0=min0(jmax,jx-1)
      jminm=max0(jmin-1,1)
!     This averaging extends diffusion coeff range to jmin-1:
      do 2 j=jminm,jmax0
        do 3 i=ilim2d(j),ilim1d(j)
          temp1(i,j)=(dc(i,j)+dc(i,j+1))*.5
 3      continue
 2    continue

      if(jmax.eq.jx) call bcast(temp1(1:iy,jx),zero,iy)

!..................................................................
!     Set dc=0 at pi and 0
!     Effectively implements bc dF/d(theta) = 0 at pi and 0
!..................................................................
        do 20 j=jminm,jmax
          temp1(1,j)=0.
          temp1(iy,j)=0.
 20     continue

!     ilim2d, ilim1d are extended by 1 in i-direction.
      do 6 j=jminm,jmax
        do 7 i=ilim2d(j),ilim1d(j)
          dc(i,j)=temp1(i,j)
 7      continue
 6    continue

      return
      end
end module urfmidv_mod
