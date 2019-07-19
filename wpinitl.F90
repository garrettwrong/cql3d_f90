module wpinitl_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use bcast_mod, only : ibcast
  use r8subs_mod, only : dcopy
  use wpwrng_mod, only : wpwrng

  !---END USE

!
!

contains

      subroutine wpinitl
      use param_mod
      use cqlcomm_mod
      use r8subs_mod, only : cvmgt, dcopy

      implicit integer (i-n), real(c_double) (a-h,o-z)

!..............................................................
!     Initialize some arrays for CQLP case (transport along B).
!     Computes the (Chang-Cooper) weights.
!..............................................................


!.......................................................................
!     1. Weights for defining f(s+1/2)= f(s+1)*(1-dls(s)) + f(s)*dls(s)
!.......................................................................

      do 100 k=1,ngen
        do 110 j=1,jx
          do 120 l=0,setup0%ls+1
            do 130 i=1,iy_(setup0%ls)
              dls(i,j,k,l)=.5
 130        continue
 120      continue
 110    continue
 100  continue

!.......................................................................
!     Ensures correct differentiation at end of intervals
!.......................................................................

      if (sbdry .ne. "periodic") then
        do 140 k=1,ngen
          do 141 j=1,jx
            do 142 i=1,iy_(setup0%ls)
              dls(i,j,k,0)=0.
              dls(i,j,k,setup0%ls)=1.
 142        continue
 141      continue
 140    continue
      endif

!.......................................................................
!     2. If the s mesh is assumed periodic, the numeration of the mesh
!     points s(l) is changed, before solving the matrix problem, in
!     order to keep a band-type matrix. The new points are indexed
!     as follows:
!     1, 2, setup0%ls, 3, setup0%ls-1, ..., setup0%ls/2, setup0%ls-(setup0%ls/2-2), setup0%ls/2+1
!     Thus, the band-width is only twice the natural width. The following
!     arrays define the above transformation between the periodic and the
!     bounded mesh and vice versa (like points anti-clockwise on circle):
!     lsbtopr(1:setup0%ls)= 1, 2, 4, 6, .., setup0%ls, setup0%ls-1, setup0%ls-3, .., 5, 3
!     lsprtob(1:setup0%ls)= 1, 2, setup0%ls, 3, setup0%ls-1, ..., setup0%ls/2, setup0%ls-(setup0%ls/2-2), setup0%ls/2+1
!
!     If 10<nummods<21, then solve a periodic system with 2(setup0%ls-1) points
!.......................................................................

      if (sbdry .eq. "periodic") then
        do 200 l=2,setup0%ls/2
          lsbtopr(l)=(l-1)*2
          lsbtopr(setup0%ls+2-l)=lsbtopr(l)+1
          lsprtob((l-1)*2)=l
          lsprtob((l-1)*2+1)=setup0%ls+2-l
 200    continue
        lsbtopr(1)=1
        lsbtopr(setup0%ls/2+1)=setup0%ls
        lsprtob(1)=1
        lsprtob(setup0%ls)=setup0%ls/2+1
!     extra points: s(0)=s(setup0%ls) and s(setup0%ls+1)=s(1)
        lsbtopr(0)=3
        lsbtopr(setup0%ls+1)=1
      else
!     bounded case: one-to-one correspondance
        do 220 l=0,setup0%ls+1
          lsbtopr(l)=l
          lsprtob(l)=l
 220    continue
      endif

!.......................................................................
!l    3. Define lower and upper limit for l mesh at each value of theta
!.......................................................................

!.......................................................................
!l    3.1 Define l_lower(i); the lowest index l for which pitch angle
!     index i is defined on the local transport pitch angle mesh.
!.......................................................................

      call ibcast(l_lower,1,iymax)
      if (iy_(1).ne.iy_(setup0%ls) .or. iy_(1).ne.iy_(setup0%ls/2)) then
!     cannot have periodic mesh and not full l mesh
!%OS  if (sbdry .eq. "periodic") call wpwrng(2)
        illend=lrors
        if (sbdry .eq. "periodic") illend=lrors/2+1
        do 310 l=2,illend
          do 311 i=iyh_(l-1)+1,iyh_(l)
            l_lower(i)=l
            l_lower(iymax+1-i)=l
 311      continue
 310    continue
      endif

!.......................................................................
!l    3.2 Define l_upper(i); the highest index l for which pitch angle
!     index i is defined on the local transport pitch angle mesh.
!.......................................................................

      call ibcast(l_upper,setup0%ls,iymax)
      if (iy_(1).ne.iy_(setup0%ls) .or. iy_(1).ne.iy_(setup0%ls/2)) then
        illend=lrors
        if (sbdry .eq. "periodic") illend=lrors/2+1
        do 320 l=illend-1,1,-1
          do 321 i=iyh_(l+1)+1,iyh_(l)
            l_upper(i)=l
            l_upper(iymax+1-i)=l
 321      continue
 320    continue
      endif

!.......................................................................
!l    3.3 Check values of l_lower and l_upper
!     So far one assumes that l=1 is where all the particles pass (min. B)
!     => l_lower(i)=1 for all i
!     For l_upper, one should have i=iyh_(l) at l=l_upper(i), i>itl_(1)
!     and l_upper(i)=setup0%ls for i<=itl_(1)
!.......................................................................

      if (iymax .ne. iy_(1)) call wpwrng(11)
      do 330 i=1,iymax
        if (l_lower(i) .ne. 1) call wpwrng(12)
 330  continue

      do 335 i=1,itl_(1)
        if (l_upper(i) .ne. setup0%ls) call wpwrng(13)
        if (l_upper(iymax+1-i) .ne. setup0%ls) call wpwrng(13)
 335  continue
      do 336 i=itl_(1)+1,iyh_(1)
        if (l_upper(min(i+1,iyh_(1))) .ne. l_upper(i)) then
          if (i .ne. iyh_(l_upper(i))) call wpwrng(14)
          if (iy_(l_upper(i))+1-i .ne. iyh_(l_upper(iymax+1-i))+1) &
            call wpwrng(14)
        endif
 336  continue

!.......................................................................
!l    4. Initialize some variables
!.......................................................................
!l    4.1 fnp1 (which may be used before first time-step in wparsou)
!.......................................................................

      call dcopy(iyjx2*ngen*lrors,f(0:iy+1,0:jx+1,1:ngen,1:lrors),1, &
                               fnp0(0:iy+1,0:jx+1,1:ngen,1:lrors),1)
      call dcopy(iyjx2*ngen*lrors,f(0:iy+1,0:jx+1,1:ngen,1:lrors),1, &
                               fnp1(0:iy+1,0:jx+1,1:ngen,1:lrors),1)
      if (sbdry .eq. "periodic") then
        call dcopy(iyjx2*ngen, fnp0(0:iy+1,0:jx+1,1:ngen,1      ),1, &
                               fnp0(0:iy+1,0:jx+1,1:ngen,lrors+1),1)
        call dcopy(iyjx2*ngen, fnp0(0:iy+1,0:jx+1,1:ngen,lrors  ),1, &
                               fnp0(0:iy+1,0:jx+1,1:ngen,0      ),1)
        call dcopy(iyjx2*ngen, fnp1(0:iy+1,0:jx+1,1:ngen,1      ),1, &
                               fnp1(0:iy+1,0:jx+1,1:ngen,lrors+1),1)
        call dcopy(iyjx2*ngen, fnp1(0:iy+1,0:jx+1,1:ngen,lrors  ),1, &
                               fnp1(0:iy+1,0:jx+1,1:ngen,0      ),1)
      endif

!.......................................................................
!l    4.2 electrostatic electric field
!.......................................................................

      call bcast(elparnw(0:setup0%ls+1),elpar0,setup0%ls+2)

!.......................................................................
!l    4.3 Effective value of l-1,l,l+1 for a given l (ex:setup0%ls+1->1 if period.)
!.......................................................................

      do 430 l=2,setup0%ls-1
        lpm1eff(l,-1)=l-1
        lpm1eff(l, 0)=l
        lpm1eff(l,+1)=l+1
 430  continue
      if (sbdry .eq. "periodic") then
        lpm1eff(0,-1)=setup0%ls-1
        lpm1eff(0, 0)=setup0%ls
        lpm1eff(0,+1)=1
        lpm1eff(1,-1)=setup0%ls
        lpm1eff(1, 0)=1
        lpm1eff(1,+1)=2
        lpm1eff(setup0%ls,-1)=setup0%ls-1
        lpm1eff(setup0%ls, 0)=setup0%ls
        lpm1eff(setup0%ls,+1)=1
        lpm1eff(setup0%ls+1,-1)=setup0%ls
        lpm1eff(setup0%ls+1, 0)=1
        lpm1eff(setup0%ls+1,+1)=2
      else
!     assumes 0->1 and setup0%ls+1->setup0%ls
        lpm1eff(0,-1)=1
        lpm1eff(0, 0)=1
        lpm1eff(0,+1)=2
        lpm1eff(1,-1)=1
        lpm1eff(1, 0)=1
        lpm1eff(1,+1)=2
        lpm1eff(setup0%ls,-1)=setup0%ls-1
        lpm1eff(setup0%ls, 0)=setup0%ls
        lpm1eff(setup0%ls,+1)=setup0%ls
        lpm1eff(setup0%ls+1,-1)=setup0%ls-1
        lpm1eff(setup0%ls+1, 0)=setup0%ls
        lpm1eff(setup0%ls+1,+1)=setup0%ls
      endif

!.......................................................................
!l    4.4 Value ieff of theta at l=+-1 such that one stays the same i=cst line.
!     Mainly for fixed_mu case when following the i=cst line for i>iyh.
!     Set value to -999 if mesh point does not exist.
!.......................................................................

      dp999=-999
      do 440 l=0,setup0%ls+1
        do 441 i=1,iyh_(lpm1eff(l,0))
          dpi=i
          ilpm1ef(i,l,-1)=cvmgt(dpi,-dp999,i.le.iyh_(lpm1eff(l,-1)))
          ilpm1ef(i,l, 0)=i
          ilpm1ef(i,l,+1)=cvmgt(dpi,-dp999,i.le.iyh_(lpm1eff(l,+1)))
          ii=iy_(lpm1eff(l,0))+1-i
          dpi=iy_(lpm1eff(l,-1))+1-i
          ilpm1ef(ii,l,-1)=cvmgt(dpi,-dp999, &
            i.le.iyh_(lpm1eff(l,-1)))
          ilpm1ef(ii,l, 0)=ii
          dpi=iy_(lpm1eff(l,+1))+1-i
          ilpm1ef(ii,l,+1)=cvmgt(dpi,-dp999, &
            i.le.iyh_(lpm1eff(l,+1)))
 441    continue

!     i=0 and iy_(l)+1: set to i=1 and iy_(l) respectively such that df/dtheta
!     gives derivative at i=1/2 and iy_(l)-1/2 respectively
        ilpm1ef(0,l,-1)=ilpm1ef(1,l,-1)
        ilpm1ef(0,l, 0)=1
        ilpm1ef(0,l,+1)=ilpm1ef(1,l,+1)
        ii=iy_(lpm1eff(l,0))+1
        ilpm1ef(ii,l,-1)=ilpm1ef(iy_(lpm1eff(l,0)),l,-1)
        ilpm1ef(ii,l, 0)=iy_(lpm1eff(l,0))
        ilpm1ef(ii,l,+1)=ilpm1ef(iy_(lpm1eff(l,0)),l,+1)
 440  continue

      return
      end subroutine wpinitl


end module wpinitl_mod
