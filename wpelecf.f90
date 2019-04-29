module wpelecf_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use iso_c_binding, only : c_double
  use r8subs_mod, only : dcopy

  !---END USE

  real(c_double) :: zdaij
  real(c_double) :: zdaijm1
  real(c_double) :: zddij
  real(c_double) :: zddim1j

contains

      real(c_double) function ghelec(i,j,k)
      use param_mod
      use comm_mod
      use advnce_mod

      ghelec = qz(j)*(zdaij*fpj(i,j,k)-zdaijm1*fpj(i,j-1,k)) &
        + ry(i,j)*(zddij*fpi(i,j)-zddim1j*fpi(i-1,j))
      end function ghelec

!
      subroutine wpelecf(kopt)
      use param_mod
      use comm_mod
      use advnce_mod
      use r8subs_mod, only : dcopy
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..............................................................
!     Computes the parallel component of the poloidal electric field due
!     to charge density along the magnetic field line, assuming toroidal
!     symmetry. The analytical solution of the Poisson equation is:
!     (with E_parallel = E_pol * B_pol/B)
!     E_parallel(s) * B(S)/B_pol(s)**2 = [E_par*B/B_pol**2](s=s_0)
!     + 4*pi* int[s_0,s] (ds'*rho(s')/B(s'))
!     where rho(s')=sum over k of q_k*reden(k,l_) is the local charge density.
!
!     kopt =  1: just compute the new E_parallel
!     2: as kopt=1 and adjust the source term velsou accordingly
!     11: same as 1, but assumes reden is already updated (to f)
!     12: same as 2, but assumes reden is already updated (to f)
!
!     This routine assumes that f is the updated distribution function
!..............................................................

      dimension z4pirho(lsa)
!.......................................................................

      if (n.lt.nonelpr .or. n.gt.noffelpr) return

!.......................................................................
!l    1. Computes the charge density.
!.......................................................................

      do l=1,ls
        z4pirho(l) = 0.0
      end do

!.......................................................................
!l    1.1 Integrate over distribution function if needed
!.......................................................................

      if (kopt .le. 10) then
        do 110 k=1,ngen
          ztra1=charge*bnumb(k)*4.*pi
          do 111 j=1,jx
            ztra2=ztra1*cint2(j)
            do 112 l=1,ls
              do 113 i=1,iy_(l)
                z4pirho(l)=z4pirho(l)+ztra2*f(i,j,k,l)*cynt2(i,l)
 113          continue
 112        continue
 111      continue
 110    continue
      endif

!.......................................................................
!     1.2 Add other species, whose density is already updated
!     nkconro(i) gives the species indices which contribute to the charge
!     density.
!.......................................................................

      do 120 ik=1,ntotal
        k=nkconro(ik)
        if (k .eq. 0) go to 121
        if (k.le.ngen .and. kopt.le.10) go to 120
        ztra1=bnumb(k)*charge*4.*pi
        do 122 l=1,ls
          z4pirho(l)=z4pirho(l)+ztra1*denpar(k,lsindx(l))
 122    continue
 120  continue
 121  continue

!.......................................................................
!l    2. Compute the new parallel component of the electric field
!.......................................................................

      call dcopy(ls+2,elparnw(0:ls+1),1,elparol(0:ls+1),1)
      zsumrho=0.0
      elparnw(1)=elpar0
      zel0cof=elparnw(1)/psipols(1)**2
      do 200 l=2,ls
        zsumrho=zsumrho+0.5*dszm5(l)*(z4pirho(l-1)/psis(l-1)+ &
          z4pirho(l)  /psis(l))
        elparnw(l)=psipols(l)**2/psis(l)*(zel0cof+zsumrho)
 200  continue

      if (sbdry .eq. "periodic") then
        elparnw(0)=elparnw(ls)
        elparnw(ls+1)=elparnw(1)
      else
        elparnw(0)=0.0
        elparnw(ls+1)=0.0
      endif

!%OS  if (kopt.eq.1 .or. kopt.eq.11) return
      return
!%OS

!.......................................................................
!l    3. Adjust velsou according to the new value of E_parallel. This
!     is needed if E_parallel is recomputed at half time-steps, in
!     between the velocity split and the spatial split.
!.......................................................................

      do 300 k=1,ngen
        ztra1=bnumb(k)*charge/fmass(k)/vnorm
        do 310 l=1,ls
          zdepar=ztra1*(elparnw(l)-elparol(l))*0.25
          do 320 j=1,jx
            zdacofp=-zdepar*(x(j)+x(j+1-1/(jx+1-j)))**2
            zdacofm=-zdepar*(x(j-1+1/j)+x(j))**2
            zddcof=zdepar*x(j)**2
!DIR$ NOVECTOR
            do 330 i=1,iy_(l)
              zdaij=zdacofp*coss(i,l)
              zdaijm1=zdacofm*coss(i,l)
              zddij=zddcof*(sinn(i,l)+sinn(i+1-1/(iy_(l)+1-i),l))**2
              zddim1j=zddcof*(sinn(i-1+1/i,l)+sinn(i,l))**2
              velsou(i,j,k,l)=velsou(i,j,k,l)+ghelec(i,j,k)
 330        continue
 320      continue
 310    continue
 300  continue

      return
      end
end module wpelecf_mod
