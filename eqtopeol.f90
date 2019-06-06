module eqtopeol_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use r8subs_mod, only : dscal
  use zcunix_mod, only : coeff1
  use zcunix_mod, only : terp1

  !---END USE

!
!

contains

      subroutine eqtopeol
      use param_mod
      use cqlcomm_mod
      use r8subs_mod, only : dscal
      implicit integer (i-n), real(c_double) (a-h,o-z)
      integer tnrza
      parameter(tnrza=3*(nnra+nnza)+1)

      dimension br(nnra,nnza),bt(nnra,nnza),bz(nnra,nnza)
      dimension workk(tnrza),temz(nnza),d2bz(nnza),d2br(nnra)
      character*15 char

!..................................................................
!     This routine converts B data read in from CULHAM TOPEOL
!     equilibrium code for use in equilibrium package
!     (Br,Bz)===>psi        R*Bt===>f(psi)
!..................................................................

      ilen=0
      open(unit=10,file='topeol',status='old')

!..................................................................
!     read in topeol data
!     er(1...nnr) contains the major radius points of the grid and
!     ez(1...nnz) the vertical coords of the grid. ez is not to be
!     confused with the z (field line distance) used by mccoy and o'brien.
!     rmincon and rmaxcon are the inner and outer major radii of the
!     LCFS at ez=0.
!     br, bt and bz are the er, toroidal and ez components of the
!     magnetic field as functions of er and ez.
!     Up-down symmetry is assumed.
!..................................................................
!
      read(10,2) char
      read(10,*) nnr,nnz
      write(*,*)
      write(*,*)'eqtopeol: nnr,nnz =', nnr,nnz
      write(*,*)
 3    format(2i5)
      read(10,2) char
      read(10,*) rmincon,rmaxcon
      write(*,*) 'eqtopeol: rmincon,rmaxcon =',rmincon,rmaxcon
      read(10,2) char
      read(10,1) (er(i),i=1,nnr)
!      write(*,*)'eqtopeol:er(i),i=1,nnr =', (er(i),i=1,nnr)
!BH070113 1    format(4(1pe15.8))
 1    format(8(1pe15.8))
      read(10,2) char
      read(10,1) (ez(i),i=1,nnz)
!      write(*,*)'eqtopeol:ez(i),i=1,nnz =', (ez(i),i=1,nnz)
      read(10,2) char
      read(10,1) ((br(i1,j1),i1=1,nnr),j1=1,nnz)
      read(10,2) char
      read(10,1) ((bt(i1,j1),i1=1,nnr),j1=1,nnz)
      read(10,2) char
      read(10,1) ((bz(i1,j1),i1=1,nnr),j1=1,nnz)
!      write(*,*)'eqtopeol: bz =', ((bz(i1,j1),i1=1,nnr),j1=1,nnz)
 2    format(a15)
      nzc=nnz/2+1

!BH   We assume nnz is an odd number, so that ez(nzc) gives
!BH   the equatorial plane.   Checking
      if (abs(ez(nzc)).gt.1.e-9) then
         write(*,*)'eqtopeol:  Problem, ez(nzc).ne.0'
         stop
      endif


!..................................................................
!     Integrate first in R direction at the equatorial plane:
!     This will be the boundary condition for the Z integration.
!     Use splines.
!     Previously (070117) integration was at Z=ez(1).  New method
!     gives calc of psi based entirely on B-fields within the
!     plasma, assuming non-reentrant plasma shape (BobH, 070117).
!..................................................................

      mstep=20
      itab(1)=1
      itab(2)=0
      itab(3)=0
      i1p(1)=4
      i1p(2)=4
      epsi(1,1)=0.
      psivl=0.
      drint=(er(2)-er(1))/mstep
      drint2=drint*.5
      jj=1
      rval=er(1)
      call coeff1(nnr,er,bz(1,nzc),d2bz,i1p,1,workk)
 5    continue
      jj=jj+1

!..................................................................
!     Integrate with mstep subintervals
!..................................................................

      do 10 j=1,mstep
        rval=rval+drint
        rvalm=rval-drint2

!..................................................................
!     Evaluate 1-D spline..
!..................................................................

        call terp1(nnr,er,bz(1,nzc),d2bz,rvalm,1,tab,itab)
        psivl=psivl-drint*rvalm*tab(1)
 10   continue
      epsi(jj,nzc)=psivl
      if (jj.lt.nnr) go to 5
!      write(*,*)'eqtopeol: epsi(i,nnz/2+1) =',(epsi(i,nnz/2+1),i=1,nnr)

!..................................................................
!     Integrate in Z direction: 2 nested loops
!..................................................................

      dzint=(ez(2)-ez(1))/mstep
      dzint2=.5*dzint
      do 30 j=1,nnr
        psivl=epsi(j,nzc)
        ii=nzc
        do 26 i=1,nnz
          temz(i)=br(j,i)
 26     continue
        call coeff1(nnz,ez,temz,d2br,i1p,1,workk)
        zval=ez(nzc)
 25     continue
        ii=ii+1
        do 20 i=1,mstep
          zval=zval+dzint
          zvalm=zval-dzint2
          call terp1(nnz,ez,temz,d2br,zvalm,1,tab,itab)
          psivl=psivl+dzint*er(j)*tab(1)
 20     continue
        epsi(j,ii)=psivl
        if (ii.le.nnz) go to 25
 30   continue

!..................................................................
!     Determine the lower portion by symmetry...
!..................................................................

      do 35 i=1,nzc-1
        ii=nnz+1-i
        do 37 j=1,nnr
          epsi(j,i)=epsi(j,ii)
 37     continue
 35   continue

!..................................................................
!     Force max of psi to be at mag axis (in accord with the convention
!     enforced in equilib.f after reading in the equilibrium data).
!..................................................................

      na=nnza*nnra
      do 40 i=1,nnr
        if (er(i).gt.rmincon) go to 41
 40   continue
 41   continue
      iin=i
      do 42 i= iin+1,nnr
        if (er(i).ge.rmaxcon) go to 43
 42   continue
 43   continue
      iout=i
      amin=ep90
      do 44 i=iin,iout
        val=abs(er(i)*bz(i,nzc))
        if (val.lt.amin) then
          amin=val
          imin=i
        endif
 44   continue
      fac=1.
      if (epsi(imin,nzc).lt. epsi(imin+1,nzc)) fac=-1.


!$$$      write(*,*)'eqtopeol: bz(i:1,nnr,nzc) =', (bz(i1,nzc),i1=1,nnr)
!$$$      write(*,*)'eqtopeol: br(imin,j:1,nnz) =',
!$$$     +                     (br(imin,j1),j1=1,nnz)
!$$$      write(*,*)'eqtopeol:iin,iout,imin,nzc,fac =',iin,iout,imin,nzc,fac

!..................................................................
!     Take toroidal current = 1. (Amp) with dirn given by bz(LCFS)
!..................................................................

      toteqd=sign(one,-bz(iout,nzc))

!..................................................................
!     Convert to cgs from mks
!..................................................................

      hundred=100.d0
      call dscal(nnr,hundred,er,1)
      call dscal(nnz,hundred,ez,1)
      convert=fac*1.e+8
      call dscal(na,convert,epsi,1)
      rbox=er(nnr)-er(1)
      zbox=ez(nnz)-ez(1)
      rboxdst=er(1)
      rmaxcon=rmaxcon*100.
      rmincon=rmincon*100.
      toteqd=toteqd*3.e+9

!..................................................................
!     Find the index imag associated with the maximum value of epsi
!     at ez=0
!..................................................................

      imag=imin

!..................................................................
!     Determine the f(psi) array as required by the "eq" package.
!..................................................................

      jq=0
      do 60 j=iout,imag,-1
        jq=jq+1
        psiar_(jq)=epsi(j,nzc)
        fpsiar_(jq)=1.e+4*bt(j,nzc)*er(j)
 60   continue
      radmaj=er(imag)
      btor=fpsiar_(1)/er(imag)
      ymideqd=0.
      rmag=radmaj
      zmag=0.                 ! Up-down symmetry

!      write(*,*)'eqtopeol: imag,iout,jq =',imag,iout,jq
      write(*,*)'eqtopeol: rbox,zbox,radmaj,rboxdst,ymideqd,btor =', &
                           rbox,zbox,radmaj,rboxdst,ymideqd,btor
!$$$      write(*,*)'eqtopeol: jq, psiar_,i=1,jq =',jq,(psiar_(i),i=1,jq)
!$$$      write(*,*)'eqtopeol: jq, fpsiar_,i=1,jq =',jq,(fpsiar_(i),i=1,jq)
!$$$      write(*,*)'eqtopeol: epsi(imag,i:1,nnz) =',(epsi(imag,i),i=1,nnz)
!$$$      write(*,*)'eqtopeol: epsi(i,nnz/2+1) =',(epsi(i,nnz/2+1),i=1,nnr)

!..................................................................
!     Regrid f(psi) to be specified on nnr equi-spaced psi points
!     from LCFS to magnetic axis (BobH, 070118).
!..................................................................

      psimag_=psiar_(jq)
      psilim_=psiar_(1)
      write(*,*)'eqtopeol,  AT exit: radmaj,psimag,psilim =', &
                                     radmaj,psimag_,psilim_
      dpsi_=(psimag_-psilim_)/(nnr-1)
      do i=1,nnr
         psiar(i)=psilim_+(i-1)*dpsi_
      enddo

      i1p(1)=4
      i1p(2)=4
      call coeff1(jq,psiar_,fpsiar_,d2fpsiar,i1p,1,workk)
      itab(1)=1
      itab(2)=0
      itab(3)=0
      do i=1,nnr
         call terp1(jq,psiar_,fpsiar_,d2fpsiar,psiar(i),1,tab,itab)
         fpsiar(i)=tab(1)
      enddo

!..................................................................
!     Set up spline array for the f =(R*BTOR) subroutine
!..................................................................

      nfp=nnr
      i1p(1)=4
      i1p(2)=4
      call coeff1(nnr,psiar,fpsiar,d2fpsiar,i1p,1,workk)
      return
      end subroutine eqtopeol

end module eqtopeol_mod
