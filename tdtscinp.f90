module tdtscinp_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE


!
!

contains

      subroutine tdtscinp
      use param_mod
      use cqlcomm_mod
      use tdeqdsk_mod, only: psimago, psilimo
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     This routine reads in input from the TSC code output file
!     tscout. This profile information (densities, temperatures
!     then override cqlinput namelist data when there would be
!     a conflict.
!     Magnetic geometry is also read in here. This supersedes eqdsk
!     files - the code will CREATE an eqdsk file later.
!..................................................................


      parameter (nra=33,nza=65)
      dimension epsinew(nra,nza),ernew(nra),eznew(nza)

!..................................................................
!     Open the TSC output file
!..................................................................

      open(unit=11,file='tscinp',status='old')

!..................................................................
!     Description of variables:
!     (note--everything in MKS if not specified)
!
!     names:   80 character description of run from TSC title card
!     npsitm: number of flux surfaces for surface averaged quantities
!     (note npsitm = npsit-1)
!     nspc:   total number of ion species
!     kcycle: TSC time cycle number (starts at 0)
!     times:  TSC problem time in sec
!     anecc:  electron density in particles/cc
!     anicc:       ion density in particles/cc
!     tekev:  electron temp in kev
!     tikev:  ion temp in kev
!     amass:  mass of each ion species in grams
!     achrg:  charge of eqch ion species (hydrogen is 1)
!     elecf: toroidal electric potential in volts
!     rho_:    sqrt of the normalized toroidal flux
!     psiar_:    poloidal flux per radian
!     pary:   pressure (MKS units)
!     ppary:  derivative of pressure wrt xsv
!     fpsiar_:   toroidal field function R*Bt (MKS)
!     gpary:  derivative of fpsiar_ wrt xsv
!     nx:     number of cartesian mesh points in x direction (for psi)
!     nz:     number of cartesian mesh points in z direction (for psi)
!     This normally begins at z=0; so epsi below spans 1/2 of
!     isym:   symmetry option,  0-no symmetry    1- up/down symmetry
!     iplim:  limiter switch   pos-plasma rests on limiter   neg-divert
!     rleft:leftmost boundary of cartesian mesh
!     rright: rightmost boundary of cartesian mesh
!     zbot_:   bottom boundary of cartesian mesh
!     ztop_:    top boundary of cartesian mesh
!     psimag,psilim:  pol flux per radian at mag axis and P/V boundary
!     raxis,zaxis:      (x,z) coordinates of magnetic axis
!     radmaj:    nominal major radius of machine
!     btor:    vacuum field strength at rgzero
!     toteqd:       plasma current in amperes
!     psep,xsep,zsep:  flux value,x and z coordinates of separatrix (2)
!     epsi:       poloidal flux per radian.  Up-down symmetry is
!                 assumed (isym=1). The psi values are given
!                 starting a the midplane (z=0) smallest major radius
!                 point(rleft) to the largest major radius(rright) for
!                 nx points, then back to rleft, z=ztop/(nz-1) to
!                 rright, z=ztop/(nz-1), etc. to top of grid.
!                 Psi values below the midplane are obtained by
!                 reflection through z=0-plane.
!
!..................................................................

      read(11,1001) (names(i),i=1,10)
      read(11,1002) npsitm,nspc,kcycle,times

      if (npsitm.gt.nrada) stop ' in tdtscinp, increase nrada'

      read(11,1003) (anecc(l),l=1,npsitm)
      read(11,1003) (tekev(l),l=1,npsitm)
      read(11,1003) ((anicc(nn,l),l=1,npsitm),nn=1,nspc)
      read(11,1003) ((tikev(nn,l),l=1,npsitm),nn=1,nspc)
      read(11,1003) (amass(nn),nn=1,nspc)
      read(11,1003) (achrg(nn),nn=1,nspc)
      read(11,1003) (elecf(l),l=1,npsitm)
      read(11,1003) (rho_(l),l=1,npsitm)
      read(11,1003) (psiar_(l),l=1,npsitm)
      read(11,1003) (psiar_(l),l=1,npsitm)
      read(11,1003) (pary(l),l=1,npsitm)
      read(11,1003) (ppary(l),l=1,npsitm)
      read(11,1003) (fpsiar_(l),l=1,npsitm)
      read(11,1003) (gpary(l),l=1,npsitm)
      read(11,1004) nx,nz,isym,iplim
      if (isym.ne.1) stop 'tdtscinp:  up-down symmtry assumed'
      read(11,1003) rleft,rright,zbot_,ztop_
      zbot_=-ztop_
      read(11,1003) psimag,psilim,raxis,zaxis
      read(11,1003) radmaj,btor,toteqd
      read(11,1003) psep,xsep,zsep
      read(11,1003) psep1,xsep1,zsep1
      nnr=nx
      nnz=2*nz-1
      if (nnr.gt.nnra.or.nnz.gt.nnza) stop 'tdtscinp: check nnr,nnz'
      read(11,1003) ((epsi(i,j),i=1,nx),j=nz,2*nz-1)

!..................................................................
!     Now expand out to full vertical mesh as required by "eq" module.
!..................................................................


      do 100 j=1,nz-1
        do 101 i=1,nnr
          epsi(i,j)=epsi(i,nnz+1-j)
 101    continue
 100  continue

!......................................................................
!     Re-represent psi on 33x65 grid using bi-linear interpolation,
!     encompassing plasma plus 2 grid points.
!......................................................................

!
!     Find region of grid containing the plasma by stepping
!     along rays issuing from the plasma magnetic axis, and
!     determining if (psi-psilim) has changed sign.

      drr=(rright-rleft)/(nnr-1)
      do 230  i=1,nnr
 230  er(i)=rleft+(i-1)*drr
      dzz=(ztop_-zbot_)/(nnz-1)
      do 231  j=1,nnz
 231  ez(j)=zbot_+(j-1)*dzz

!990131      dss=0.9*amin1(drr,dzz)
!990131      dthetap=8.*atan2(1.,1.)/201.
      dss=0.9*min(drr,dzz)
      dthetap=8.d0*atan2(one,one)/201.d0
      thetap=0.0
!990131      ismax=max1((rright-rleft)/dss,(ztop_-zbot_)/dss)
      ismax=max((rright-rleft)/dss,(ztop_-zbot_)/dss)
      jmin=nnz
      jmax=1
      iminn=nnr
      imaxx=1

      do 200  it=1,200
        thetap=thetap+dthetap
        s=0.0

        do 210  is=1,ismax
          s=s+dss
          rval=raxis+s*cos(thetap)
          zval=zaxis+s*sin(thetap)

          if (rval.gt.rright.or.rval.lt.rleft) go to 200
          if (zval.gt.ztop_.or.zval.lt.zbot_) go to 200

          kr1=(rval-rleft)/drr+1
          kz1=(zval-zbot_)/dzz+1
          kr1=min0(kr1,nnr-1)
          kz1=min0(kz1,nnz-1)
          kr2=kr1+1
          kz2=kz1+1

          jmin=min0(jmin,kz1)
          jmax=max0(jmax,kz2)
          iminn=min0(iminn,kr1)
          imaxx=max0(imaxx,kr2)

          f1=epsi(kr1,kz1)+(rval-er(kr1))* &
            (epsi(kr2,kz1)-epsi(kr1,kz1))/(er(kr2)-er(kr1))
!
          f2=epsi(kr1,kz2)+(rval-er(kr1))* &
            (epsi(kr2,kz2)-epsi(kr1,kz2))/(er(kr2)-er(kr1))
!
          val=f1+(zval-ez(kz1))*(f2-f1)/(ez(kz2)-ez(kz1))

          if ((val-psilim)*(psimag-psilim).lt.0.0) go to 200

 210    continue
 200  continue

      jmax=jmax+2
      if(jmax.ge.nnz)  go to 400
      jmin=nnz-jmax+1
      imaxx=min0(imaxx+2,nnr)
      iminn=max0(1,iminn-2)
      if (iminn.eq.1.and.imaxx.eq.nnr) go to 400

      drr=(er(imaxx)-er(iminn))/(nra-1)
      dzz=(ez(jmax)-ez(jmin))/(nza-1)
      do 232  i=1,nra
 232  ernew(i)=er(iminn)+(i-1)*drr
      do 233  j=1,nza
 233  eznew(j)=ez(jmin)+(j-1)*dzz

      kz1=1
      kz2=2
      do 250  j=1,nza
        zval=eznew(j)
 260    if (zval.le.ez(kz2)) go to 265

        kz1=kz1+1
        kz2=kz1+1
        go to 260
!
 265    continue
        kr1=1
        kr2=2
        do 270  i=1,nra
          rval=ernew(i)
 280      if (rval.le.er(kr2)) go to 285

          kr1=kr1+1
          kr2=kr1+1
          go to 280
!
 285      continue
          f1=epsi(kr1,kz1)+(rval-er(kr1))* &
            (epsi(kr2,kz1)-epsi(kr1,kz1))/(er(kr2)-er(kr1))
!
          f2=epsi(kr1,kz2)+(rval-er(kr1))* &
            (epsi(kr2,kz2)-epsi(kr1,kz2))/(er(kr2)-er(kr1))
!
          val=f1+(zval-ez(kz1))*(f2-f1)/(ez(kz2)-ez(kz1))
          epsinew(i,j)=val
 270    continue
 250  continue

!     Substitute new psi grid for old:
      nnr=nra
      nnz=nza
      do 300 j=1,nnza
        do 301 i=1,nnra
 301    epsi(i,j)=0.0
 300  continue

      do 310 j=1,nnz
        do 311 i=1,nnr
 311    epsi(i,j)=epsinew(i,j)
 310  continue

      do 320  i=1,nnra
 320  er(i)=0.0
      do 321  j=1,nnza
 321  ez(j)=0.0
      do 322  i=1,nnr
 322  er(i)=ernew(i)
      do 323  j=1,nnz
 323  ez(j)=eznew(j)
      rleft=er(1)
      rright=er(nnr)
      zbot_=ez(1)
      ztop_=ez(nnz)
      nx=nnr
      nz=(nnz+1)/2

!.......................................................................
!     Set initial values of rmag,zmag
!     (to be refined in eqrhopsi for bicubic splines of epsi).
!.......................................................................

      rmag=raxis*1.e+2
      zmag=zaxis*1.e+2

!
 400  continue
!
      iprone="disabled"
      iprote="disabled"
      iproti="disabled"
      iprozeff="disabled"
!BH070116      nnv=npsitm
      nfp=npsitm
      psimago=psimag
      psilimo=psilim
!
!BH070116      do 10 i=1,nnv
      do 10 i=1,nfp
        elecf(i)=elecf(i)/(2.*pi*radmaj*100.)
 10   continue
!
      close(unit=11)
      return
 1001 format(10a8)
 1002 format(3i5,e16.6)
 1003 format(5e16.6)
 1004 format(5i5)
      end subroutine tdtscinp


end module tdtscinp_mod
