module psif_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double
  use iso_c_binding, only : c_double

  use diagwrng_mod, only : diagwrng
  use zcunix_mod, only : terp1

  !---END USE

  !
  !
  !

contains

  real(c_double) function psif(xz)
    use param_mod
    use comm_mod
    implicit integer (i-n), real(c_double) (a-h,o-z)
    save
    !-------------------------------------------------------------------
    !
    !mnt  function psif(arclength)
    !mnt  calculates the ratio of mod B at arclength from the midplane
    !mnt  to mod B at the midplane; 1 < psi < psimx(lr_).
    !mnt  the value of psimodel controls which model is expressed:
    !mnt  options are:
    !mnt  psimodel="axitorus", "smith", "spline"
    !
    !     OK with eqsym.eq."none".  xz is measured from es(1) and in
    !       this case covers a whole poloidal turn.
    !
    !-------------------------------------------------------------------

    dimension tabl(3),itabl(3)
    data itabl(1) /1/, itabl(2) /0/, itabl(3) /0/
    if (machine.eq."toroidal") then
       if (psimodel.eq."axitorus") then
          psif=(1.+eps(lr_))/(1.+eps(lr_)*cos(pi*xz*zmaxi(lr_)))
       else if (psimodel.eq."spline") then
          call terp1(lorbit(lr_),es(1:lorbit(lr_),lr_),bpsi(1:lorbit(lr_),lr_),d2bpsi(1:lorbit(lr_),lr_), &
               xz,1,tabl,itabl)
          psif=tabl(1)
       endif
    else if (machine.eq."mirror") then
       if (psimodel.eq."axitorus") then
          psif=(1.+eps(lr_))/(1.+eps(lr_)*cos(pi*xz*zmaxi(lr_)))
       else if (psimodel.eq."smith") then
          gsexm=exp(-((xz-gszb)/gslb)**2)
          gsexp=exp(-((xz+gszb)/gslb)**2)
          psif=(1.+(xz/gsla)**2+gsb*(gsexm+gsexp))/gsnm
       else if (psimodel.eq."spline") then
          call terp1(lorbit(lr_),es(1:lorbit(lr_),lr_),bpsi(1:lorbit(lr_),lr_),d2bpsi(1:lorbit(lr_),lr_), &
               xz,1,tabl,itabl)
          psif=tabl(1)
       endif
    endif
    if (trapmod.eq."enabled") then
       !        Reduce B/B_min by the amount trapredc*(B/B_min-1.)
       !        This is for purposes of a heuristic check of trapping effects.
       psif=psif-trapredc*(psif-1.)
    endif
    return
  end function psif

  real(c_double) function psifp(xz)
    use param_mod
    use comm_mod
    implicit integer (i-n), real(c_double) (a-h,o-z)
    save
    !-------------------------------------------------------------------
    !
    !mnt  function psifp(arclength)
    !mnt  calculates the derivative of psi=B/B0 w.r.t. arclength, as a
    !mnt  function of arclength: psifp(s)=dpsi/ds (s)
    !mnt  model fields are as for psif.
    !
    !     OK with eqsym.eq."none".  xz is measured from es(1) and in
    !       this case covers a whole poloidal turn.  The function
    !       will be negative for xz>es_bmax(lr_).
    !
    !-------------------------------------------------------------------

    itab(1)=0
    itab(2)=1
    itab(3)=0
    if (machine.eq."toroidal") then
       if (psimodel.eq."axitorus") then
          psifp=psifpy(psif(xz))
       else if (psimodel.eq."spline") then
          call terp1(lorbit(lr_),es(1:lorbit(lr_),lr_),bpsi(1:lorbit(lr_),lr_),d2bpsi(1:lorbit(lr_),lr_), &
               xz,1,tab,itab)
          psifp=tab(2)
          if (trapmod.eq."enabled") then
             !     Reduce B/B_min by the amount trapredc*(B/B_min-1.)
             !     This is for purposes of a heuristic check of trapping effects.
             psifp=psifp-trapredc*psifp
          endif
       endif
    else if (machine.eq."mirror") then
       if (psimodel.eq."axitorus") then
          psifp=psifpy(psif(xz))
       else if (psimodel.eq."smith") then
          gszxm=xz-gszb
          gszxp=xz+gszb
          gsexm=exp(-gszxm**2/gslb2)
          gsexp=exp(-gszxp**2/gslb2)
          psifp=2.*(xz/gsla2-gsb*(gszxm*gsexm+gszxp*gsexp)/gslb2)/gsnm
       else if (psimodel.eq."spline") then
          call terp1(lorbit(lr_),es(1:lorbit(lr_),lr_),bpsi(1:lorbit(lr_),lr_),d2bpsi(1:lorbit(lr_),lr_), &
               xz,1,tab,itab)
          psifp=tab(2)
          if (trapmod.eq."enabled") then
             !     Reduce B/B_min by the amount trapredc*(B/B_min-1.)
             !     This is for purposes of a heuristic check of trapping effects.
             psifp=psifp-trapredc*psifp
          endif
       endif
    endif
    return
  end function psifp

  real(c_double) function psiinv(yval,iupdown)
    use param_mod
    use comm_mod
    implicit integer (i-n), real(c_double) (a-h,o-z)
    save
    !---------------------------------------------------------------------
    !
    !mnt  function psiinv(psi)
    !mnt  calculates the arclength s associated with given value of
    !mnt  yval=psi=B(s)/B0,  ie., performs the functional inversion of
    !mnt  psi(s,lr_):   thus psi(psiinv(psi))=psi or psiinv(psi(s))=s.
    !
    !     yval is a psi(=B(s)/B0) value.
    !     Note: Assuming u=u0 (energy conservation),
    !           we have B(s)/B0= [sin(theta)/sin(theta0)]^2
    !           For a trapped particle: B(s)/B0= [1/sin(theta0)]^2
    !           at the bouncing point.
    !
    !     iupdown=1 (upper equilibrium) or -1 (lower), for eqsym="none".
    !     Arclength is measured counter-clockwise from s=0 at the
    !     minimum |B|-pt, in either case.
    !
    !---------------------------------------------------------------------

    data toll/1.e-8/,nhalve/10/
    niter=0

    isign=+1
    if (eqsym.eq."none" .and. iupdown.eq.-1) isign=-1

    if (yval.lt.1. .or. yval.gt.psimx(lr_)) then
       write(*,*)'psiinv: lr_, yval, psimx(lr_)=', lr_, yval,psimx(lr_)
       call diagwrng(7)
    else if (yval.eq.1.) then
       psiinv=0.
       return
    else if (abs(psimx(lr_)-yval).lt.1.e-13) then
       psiinv=z_bmax(lr_)
       return
    else
       if (isign.eq.1) then
          sold=((yval-1.)/(psimx(lr_)-1.))*z_bmax(lr_)
       else
          sold=((yval-1.)/(psimx(lr_)-1.))*(zmax(lr_)-z_bmax(lr_))
       endif
    endif

    if (isign.eq.1) then
       sl=0.
       sr=zmax(lr_)
    else
       sl=zmax(lr_)
       sr=z_bmax(lr_)
    endif
15  continue
    do 5 i=1,nhalve
       yo=psif(sold)
       if (yo.lt.yval) then
          sl=sold
       else
          sr=sold
       endif
       sold=(sr+sl)*0.5
5      continue
       relerr=abs((sl-sr)/sr)
       if(relerr.lt.toll) then
          psiinv=sold
          return
       endif
       !     Newton-Raphson interation:
10     snew=sold-(psif(sold)-yval)/psifp(sold)  !psifp has sign=isign
       if (((snew.gt.sr).or.(snew.lt.sl)).and.isign.eq.1) then
          sold=(sl+sr)*0.5
          go to 15
       endif
       if (((snew.gt.sl).or.(snew.lt.sr)).and.isign.eq.-1) then
          sold=(sl+sr)*0.5
          go to 15
       endif
       relerr=abs((snew-sold)/sold)
       if (relerr.lt.toll) go to 30
       niter=niter+1
       if (niter.gt.100) call diagwrng(2)
20     continue
       sold=snew
       go to 10
30     psiinv=snew
       return

  end function psiinv

  real(c_double) function psifpy(yval)
      use param_mod
      use comm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

      !------------------------------------------------------------
      !
      !mnt  calculates the derivative of bbpsi with respect to arc
      !mnt  length along a field line-as a function of bbpsi:
      !mnt  psifpy(bbpsi)=dbbpsi/ds (bbpsi)
      !
      !BH091106:  Only called from psifp, and only for psimodel.eq."axitorus".
      !
      !------------------------------------------------------------


      iupdown=1

      if (machine.eq."toroidal") then
         if (psimodel.eq."axitorus") then
            psifpy=sqrt((yval-1.)*(1.-(yval/psimx(lr_))))*yval* &
                 pibzmax(lr_)
         else if (psimodel.eq."spline") then !But only called with axitorus
            xt=psiinv(yval,iupdown)
            psifpy=psifp(xt)
         endif
      else if (machine.eq."mirror") then
         if (psimodel.eq."axitorus") then
            psifpy=sqrt((yval-1.)*(1.-(yval/psimx(lr_))))*yval* &
                 pibzmax(lr_)
         else if (psimodel.eq."smith") then
            xz=psiinv(yval,iupdown)
            psifpy=psifp(xz)
         else if (psimodel.eq."spline") then !But only called with axitorus
            xt=psiinv(yval, iupdown) !XXXXX BUG, missing iupdown
            psifpy=psifp(xt)
         endif
      endif
      return
  end function psifpy

  real(c_double) function psifppy(yval,iupdown)
    use param_mod
    use comm_mod
    implicit integer (i-n), real(c_double) (a-h,o-z)
    save

    !------------------------------------------------------------
    !
    !mnt  function psifppy(psi)
    !mnt  calculates the second derivative of psi w/r/t arc length
    !mnt  at the midplane or the throat as a function of psi
    !mnt  psifppy(psi)=d(dpsi/ds)/ds (psi)

    !     NOTE:  Not presently used.
    !            If to be re-commisioned for some purpose,
    !            update/check for use with trapmod="enabled".
    !            Also, update for iupdown (for eqsym.eq."none").
    !
    !
    !
    !------------------------------------------------------------


    itab(1)=0
    itab(2)=0
    itab(3)=1
    if (machine.eq."toroidal") then
       if (psimodel.eq."axitorus") then
          if (yval.eq.1.) then
             psifppy=.5*(1.-1./psimx(lr_))*pibzmax(lr_)**2
          else if (abs(psimx(lr_)-yval).lt.1.e-13) then
             psifppy=-.5*psimx(lr_)*(psimx(lr_)-1.)*pibzmax(lr_)**2
          else
             call diagwrng(6)
          endif
       else if (psimodel.eq."spline") then
          if (yval.eq.1.) then
             xz=0.
          else if (abs(psimx(lr_)-yval).lt.1.e-13) then
             xz=zmax(lr_)
          else
             call diagwrng(6)
          endif
          call terp1(lorbit(lr_),es(1:lorbit(lr_),lr_),bpsi(1:lorbit(lr_),lr_),d2bpsi(1:lorbit(lr_),lr_), &
               xz,1,tab,itab)
          psifppy=tab(3)
       endif
    else if (machine.eq."mirror") then
       if (psimodel.eq."axitorus") then
          if (yval.eq.1.) then
             psifppy=.5*(1.-1./psimx(lr_))*pibzmax(lr_)**2
          else if (abs(psimx(lr_)-yval).lt.1.e-13) then
             psifppy=-.5*psimx(lr_)*(psimx(lr_)-1.)*pibzmax(lr_)**2
          else
             call diagwrng(6)
          endif
       else if (psimodel.eq."smith") then
          if (yval.eq.1.) then
             gsa3=2.*(1.+gszb2/gslb2)*gseb
             psifppy=2.*(1./gsla2+2.*gsb*gsa3/gslb2)/gsnm
          else if (abs(psimx(lr_)-yval).lt.1.e-13) then
             gsa3=(1.+gszm2/gslb2)*gsem+(1.+gszp2/gslb2)*gsep
             psifppy=2.*(1./gsla2+2.*gsb*gsa3/gslb2)/gsnm
          else
             call diagwrng(6)
          endif
       else if (psimodel.eq."spline") then
          if (yval.eq.1.) then
             xz=0.
          else if (abs(psimx(lr_)-yval).lt.1.e-13) then
             xz=zmax(lr_)
          else
             call diagwrng(6)
          endif
          call terp1(lorbit(lr_),es(1:lorbit(lr_),lr_),bpsi(1:lorbit(lr_),lr_),d2bpsi(1:lorbit(lr_),lr_), &
               xz,1,tab,itab)
          psifppy=tab(3)
       endif
    endif
    return
  end function psifppy


end module psif_mod
