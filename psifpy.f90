module psifpy_mod

  !---BEGIN USE

  !---END USE

!

contains

      real*8 function psifpy(yval)
      use param_mod
      use comm_mod
      implicit integer (i-n), real*8 (a-h,o-z)
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
          xt=psiinv(yval)
          psifpy=psifp(xt)
        endif
      endif
      return
      end
end module psifpy_mod
