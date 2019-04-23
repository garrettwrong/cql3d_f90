module psif_mod

  !---BEGIN USE

  use zcunix_mod, only : terp1

  !---END USE

!
!
!

contains

      real*8 function psif(xz)
      use param_mod
      use comm_mod
      implicit integer (i-n), real*8 (a-h,o-z)
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
          call terp1(lorbit(lr_),es(1,lr_),bpsi(1,lr_),d2bpsi(1,lr_), &
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
          call terp1(lorbit(lr_),es(1,lr_),bpsi(1,lr_),d2bpsi(1,lr_), &
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
      end
end module psif_mod
