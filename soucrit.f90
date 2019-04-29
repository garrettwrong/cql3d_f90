module soucrit_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double
  use iso_c_binding, only : c_double

  use cfpgamma_mod, only : cfpgamma

  !---END USE

!
!

contains

      subroutine soucrit
      use param_mod
      use comm_mod
      use r8subs_mod, only : luf
      implicit integer (i-n), real(c_double) (a-h,o-z)
!
!   scchiu, 9609..
!  calculate the critical momentum-per-mass.
!  Since jxcrit is used for calculation of relativistic runaways,
!   it is chosen to be no larger than for particles at 3.*clight.
!   This is to reduce a problem of definition which occurs
!   transiently for abs(E) < abs(E_crit_for_runaway).
!
      call cfpgamma
      do 10 k=1,ngen
        fack=abs(elecfld(lr_))/reden(k,lr_)*1.e16/0.0918 &
                    *18./gama(k,k)
        eoe0(k,lr_)=fack

            if((fack-1.).eq.zero)then
             write(*,*)'soucrit: (fack-1.)=0=',(fack-1.)
             !pause
            endif

        fack1=1./(fack-1.)
        if (fack1.le.0.d0) then
          ucrit(k,lr_)=1.e20
        else
          ucrit(k,lr_)=clight*sqrt(fack1)/vnorm
        endif
!  Take runaway electrons to have momentum per mass beyond the
!  minimum of 3.*clight or ucrit:
!990131        xcrit=amin1(3.*clight/vnorm,ucrit(k,lr_))
        xcrit=min(3.*clight/vnorm,ucrit(k,lr_))
        jxcrit(k,lr_)=luf(xcrit,x,jx)
10    continue
      return
      end
end module soucrit_mod
