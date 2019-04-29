module cfpgamma_mod

  !---BEGIN USE

  use comm_mod
  use iso_c_binding, only : c_double
  use param_mod

  !---END USE

contains

  subroutine cfpgamma
    !implicit integer (i-n), real*8 (a-h,o-z)
    implicit none
    real(c_double) :: deby
    real(c_double) :: dense
    real(c_double) :: energ
    real(c_double) :: flam
    real(c_double) :: flamcql
    real(c_double) :: flamrp1
    real(c_double) :: gam3
    real(c_double) :: omegape
    real(c_double) :: sf
    real(c_double) :: si
    real(c_double) :: sk
    real(c_double) :: vikf
    integer :: i
    integer :: j
    integer :: k
    integer :: km
    

    !.......................................................................
    !     calculate the Coulomb logarithm and an energy dependent factor.
    !     BH180807: Needs reviewing.  Most often gamaset=fixed value (16.
    !     or 17.) is used.  Should re-investigate ln(Lambda) expressions.
    !     Compare/code NRL formulary expressions.
    !.......................................................................



    !..................................................................
    !     If log is set by fiat,  set it and move to gamefac
    !..................................................................

    if (gamaset.gt.zero) then
       do 12 i=1,ntotal
          do 13 k=1,ntotal
             gama(i,k)=gamaset
13        end do
12     end do
       goto 99
    endif


    energ=ergtkev*energy(kelec,lr_)
    dense=reden(kelec,lr_)
    if (cqlpmod .eq. "enabled") then
       energ=ergtkev*enrgypa(kelec,ls_)
       dense=denpar(kelec,ls_)
    endif

    !..................................................................
    !     If two electron species exist - do an average for deby calc.
    !..................................................................
    
    if (colmodl.eq.1 .or. (colmodl.eq.3 .and. kelecm.ne.0)) then
       dense=reden(kelecm,lr_)
       energ=reden(kelecm,lr_)*energy(kelecm,lr_)*ergtkev/dense
       if (cqlpmod .eq. "enabled") then
          dense=denpar(kelecm,ls_)
          energ=enrgypa(kelecm,ls_)*ergtkev
       endif

    elseif (kelecg.gt.1 .and. kelecm.gt. 1) then
       dense=reden(kelec,lr_)+reden(kelecm,lr_)
       km=kelecm
       energ=(reden(kelecg,lr_)*energy(kelecg,lr_)+reden(km,lr_)* &
            energy(km,lr_))*ergtkev/dense
       if (cqlpmod .eq. "enabled") then
          dense=denpar(kelec,ls_)+denpar(kelecm,ls_)
          energ=(denpar(kelecg,ls_)*enrgypa(kelecg,ls_)+denpar(km,ls_)* &
               enrgypa(km,ls_))*ergtkev/dense
       endif
    endif
    deby=sqrt(energ/(dense*6.*pi*charge**2))
    do 10 i=1,ntotal
       si=energy(i,lr_)/fmass(i)
       if (cqlpmod .eq. "enabled") si=enrgypa(i,ls_)/fmass(i)
       do 11 k=1,ntotal
          sk=energy(k,lr_)/fmass(k)
          if (cqlpmod .eq. "enabled") sk=enrgypa(k,ls_)/fmass(k)
          !990131          sf=amax1(si,sk)
          sf=max(si,sk)
          vikf=sqrt(sf*ergtkev*2.)          
          gam3=gamt(i,k)*deby*vikf
          !990131          gama(i,k)=(alog(gam3)-.5)
          gama(i,k)=(log(gam3)-.5)
11     end do
10  end do


99  continue

    !..................................................................
    !     If kelecg.eq.1 and gamafac.eq."enabled", set energy
    !       dependant factor gamefac(j) for Coulomb logrithm.
    !     The energy dependent Coulomb log will be
    !       alog(flamcql + min(sqrt(gamma-1),1)*(flamrp-flamcql)),
    !       where flamcql is the argument of cql e-e Coulomb log,
    !       and flamrp is the argument of the Rosenbluth-Putvinski
    !       (Nucl. Fus. 1997) high energy electron Coulomb log.
    !       The sqrt(gamma-1)-factor is chosen in accord with
    !       the energy factor in the the CQL Coulomb log
    !       (cf., Killeen et al. book).  The CQL Coulomb log
    !       reduces to the NRL, Te.gt.10eV expression.
    !     Some further incites into this factor can be obtained
    !       from Physics Vade Mecum, H.L. Anderson, Ed.
    !       2nd edition, AIP, NY (1989); Sect. 16.07.F, "energy
    !       loss due to scattering from atomic electrons ... given
    !       by Moller scattering, with I=(electronic charge)**2/(
    !       Debye length) [Rosenbluth, personal comm. 1996].
    !       (BH, 980505)
    !..................................................................
    
    do 14 j=1,jx
       gamefac(j)=1.0
14  end do

    if (kelecg.eq.1 .and. ngen.eq.1 .and. gamafac.eq."enabled") then
       
       flamcql=exp(gama(1,1))
       omegape=5.64e4*sqrt(reden(1,lr_))
       flamrp1=(2**.25)*fmass(1)*clite2/(1.0546e-27*omegape)
       
       do  j=1,jx
          !990131            flam=flamcql+min(sqrt(gamma(j)-1.),1.)*
          flam=flamcql+min(sqrt(gamma(j)-1.),one)* &
               (flamrp1*gamma(j)**1.5 - flamcql)
          !990131            gamefac(j)=alog(flam)/gama(1,1)
          gamefac(j)=log(flam)/gama(1,1)
       enddo
    endif
    
    
    return
  end subroutine cfpgamma

end module cfpgamma_mod
