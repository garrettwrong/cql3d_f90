module tdtoarad_mod

!
!

contains

      subroutine tdtoarad
      use param_mod
      use comm_mod
      implicit integer (i-n), real*8 (a-h,o-z)

!     Storing values on whole mesh lrzmax and equilibrium quantities
!     recall:
!     lr_     : radial index
!     l_      : spatial variable index
!     lmdpln_ : midplane at l_=lmdpln_
!


!..................................................................
!     BH091022:  Looks like a somewhat orphaned idea.  The "-z"
!                variable here a not used much elsewhere, and
!                original (non-z) variables could be used.
!                Seems to create confusion of variable names.
!                Clean it up sometime.....
!                No references to l_ or lmdpln_?
!..................................................................

!..................................................................
!     Copies flux surface dependent quantities into radial arrays
!..................................................................

!..................................................................
!     rfpwrz(k,lr_) is the local rf power deposited in species k
!     at flux surface lr_ in Watts/cc (averaged over the orbit)
!     currz(k,lr_) is the current of species k at lr_ in Amps (averaged
!     currtz(lr_) is the sum of currz(k,lr_) over ionic general species.
!     currmtz(lmdpln_) is the sum of the midplane currents over ionic species
!     curr(m)tpz(lr_) is curr(m)tz plus the electron contribution to
!     the current.
!..................................................................


      do 21 k=1,ngen
        if (nso.gt.0) then
          numsrce=nso
        else
          numsrce=1
        endif
        do 23 m=1,numsrce
          asorz(k,m,lr_)=asor(k,m,lr_)
 23     continue
 21   continue

      if (eqmod.eq."enabled") then
        area(lr_)=areacon(lr_)
        vol(lr_)=volcon(lr_)
        equilpsi(lr_)=epsicon(lr_)
        bmdplne(lr_)=bmidplne(lr_)
        onovrpz(lr_,1)=onovrp(1,lr_)
        onovrpz(lr_,2)=onovrp(2,lr_)
        aspin(lr_)=(rpcon(lr_)-rmcon(lr_))/(rpcon(lr_)+rmcon(lr_))
        rmconz(lr_)=rmcon(lr_)
        rpconz(lr_)=rpcon(lr_)
        fpsiz(lr_)=fpsi(lr_)
        bpolsqaz(lr_)=bpolsqa(lr_)

!..................................................................
!     Compute the plasma pressure (cgs)
!     if (kpress(k).ne. "enabled") this species is not used
!     in calculation of total pressure.
!..................................................................

        phot=0.
        do 10 k=1,ngen
          if (kpress(k).eq."disabled") go to 10
          phot=phot+reden(k,lr_)*energy(k,lr_)*2./3.*ergtkev
 10     continue
        prest(lr_)=phot
        do 20 k=ngen+1,ntotal
          if (kpress(k).eq."disabled") go to 20
          prest(lr_)=prest(lr_)+reden(k,lr_)*energy(k,lr_)*2./3.* &
            ergtkev
 20     continue
 30     continue
      endif
      return
      end
end module tdtoarad_mod