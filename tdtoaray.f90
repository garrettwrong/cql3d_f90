module tdtoaray_mod

!
!

contains

      subroutine tdtoaray
      use param_mod
      use comm_mod
      implicit integer (i-n), real*8 (a-h,o-z)

!     Storing values on limited lrz mesh and time-dependent quantities
!     recall:
!     lr_     : radial index
!     l_      : spatial variable index
!     lmdpln_ : midplane at l_=lmdpln_
!


      if (kelecg .gt. 0) vfluxz(l_)=vflux(jx-1,1,l_)
      currmtz(l_)=currmt(l_)/3.e+9
      currmtpz(l_)=currmtp(l_)/3.e+9

!..................................................................
!     Copies flux surface dependent quantities into radial arrays
!..................................................................

      if (l_ .ne. lmdpln_) return

      do 5 k=1,ngen
        pegyz(k,lr_)=entr(k,12,l_)
        pplossz(k,lr_)=entr(k,6,l_)
 5    continue
!BH081202      psyncz(lr_)=entr(kelecg,11,l_)
      if (kelecg.ne.0) then
         psyncz(lr_)=entr(kelecg,11,l_)
      else
         psyncz(lr_)=zero
      endif

!..................................................................
!     rfpwrz(k,lr_) is the local rf power deposited in species k
!     at flux surface lr_ in Watts/cc (averaged over the orbit)
!     currz(k,lr_) is the current density of species k at lr_ in A/cm**2
!       (averaged over incremental toroidal-area at given flux surface).
!     currtz(lr_) is the sum of currz(k,lr_) over ionic general species.
!     currmtz(lmdpln_) is the sum of the midplane currents over ionic species
!     curr(m)tpz(lr_) is curr(m)tz plus the electron contribution to
!     the current.
!     1 Amp = 3.e9 statAmps
!..................................................................


      do 21 k=1,ngen
        currz(k,lr_)=curr(k,lr_)/3.e+9
        rfpwrz(k,lr_)=entr(k,3,l_)
        if (nso.gt.0) then
          numsrce=nso
        else
          numsrce=1
        endif
        do 23 m=1,numsrce
          asorz(k,m,lr_)=asor(k,m,lr_)
 23     continue
 21   continue

      currtz(lr_)=currt(lr_)/3.e+9
      currtpz(lr_)=currtp(lr_)/3.e+9

!..................................................................
!     Compute the plasma pressure (cgs)
!     if (kpress(k).ne. "enabled") this species is not used
!     in calculation of total pressure.
!..................................................................

      if (eqmod.eq."enabled") then
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
end module tdtoaray_mod