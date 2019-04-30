module tdtscout_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use tdinterp_mod, only : tdinterp

  !---END USE

!
!

contains

      subroutine tdtscout
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     This routine creates a file 'tscout'; data
!     is read in from the file created by the TSC code. The file
!     created here has the power and the current as a function of
!     rho_, the radial coordinate given by TSC.
!..................................................................

!MPIINSERT_INCLUDE

      dimension powtsc(nrada),currtsc(nrada)

!MPIINSERT_IF_RANK_NE_0_RETURN

      if (eqsource.ne."tsc") return

!..................................................................
!     Create output file.
!..................................................................

      open(unit=12,file='tscout',delim='apostrophe',status='unknown')

!..................................................................
!     interpolate to TSC mesh
!..................................................................

      do 10 l=1,lrzmax
        tr(l)=sorpw_rf(kelecg,l)
 10   continue
      call tdinterp("free","free",rya(1),tr(1),lrzmax,rho_,powtsc, &
        npsitm)
      call tdinterp("free","free",rya(1),currtpz(1),lrzmax,rho_,currtsc &
        ,npsitm)

!..................................................................
!     Now write to disk.
!     powtsc is the power in Watts/cc
!     currtsc is the current in Amps/cm**2
!..................................................................

      write(12,100) (powtsc(l),l=1,npsitm)
      write(12,100) (currtsc(l),l=1,npsitm)
 100  format(5e16.6)
      close(unit=12)
      return
      end
end module tdtscout_mod
