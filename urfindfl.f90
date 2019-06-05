module urfindfl_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine urfindfl
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     This routine sets defaults for the "urf" module
!..................................................................



!..................................................................
!     Set defaults
!..................................................................

      one=1.d0

      nrayn=1 !-YuP 101122: added here, instead of nraya in param.h;
              ! will be re-defined in urfsetup.
              ! Will be used to allocate arrays, mostly in urfalloc
      nrayelts=1 !-YuP: added here, instead of nrayelta in param.h;
              ! will be re-defined in urfsetup.
              ! will be used to allocate arrays, mostly in urfalloc.

!..................................................................
!     ieqbrurf designates the source of equilibrium data to be used by
!     xbram.  Appropriate values are: ieqbrurf=0 to use Brambilla
!     ieqbrurf=0, to use Brambilla analytic "equilibria",
!     =3, to use standard eqdsk input.
!     =4, to use extended eqdsk, including density, and
!     temperature profiles,...[output by GA ONETWO transport code].
!     If eqsource="tsc", ieqbrurf is reset to = 4.
!..................................................................

      ieqbrurf=1

!..................................................................
!     Step interval at which to recalc diffusion coefficients,
!     as a fraction of ncoef.
!..................................................................

      urfncoef=one

!..................................................................
!     number of elements in Bessel table at setup time
!..................................................................

      nbssltbl=10000

!..................................................................
!     damping on at step nondamp
!..................................................................

      nondamp=0

!..................................................................
!     maximum number of ray elements per ray for first call to R.T code
!     [Not presently implemented, BH060314]
!..................................................................

      do i=1,nmodsa
         nrfstep1(i)=100
      enddo
!..................................................................
!     scaling of power and nparallel- width in rays
!..................................................................

      do 20  i=1,nmodsa
        pwrscale(i)=one
        wdscale(i)=one
 20   continue

!..................................................................
!     nharms.gt.0, then calculate damping for harmonics
!     nharm1 to nharm1+(nharms-1), rather than according to nharm
!     in the rayop file.
!     This option is only viable for one rf mode, i.e., only one
!     of lh, ech, or fw may be enabled.
!     060214: Extended to multi-rf-wave types.
!..................................................................

      do i=1,nmodsa
         nharms(i)=0
         nharm1(i)=0
      enddo

!..................................................................
!     Additional time-dependent scaling of power
!..................................................................

      nurftime=0
      do  i=1,nbctimea
        pwrscale1(i)=one
        urftime(i)=0.0
      enddo

!..................................................................
!     number of steps in power ramp-up
!..................................................................

      nrfpwr=3

!..................................................................
!     number of iterations at full power, for first ray elements
!..................................................................

      nrfitr1=1

!..................................................................
!     number of iterations after additional ray data
!..................................................................

      nrfitr2=1

!..................................................................
!     number of steps adding new ray data
!..................................................................

      nrfitr3=2

!..................................................................
!     number of ray elements added at each step in addition of new data.
!     [Not presently implemented, BH060314]
!..................................................................

      nrfstep2=50

!..................................................................
!     scaleurf="enabled" rescale contribution to B so that a particular
!     ray does not "overdamp" on a given flux surface volume.
!..................................................................

      scaleurf="enabled"

!..................................................................
!     iurfcoll (iurfl) indicate use of collisional (additional linear)
!     absorption coeff. Passed in ray data files.
!..................................................................

      do i=1,nmodsa
         iurfcoll(i)="disabled"
         iurfl(i)="disabled"
      enddo

!..................................................................
!     lh,ech and fw determine which wave modes are utilized
!     Alternatively(060314), this data can be input through rftype()
!..................................................................

      lh="disabled"
      ech="disabled"
      fw="disabled"
      do i=1,nmodsa
         rftype(i)="notset"
         rffile(i)="notset"
         nrfspecies(i)=1
      enddo

!..................................................................
!     Setting nmods=nmodsa, for the time being. Generalize later.
!..................................................................

      nmods=nmodsa ! YuP-101220: should be mrfn, but not known yet

!..................................................................
!     Variables determine which, if any,  ray tracing code is called.
!..................................................................

      call_lh="disabled"
      call_ech="disabled"
      call_fw="disabled"

      rfread="text"

!..................................................................
!     urfdmp="secondd" means utilize the "second derivative" damping
!     diagnostic in computing the damping due to each ray as
!     it moves through the plasma. If urfdmp .ne. "secondd" the
!     diagnostic uses an integration by parts technique to compute
!     the damping. We highly recommend "secondd" because convergence
!     is the agreement between the sum of the damping of all rays
!     and the absorbed power as computed from dF/dt. This latter
!     diagnostic utilizes the "second derivative" approach so
!     consistency demands "secondd" for the rays.
!..................................................................

      urfdmp="secondd"
      return
      end subroutine urfindfl
      
end module urfindfl_mod
