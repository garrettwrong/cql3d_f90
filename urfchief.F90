module urfchief_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use tdnflxs_mod, only : tdnflxs
  use urfavg_mod, only : urfavg
  use urfb0_mod, only : urfb0
  use urfbes_mod, only : urfbes
  use urfbplt_mod, only : urfbplt
  use urfdamp0_mod, only : urfdamp0
  use urffflx_mod, only : urffflx
  use urfpack_mod, only : urfpack
  use urfrays_mod, only : urfrays
  use urfread_mod, only : urfread
  use urfwrite_mod, only : urfwrite
  !---END USE


!
!

contains

  subroutine urfchief
    use cqlconf_mod, only : setup0
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!..................................................................
!     This routine controls the urf module
!     The following control variables are used in determining the
!     action of this subroutine:
!     nrfstep1, nrfstep2, nrfpwr, nrfitr1, nrfitr2, nrfitr3, urfncoef
!
!     nurf=a counter for the number of calls to urfchief which have
!     resulted in calculation or recalculation of the diffusion
!     coefficients. (Initialized to 0).
!     The nurf variable is incremented (at end of urfchief) each time
!     n/integer(urfncoef*ncoef)*integer(urfncoef*ncoef).eq.n.
!
!     The sequence of actions (as a function of
!     increasing nurf.ge.0 at each call) is given by the following steps
!     (after each diffusion coeff calc, control returns to the calling
!     subroutine):
!     1. For nurf=0, calc or read ray data for nrfstep1 steps along ray.
!     2. Calc. damping of ray data, and then resulting quasilinear
!     diffusion coeffs, using a fraction of the input power
!     = (1/2)**nrfpwr.
!     return.
!     3. Repeat step 2 for next nrfpwr calls, but with fractional input
!     power (1/2)**(nrfpwr-1), (1/2)**(nrfpwr-2),.... (1/2)**0.
!     (This step is a no-op if nrfpwr=0).
!     4. Iterate step 2 with full input power for next nrfitr1 calls.
!     (This step is a no-op if nrfitr1=0).
!     5. Extend extendable rays by nrfstep2 steps.
!     6. Re-calc damping from ray data and then quasilinear diffusion
!     coeffs.  Iterate this step nrfitr2 additional calls.
!     7. Steps 5 and 6 are carried out nrfitr3 times.
!
!     Thus choose
!     nstop=
!     (nrfpwr+1+nrfitr1+nrfitr3*(nrfitr2+1))*integer[urfncoef*ncoef]
!     if the above sequence is to be completed.
!
!     This urf module treats multiple wave types, multiple cyclotron
!     harmonics for each wave type, and there may be multiple
!     rf diffused general (FP'd) species.   We refer to the sum
!     of the wave harmonics over the wave types, as the wave modes,
!     that is, one wave mode for each harmonic of each wave type.
!     The following variables are used in the urf*.f routines to
!     direct the logic:
!        mrf= number of rf "types"
!        mrfn= number of rf "modes" (sum over mrf of the nharms())
!        irfn(1:mrf)=  rf mode index (in 1:mrfn) of the lowest
!                      harmonic for each wave type
!        krfn(1:mrfn)=  wave type index (in 1:mrf) for each rf mode
!        nharm1(1:mrf)= lowest cyclotron harmonic, for each rf type
!        nharms(1:mrf)= number of cyclotron harmonics, for each rf type
!        nharm(1:mrfn)= harmonic number for each "mode".
!        nrfspecies(1:mrf) = general species index which each rf
!                            type is applied to (a nml variable)
!
!..................................................................

#ifdef __MPI
      include 'mpilib.h'
#endif

!...................................................................
!     Control logic
!...................................................................

      imprf=0
      iurfncf=urfncoef*ncoef
      if(n/iurfncf*iurfncf.ne.n)  return

!.......................................................................
!     urfrstrt.eq."enabled"  => "do not update delpwr" option
!     For convergence studies with previously calculated ray data files...
!.......................................................................

      iopt=1
      if (urfrstrt.eq."enabled") iopt=2

!pu      call cpu_time(t_urf1) !-.-.-.-.-.-.-.-.-.-.-.-.-.

!...................................................................
!     Get a time average of the distribution function for purposes
!     of computing wave deposition and damping.
!...................................................................
      call urfavg


      if(nurf.eq.0)  then

!...................................................................
!     Initialize counters for control logic
!...................................................................

        irfpwr=nrfpwr
        irfitr1=nrfitr1
        irfitr2=0
        irfitr3=nrfitr3

!...................................................................
!     Calculate the ray data, depending on call_lh,call_ech,call_fw
!...................................................................

        initrf=0
!Following needs reconfiguring for arbitrary nmodsa [BH, 060314]
!$$$        nraypts1=nrfstep1(1)
!$$$        nraypts2=nrfstep1(2)
!$$$        nraypts3=nrfstep1(3)
!$$$        call urfrays(initrf,nraypts1,nraypts2,nraypts3)

!...................................................................
!     Read in ray data
!...................................................................
        write(*,*)'call urfread: n, nurf=',n,nurf
        call urfread

!...................................................................
!     Obtain segregation of ray data by flux surface
!...................................................................

        call urffflx

!...................................................................
!     Save initial power in the rays
!...................................................................

        do 10 krf=1,mrf
          kk=irfn(krf)       !irfn setup in urfread
          do 11  iray=1,nray(kk)
            delpwr0(iray,kk)=delpwr(1,iray,kk)
 11       continue
 10     continue

!...................................................................
!     Initialize Bessel function table
!...................................................................

        call urfbes

!...................................................................
!     Initialize packed data and do damping
!...................................................................

        call urfpack

        call urfdamp0(irfpwr,iopt)

      elseif (iopt .eq. 2) then

!.......................................................................
!     compute damping and power but do not update delpwr according
!     to what has been absorbed before => compatible with setup0%lrz<setup0%lrzmax
!.......................................................................

        call urfdamp0(irfpwr,iopt)

      elseif (irfpwr.gt.0)  then
        irfpwr=irfpwr-1

!...................................................................
!     do damping
!...................................................................

        call urfdamp0(irfpwr,iopt)

      elseif (irfitr1.gt.0)  then
        irfitr1=irfitr1-1

!...................................................................
!     do damping
!...................................................................

!pu      call cpu_time(t_urf2) !-.-.-.-.-.-.-.-.-.-.-.-.-.
!pu      write(*,*)'urfchief t2=', t_urf2-t_urf1

        call urfdamp0(0,iopt)

!pu      call cpu_time(t_urf3) !-.-.-.-.-.-.-.-.-.-.-.-.-.
!pu      write(*,*)'urfchief t3=', t_urf3-t_urf1

      elseif (irfitr3.gt.0 .and. irfitr2.eq.0)  then
        irfitr3=irfitr3-1
        irfitr2=nrfitr2

!...................................................................
!     Extend rays
!...................................................................

        initrf=1
        nraypts1=nraypts1+nrfstep2
        nraypts2=nraypts2+nrfstep2
        nraypts3=nraypts3+nrfstep2
        call urfwrite
        call urfrays(initrf,nraypts1,nraypts2,nraypts3)
                write(*,*)'call urfread: n, nurf,irfitr2,irfitr3=' &
               ,n,nurf,irfitr2,irfitr3
        call urfread

!...................................................................
!     Obtain re-segregation by flux surface of ray data
!...................................................................

        call urffflx

!...................................................................
!     Extend packed data and do damping
!...................................................................

        call urfpack
        call urfdamp0(0,iopt)

      elseif(irfitr3.ge.0 .and. nrfitr3.gt.0 .and. irfitr2.gt.0)  then
        irfitr2=irfitr2-1

!...................................................................
!     do damping
!...................................................................

        call urfdamp0(0,iopt)
      endif


!..................................................................
!     Call the coefficient generator.
!     Radial surfaces: setup0%lrz (and not lrors)
!     Generate B,C,E and F
!..................................................................
      call cpu_time(t_urf1) !-.-.-.-.-.-.-.-.-.-.-.-.-.
      call urfb0  ! YuP-110222: Now includes all lr_ internally
      call cpu_time(t_urf2) !-.-.-.-.-.-.-.-.-.-.-.-.-.

#ifdef __MPI
      if(mpirank.eq.0) then
#endif
      WRITE(*,'(a,f10.3)')'tcpu_URFB0=',t_urf2-t_urf1
#ifdef __MPI
      endif  ! for if(mpirank.eq.***)
#endif

      imprf=1


      do 50 ll=1,setup0%lrz
!..................................................................
!     Read flux surface dependent data from disk.
!     for midplane parameter values
!..................................................................
        call tdnflxs(lmdpln(ll))
!%os
        if (indxlr_ .ne. ll) then
          print *,' error in defining l_,indxlr_,... in urfchief'
          stop 'urfchief: error in defining l_,indxlr_,...'
        endif
!%os
!..................................................................
!     Plot out diffusion coefficients
!..................................................................
        if (setup0%noplots.ne."enabled1") call urfbplt !only when plturfb="enabled"
                                                !or  plturfb='color'
 50   continue

      nurf=nurf+1

      return
      end subroutine urfchief


end module urfchief_mod
