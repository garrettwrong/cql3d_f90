module urfread__mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use urfwrong_mod, only : urfwrong

  !---END USE

!
!

contains

      subroutine urfread_(krf,i_)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)


!..................................................................
!     Up to several files, corresponding to ray data from
!     different wave modes/types, may be read.
!     For each mode considered (krf):
!     Read total number of rays
!     Read harmonic for ql diffusion (0, for Landau/TTMP)
!     Read in the wave antenna frequency (Hz).
!..................................................................
      io1=0
      io2=0
      read(i_,1) nray(krf),nharm(krf),freqcy(krf)
      WRITE(*,*)'urfread_:nray,nharm,freqcy=', &
                 nray(krf),nharm(krf),freqcy(krf)
      omega(krf)=2.*pi*freqcy(krf)

!..................................................................
!     Send message if nray.gt.nrayn
!..................................................................

      if (nray(krf).gt.nrayn) then
        nr=nrayn
        WRITE(*,200) nray(krf),nr
        nray(krf)=nrayn
      endif
 200  format("Ray tracing code provides more rays than CQL parameter",/, &
        "nrayn is equipped to handle. Recompile CQL with nrayn=",i5,/, &
        "CQL will proceed with the calculation with nray reset to",i5)

!..................................................................
!     All units are cgs.
!
!     (Referring to following read statements:)
!     Read in ray data, ray by ray.  For each ray read:
!     nrayelt(iray,k)=number of ray elements (i.e., data) for ray iray,k
!     The next nine integers in the read are the final states of the
!     mode, number of refl., integ.steps, number of electron and ion
!     power calcs and the restart status of the ray.
!     (For EC, can be set equal to 0)
!     (Also can be 0, if no restart needed).
!     For EC, next 4 real numbers may be zero).
!     sxxrt(iray,k) - final value of integ. variable
!     skpsi(iray,k) - canonical n_psi
!     skth(iray,k) - canonical n_th
!     skphi - canonical n_phi
!
!     The following data is given for is=1,nrayelt(iray), for each ray:
!     ws - poloidal plane projection of distance along ray.
!     seikon - eikon    (wave length count.  But 0.0's OK).
!     spsi - rho   (any radial-like variable.  Used for plotting).
!     wr - R    (major radius)
!     wphi - toroidal angle phi (radians, used for plotting).
!     wz - Z     (height above midplane)
!     wnpar - n_parallel (sign changed below,
!                         if partner="lsc" or bsign.lt.0.)
!     wnper - n_perp
!     delpwr - power flowing in the ray channel (ergs/sec).
!     (For EC, Initial value, then 0.'s, is good).
!     (Can be modified with namelist variable pwrscale.)
!     sdpwr - dpwri (power lost due to ion absorption, if included).
!     (EC:  0.0's). (Presently not used for calculations in the code).
!     wdnpar - delta n_parallel (n_parallel width locally)
!     (In toroidal system, this is chosen so that
!     the delta n_phi x R, is constant along the ray).
!     (Can be modified with namelist variable wdscale).
!     cwexde - the complex*16 polarization E_x/E
!     cweyde - the complex*16 polarization E_y/E
!     cwezde - the complex*16 polarization E_z/E
!     fluxn - the normalized flux (electric field=1)
!       .5*[B.B^* + E.d(omega(krf) K_h)/d(omega(krf)).E^*]*
!                                                  (v_group_pol/clight)
!       with abs(E)=1., n x n x B = n x E.
!       v_group_pol is group vel in poloidal plane,
!       K_h is the Hermitian part of the dielectric tensor.
!       cf. Stix, Waves in Plasmas (1992), p. 70-77.
!     sbtot - total B field (gauss)
!     sene - n_e(/cc)
!     salphac - perp.(i.e., appox in poloidal plane ) coll.
!     energy absorption coefficient (/cm).
!     EC:  0.0's or other non-cql3d damping, is ok
!     salphal - linear energy perp absorption coeff. (rf:on ions) (/cm).
!     (EC: 0.0's, or a linear cyclotron rate for comparison with cql3d).
!..................................................................


      iformat1=0   !To handle format 3 or 4 for integers.
      iformat2=0   !BH091116: To handle complex cnperp variable in
                   !toray data, added into toray 1.8 (and perhaps
                   !earlier).

      do 10 iray=1,nray(krf)
!        write(*,*)'urfread_:Here0, iray, iformat1',iray,iformat1
        if (iformat1.eq.1) go to 7
        read(i_,3,iostat=io1) nrayelt(iray,krf),jslofas(iray,krf), &
              nurefls(iray,krf),keiks(iray,krf),jpes(iray,krf), &
              jpis(iray,krf),istarts(iray,krf), &
              iprmt5(iray,krf),jhlfs(iray,krf)
        if (iray.eq.1 .and. io1.ne.0) go to 5
        if (iray.ge.2 .and. io1.ne.0) then
           iformat2=1
           go to 100
        endif
!        write(*,*)'urfread_:Here0.1, iray, iformat2',iray,iformat2

!        write(*,3) nrayelt(iray,krf),jslofas(iray,krf),
!     +        nurefls(iray,krf),keiks(iray,krf),jpes(iray,krf),
!     +        jpis(iray,krf),istarts(iray,krf),
!     +        iprmt5(iray,krf),jhlfs(iray,krf)
!        write(*,*)'urfread_:Here1, iray, iformat1',iray,iformat1
        go to 6

!     Alternate format for this read:
 5      iformat1=1
        rewind i_
        read(i_,1) nray(krf),nharm(krf),freqcy(krf)

 7      read(i_,4,iostat=io2) nrayelt(iray,krf),jslofas(iray,krf), &
              nurefls(iray,krf),keiks(iray,krf),jpes(iray,krf), &
              jpis(iray,krf),istarts(iray,krf), &
              iprmt5(iray,krf),jhlfs(iray,krf)

!BH091116  We assume that if can't find these integers in the
!BH091116  standard place in the file, then a new toray file
!BH091116  with cnper is being read.
        if (io2.ne.0) then
           iformat2=1
           go to 100
        endif

!        write(*,3) nrayelt(iray,krf),jslofas(iray,krf),
!     +        nurefls(iray,krf),keiks(iray,krf),jpes(iray,krf),
!     +        jpis(iray,krf),istarts(iray,krf),
!     +        iprmt5(iray,krf),jhlfs(iray,krf)
!        write(*,*)'urfread_:nray(krf),krf,iray,nrayelt(iray,krf)=',
!     +                                krf,iray,nrayelt(iray,krf)

 6      read(i_,2) sxxrt(iray,krf),skpsi(iray,krf),skth(iray,krf), &
          skphi(iray,krf)

!..................................................................
!     Abort if the ray tracing code provides more ray elements than
!     CQL is equipped to handle.
!     YuP 101122: Should never happen now - nrayelts
!     is determined in urfsetup from reading the rays' data file.
!..................................................................

        if (nrayelt(iray,krf).gt.nrayelts) then
          call urfwrong(6)
        endif
        if (nrayelt(iray,krf).eq.0) then
           WRITE(*,*) 'urfread_: nrayelt(iray,krf)=0, iray,krf= ', &
                       iray,krf
           go to 10
        endif
        ! for iformat2=0:
        read(i_,2) (ws(is,iray,krf),is=1,nrayelt(iray,krf))      ! #1
        read(i_,2) (seikon(is,iray,krf),is=1,nrayelt(iray,krf))  ! #2
        read(i_,2) (spsi(is,iray,krf),is=1,nrayelt(iray,krf))    ! #3
        read(i_,2) (wr(is,iray,krf),is=1,nrayelt(iray,krf))      ! #4
        read(i_,2) (wphi(is,iray,krf),is=1,nrayelt(iray,krf))    ! #5
        read(i_,2) (wz(is,iray,krf),is=1,nrayelt(iray,krf))      ! #6
        read(i_,2) (wnpar(is,iray,krf),is=1,nrayelt(iray,krf))   ! #7
        read(i_,2) (wnper(is,iray,krf),is=1,nrayelt(iray,krf))   ! #8
        read(i_,2) (delpwr(is,iray,krf),is=1,nrayelt(iray,krf))  ! #9
        read(i_,2) (sdpwr(is,iray,krf),is=1,nrayelt(iray,krf))   ! #10
        read(i_,2) (wdnpar(is,iray,krf),is=1,nrayelt(iray,krf))  ! #11
        read(i_,2) (cwexde(is,iray,krf),is=1,nrayelt(iray,krf))  ! #12 complex
        read(i_,2) (cweyde(is,iray,krf),is=1,nrayelt(iray,krf))  ! #13 complex
        read(i_,2) (cwezde(is,iray,krf),is=1,nrayelt(iray,krf))  ! #14 complex
        read(i_,2) (fluxn(is,iray,krf),is=1,nrayelt(iray,krf))   ! #15
        read(i_,2) (sbtot(is,iray,krf),is=1,nrayelt(iray,krf))   ! #16
        read(i_,2) (sene(is,iray,krf),is=1,nrayelt(iray,krf))    ! #17
        read(i_,2) (salphac(is,iray,krf),is=1,nrayelt(iray,krf)) ! #18
        read(i_,2) (salphal(is,iray,krf),is=1,nrayelt(iray,krf)) ! #19

        do  is=1,nrayelt(iray,krf)
           wdnpar(is,iray,krf)=wdscale(krf)*abs(wdnpar(is,iray,krf))
        enddo

!     ON NEGATIVE SBTOT:          (BobH prompted by Lin-Liu, 981016)
!     =================
!     CQL3D RF routines have not been set up to use negative sbtot.
!     (btor and f=R*B_phi are also made positive in subroutine equilib).
!     (For negative toroidal field in the eqdsk, bsign must be < 0.
!      in the namelist input for CQL3D to successfully execute).
!
!     n_parallel/k_parallel in CQL3D is understood to be positive if it
!       has a component in the positive (counter-clockwise from above)
!       toroidal direction.  k_parallel may be positive or negative.
!       (Maybe this should be re-examined: BobH, 991014).
!
!     In TORAY: k_parallel is defined according to the usual
!              plasma physics convention positive in the
!              direction of the magnetic field.
!              Polarizations are also defined w.r.t. the dirn of B.
!              (c.f., Stix, Waves in Plasmas (1992), p. 10. Also p499).
!              We observe in TORAY, when k_parallel ==> -k_parallel,
!              Mazzucato polarizations (cqldat=.true.) transform as
!                   cwexde ==> -cwexde
!                   cweyde ==> -cweyde
!                   cwezde ==> +cwezde.
!              This k_parallel-transformation has the effect of
!              reversing the parallel velocity direction for resonance
!              in CQL3D but keeps the magnitude of
!              damping constant (in accord with symmetry under
!              change of dirn w.r.t. B).
!
!    In GENRAY:k_parallel is also defined to be positive
!              in the direction of the magnetic field.
!              But sbtot is taken to be the magnitude of the
!              total magnetic field, with sign of vacuum magnetic fld.
!    Changed as of 040816:     total magnetic field (always positive).
!
!    Therefore, for the case of negative vacuum toroidal field
!      in the eqdsk, we reverse the sign of n_parallel (wnpar)
!      as well as in cwexde, and cweyde.  This maintains the CQL3D
!      defn of positive n_parallel as having a positive component
!      in the positive toroidal direction.  This is also the direction
!      of positive v_parallel.  Reversing sign of cwexde and cweyde
!      keeps the the magnitude of damping constant with the sign
!      change of n_parallel.
!    Negative sbtot values (as from TORAY or GENRAY for negative
!      vacuum toroidal B-field) are sign-reversed.
!
!
!     Change sign of sbtot(if neg)  and wnpar, if bsign.lt.0.
!     (see cqlinput_help)
!     However, GENRAY and TORAY can have different conventions as above)
!     for the sign of sbtot.  cql3d thus detects sign of first
!     element of sbtot and saves it in bsign1, so the sign convention
!     of the input file can be maintained for the output ray data file
!     from cql3d.
        if (eqsource.eq."eqdsk" .and. bsign.lt.0.) then
           sbsign=sign(one,bsign)
           if (iray.eq.1) then
              bsign1(krf)=sign(one,sbtot(1,iray,krf))
              WRITE(*,*)'urfread_:eqsource,bsign,bsign1 = ', &
                   eqsource,bsign,bsign1
           endif
           do is=1,nrayelt(iray,krf)
!              if(is.eq.1) write(*,*)
!     1                  'urfread_:is,iray,krf, sbtot((is,iray,krf) = ',
!     1                  is,iray,krf,sbtot(is,iray,krf)
              sbtot(is,iray,krf)=bsign1(krf)*sbtot(is,iray,krf)
              wnpar(is,iray,krf)=sbsign*wnpar(is,iray,krf)
              cwexde(is,iray,krf)=sbsign*cwexde(is,iray,krf)
              cweyde(is,iray,krf)=sbsign*cweyde(is,iray,krf)
           enddo
        endif

!     fluxn is renormalized to be as in Stix or Bekefi:
        do is=1,nrayelt(iray,krf)
           fluxn(is,iray,krf)=fluxn(is,iray,krf)*clight/(8.*pi)
        enddo

!     shift z-position of ray data, if eqdsk has been
!     vertically shifted in subroutine equilib:
        if (zshift.ne.zero) then
           do is=1,nrayelt(iray,krf)
              wz(is,iray,krf)=wz(is,iray,krf)-zshift
           enddo
        endif

 10   continue


!BH091116:  Alternative ray data read, adjusting for additional
!BH091116:  complex cnper variable in ray data.  Use cwexde memory
!BH091116:  to read it, then read cwexde in over it.
!BH091116:  iformat1, iformat2 have been set above

 100  if (iformat2.eq.1) then

      rewind i_
      read(i_,1) nray(krf),nharm(krf),freqcy(krf)

      do 110 iray=1,nray(krf)
!        write(*,*)'urfread_:Here0.2, iray, iformat1',iray,iformat1
        if (iformat1.eq.1) go to 107
        read(i_,3,iostat=io1) nrayelt(iray,krf),jslofas(iray,krf), &
              nurefls(iray,krf),keiks(iray,krf),jpes(iray,krf), &
              jpis(iray,krf),istarts(iray,krf), &
              iprmt5(iray,krf),jhlfs(iray,krf)

!        write(*,3) nrayelt(iray,krf),jslofas(iray,krf),
!     +        nurefls(iray,krf),keiks(iray,krf),jpes(iray,krf),
!     +        jpis(iray,krf),istarts(iray,krf),
!     +        iprmt5(iray,krf),jhlfs(iray,krf)
!        write(*,*)'urfread_:Here1.2,iray,io1,iformat1',iray,io1,iformat1
        go to 106

 107    read(i_,4,iostat=io2) nrayelt(iray,krf),jslofas(iray,krf), &
              nurefls(iray,krf),keiks(iray,krf),jpes(iray,krf), &
              jpis(iray,krf),istarts(iray,krf), &
              iprmt5(iray,krf),jhlfs(iray,krf)

!        write(*,3) nrayelt(iray,krf),jslofas(iray,krf),
!     +        nurefls(iray,krf),keiks(iray,krf),jpes(iray,krf),
!     +        jpis(iray,krf),istarts(iray,krf),
!     +        iprmt5(iray,krf),jhlfs(iray,krf)
!        write(*,*)'urfread_:nray(krf),krf,iray,nrayelt(iray,krf)=',
!     +                                krf,iray,nrayelt(iray,krf)

 106    read(i_,2) sxxrt(iray,krf),skpsi(iray,krf),skth(iray,krf), &
          skphi(iray,krf)

!..................................................................
!     Abort if the ray tracing code provides more ray elements than
!     CQL is equipped to handle.
!     YuP 101122: Should never happen now - nrayelts
!     is determined in urfsetup from reading the rays' data file.
!..................................................................

        if (nrayelt(iray,krf).gt.nrayelts) then
          call urfwrong(6)
        endif
        if (nrayelt(iray,krf).eq.0) then
           WRITE(*,*) 'urfread_: nrayelt(iray,krf)=0, iray,krf= ', &
                       iray,krf
           go to 110
        endif
        ! for iformat2=1:
        read(i_,2) (ws(is,iray,krf),is=1,nrayelt(iray,krf))      ! #1
        read(i_,2) (seikon(is,iray,krf),is=1,nrayelt(iray,krf))  ! #2
        read(i_,2) (spsi(is,iray,krf),is=1,nrayelt(iray,krf))    ! #3
        read(i_,2) (wr(is,iray,krf),is=1,nrayelt(iray,krf))      ! #4
        read(i_,2) (wphi(is,iray,krf),is=1,nrayelt(iray,krf))    ! #5
        read(i_,2) (wz(is,iray,krf),is=1,nrayelt(iray,krf))      ! #6
        read(i_,2) (wnpar(is,iray,krf),is=1,nrayelt(iray,krf))   ! #7
        read(i_,2) (wnper(is,iray,krf),is=1,nrayelt(iray,krf))   ! #8
        read(i_,2) (delpwr(is,iray,krf),is=1,nrayelt(iray,krf))  ! #9
        read(i_,2) (sdpwr(is,iray,krf),is=1,nrayelt(iray,krf))   ! #10
        read(i_,2) (wdnpar(is,iray,krf),is=1,nrayelt(iray,krf))  ! #11

!BH091116:  Next line is dummy read of extra toray data (cnper).
        read(i_,2) (cwexde(is,iray,krf),is=1,nrayelt(iray,krf))  ! #12 complex

        read(i_,2) (cwexde(is,iray,krf),is=1,nrayelt(iray,krf))  ! #13 complex
        read(i_,2) (cweyde(is,iray,krf),is=1,nrayelt(iray,krf))  ! #14 complex
        read(i_,2) (cwezde(is,iray,krf),is=1,nrayelt(iray,krf))  ! #15 complex
        read(i_,2) (fluxn(is,iray,krf),is=1,nrayelt(iray,krf))   ! #16
        read(i_,2) (sbtot(is,iray,krf),is=1,nrayelt(iray,krf))   ! #17
        read(i_,2) (sene(is,iray,krf),is=1,nrayelt(iray,krf))    ! #18
        read(i_,2) (salphac(is,iray,krf),is=1,nrayelt(iray,krf)) ! #19
        read(i_,2) (salphal(is,iray,krf),is=1,nrayelt(iray,krf)) ! #20

!BH091116:  Next line is dummy read of extra toray data (salphal_vg).
        read(i_,2) (urftmp(is),is=1,nrayelt(iray,krf))           ! #21


        do  is=1,nrayelt(iray,krf)
           wdnpar(is,iray,krf)=wdscale(krf)*abs(wdnpar(is,iray,krf))
        enddo

!     Change sign of sbtot(if neg)  and wnpar:  See above comments.
        if (eqsource.eq."eqdsk" .and. bsign.lt.0.) then
           sbsign=sign(one,bsign)
           if (iray.eq.1) then
              bsign1(krf)=sign(one,sbtot(1,iray,krf))
              WRITE(*,*)'urfread_:eqsource,bsign,bsign1 = ', &
                   eqsource,bsign,bsign1
           endif
           do is=1,nrayelt(iray,krf)
!              if(is.eq.1) write(*,*)
!     1                  'urfread_:is,iray,krf, sbtot((is,iray,krf) = ',
!     1                  is,iray,krf,sbtot(is,iray,krf)
              sbtot(is,iray,krf)=bsign1(krf)*sbtot(is,iray,krf)
              wnpar(is,iray,krf)=sbsign*wnpar(is,iray,krf)
              cwexde(is,iray,krf)=sbsign*cwexde(is,iray,krf)
              cweyde(is,iray,krf)=sbsign*cweyde(is,iray,krf)
           enddo
        endif

!     fluxn is renormalized to be as in Stix or Bekefi:
        do is=1,nrayelt(iray,krf)
           fluxn(is,iray,krf)=fluxn(is,iray,krf)*clight/(8.*pi)
        enddo

!     shift z-position of ray data, if eqdsk has been
!     vertically shifted in subroutine equilib:
        if (zshift.ne.zero) then
           do is=1,nrayelt(iray,krf)
              wz(is,iray,krf)=wz(is,iray,krf)-zshift
           enddo
        endif

 110   continue

      endif

      WRITE(*,*)'urfread_:iformat1,iformat2,io1,io2=', &
                          iformat1,iformat2,io1,io2

 1    format(2i5,1pe16.9)
 2    format(5(1pe16.9))
 3    format(9i6)
 4    format(9i5)
!cc 3    format (9i5) ! YuP:  as in GENRAY
!cc 4    format (9(' ',i4)) ! YuP:  as inGENRAY

      return
      end subroutine urfread_





!=======================================================================
!=======================================================================
      subroutine urfread_i(krf,i_)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

      ! working array for dummy read
      real(c_double), allocatable, dimension(:) :: urfwk

!===> YuP 101122: A simplified version of urfread_ above.
!     Called by urfsetup
!     to read the values of nray(krf) and nrayelt(iray,krf)

      read(i_,1) nray(krf),nharm(krf),freqcy(krf)

      iformat1=0   !To handle format 3 or 4 for integers.
      iformat2=0   !BH091116: To handle complex cnperp variable in
                   !toray data, added into toray 1.8 (and perhaps
                   !earlier).
      io1=0
      io2=0

      do 10 iray=1,nray(krf)
        if (iformat1.eq.1) go to 7
        read(i_,3,iostat=io1) nrayelt(iray,krf),jslofas(iray,krf), &
              nurefls(iray,krf),keiks(iray,krf),jpes(iray,krf), &
              jpis(iray,krf),istarts(iray,krf), &
              iprmt5(iray,krf),jhlfs(iray,krf)
        if (iray.eq.1 .and. io1.ne.0) go to 5
        if (iray.ge.2 .and. io1.ne.0) then
           iformat2=1
           go to 100
        endif
        go to 6

!     Alternate format for this read:
 5      iformat1=1
        rewind i_
        read(i_,1) nray(krf),nharm(krf),freqcy(krf)

 7      read(i_,4,iostat=io2) nrayelt(iray,krf),jslofas(iray,krf), &
              nurefls(iray,krf),keiks(iray,krf),jpes(iray,krf), &
              jpis(iray,krf),istarts(iray,krf), &
              iprmt5(iray,krf),jhlfs(iray,krf)

!BH091116  We assume that if can't find these integers in the
!BH091116  standard place in the file, then a new toray file
!BH091116  with cnper is being read.
        if (io2.ne.0) then
           iformat2=1
           go to 100
        endif

 6      read(i_,2) sxxrt(iray,krf),skpsi(iray,krf),skth(iray,krf), &
          skphi(iray,krf)

        if (nrayelt(iray,krf).eq.0) then
           WRITE(*,*)'urfread_i in do10: nrayelt(iray,krf)=0,iray,krf=', &
                       iray,krf
           go to 10
        endif

        nrayelti= nrayelt(iray,krf)
        allocate(urfwk(nrayelti),STAT=istat) ! real(c_double)
        ! dummy read for iformat2=0:
        read(i_,2) (urfwk(is),is=1,nrayelti) ! #1
        read(i_,2) (urfwk(is),is=1,nrayelti) ! #2
        read(i_,2) (urfwk(is),is=1,nrayelti) ! #3
        read(i_,2) (urfwk(is),is=1,nrayelti) ! #4
        read(i_,2) (urfwk(is),is=1,nrayelti) ! #5
        read(i_,2) (urfwk(is),is=1,nrayelti) ! #6
        read(i_,2) (urfwk(is),is=1,nrayelti) ! #7
        read(i_,2) (urfwk(is),is=1,nrayelti) ! #8
        read(i_,2) (urfwk(is),is=1,nrayelti) ! #9
        read(i_,2) (urfwk(is),is=1,nrayelti) ! #10
        read(i_,2) (urfwk(is),is=1,nrayelti) ! #11
        read(i_,2) (urfwk(is),urfwk(is),is=1,nrayelti) ! #12 complex
        read(i_,2) (urfwk(is),urfwk(is),is=1,nrayelti) ! #13 complex
        read(i_,2) (urfwk(is),urfwk(is),is=1,nrayelti) ! #14 complex
        read(i_,2) (urfwk(is),is=1,nrayelti) ! #15
        read(i_,2) (urfwk(is),is=1,nrayelti) ! #16
        read(i_,2) (urfwk(is),is=1,nrayelti) ! #17
        read(i_,2) (urfwk(is),is=1,nrayelti) ! #18
        read(i_,2) (urfwk(is),is=1,nrayelti) ! #19
        deallocate(urfwk,STAT=istat)

 10   continue ! iray=1,nray(krf)


!BH091116:  Alternative ray data read, adjusting for additional
!BH091116:  complex cnper variable in ray data.  Use cwexde memory
!BH091116:  to read it, then read cwexde in over it.
!BH091116:  iformat1, iformat2 have been set above

 100  if (iformat2.eq.1) then

      rewind i_
      read(i_,1) nray(krf),nharm(krf),freqcy(krf)

      do 110 iray=1,nray(krf)
        if (iformat1.eq.1) go to 107
        read(i_,3,iostat=io1) nrayelt(iray,krf),jslofas(iray,krf), &
              nurefls(iray,krf),keiks(iray,krf),jpes(iray,krf), &
              jpis(iray,krf),istarts(iray,krf), &
              iprmt5(iray,krf),jhlfs(iray,krf)

        go to 106

 107    read(i_,4,iostat=io2) nrayelt(iray,krf),jslofas(iray,krf), &
              nurefls(iray,krf),keiks(iray,krf),jpes(iray,krf), &
              jpis(iray,krf),istarts(iray,krf), &
              iprmt5(iray,krf),jhlfs(iray,krf)

 106    read(i_,2) sxxrt(iray,krf),skpsi(iray,krf),skth(iray,krf), &
          skphi(iray,krf)

        if (nrayelt(iray,krf).eq.0) then
          WRITE(*,*)'urfread_i in do110: nrayelt(iray,krf)=0,iray,krf=', &
                       iray,krf
          go to 110
        endif

        nrayelti= nrayelt(iray,krf)
        allocate(urfwk(nrayelti),STAT=istat) ! real(c_double)
        ! dummy read for iformat2=1:
        read(i_,2) (urfwk(is),is=1,nrayelti) ! #1
        read(i_,2) (urfwk(is),is=1,nrayelti) ! #2
        read(i_,2) (urfwk(is),is=1,nrayelti) ! #3
        read(i_,2) (urfwk(is),is=1,nrayelti) ! #4
        read(i_,2) (urfwk(is),is=1,nrayelti) ! #5
        read(i_,2) (urfwk(is),is=1,nrayelti) ! #6
        read(i_,2) (urfwk(is),is=1,nrayelti) ! #7
        read(i_,2) (urfwk(is),is=1,nrayelti) ! #8
        read(i_,2) (urfwk(is),is=1,nrayelti) ! #9
        read(i_,2) (urfwk(is),is=1,nrayelti) ! #10
        read(i_,2) (urfwk(is),is=1,nrayelti) ! #11
        read(i_,2) (urfwk(is),urfwk(is),is=1,nrayelti) ! #12 complex
        read(i_,2) (urfwk(is),urfwk(is),is=1,nrayelti) ! #13 complex
        read(i_,2) (urfwk(is),urfwk(is),is=1,nrayelti) ! #14 complex
        read(i_,2) (urfwk(is),urfwk(is),is=1,nrayelti) ! #15 complex
        read(i_,2) (urfwk(is),is=1,nrayelti) ! #16
        read(i_,2) (urfwk(is),is=1,nrayelti) ! #17
        read(i_,2) (urfwk(is),is=1,nrayelti) ! #18
        read(i_,2) (urfwk(is),is=1,nrayelti) ! #19
        read(i_,2) (urfwk(is),is=1,nrayelti) ! #20
        read(i_,2) (urfwk(is),is=1,nrayelti) ! #21
        deallocate(urfwk,STAT=istat)

 110  continue ! iray=1,nray(krf)  for iformat2 = 1

      endif ! iformat2=1

      WRITE(*,*)'urfread_i: iformat1,iformat2,io1,io2=', &
                            iformat1,iformat2,io1,io2

 1    format(2i5,1pe16.9)
 2    format(5(1pe16.9))
 3    format(9i6)
 4    format(9i5)
!cc 3    format (9i5) ! YuP:  as in GENRAY
!cc 4    format (9(' ',i4)) ! YuP:  as inGENRAY
!MG added 11/13/2017
      if(allocated(urfwk)) deallocate(urfwk)
!MG end added 11/13/2017
      return
      end subroutine urfread_i
      
end module urfread__mod
