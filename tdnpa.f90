module tdnpa_mod

  !---BEGIN USE

  use bcast_mod, only : bcast
  use tdfinterp_mod, only : tdfinterp
  use tdnpalam_mod, only : lam

  !---END USE

!
!

contains

      subroutine tdnpa(vel,pitch,ibinn,thpol,sigmanpa,emitnpa)
      use param_mod
      use comm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!vt....Modified from SRX routine tdsxr for NPA calculations.......
!vt....V. Tang, 9/29/05...................
!vt....ouputs emitnpa, which is the neutral emittance as a function
!vt....of energy before attenuation
!vt....needs sigmanpa, which is calculated in tdsetnpa.f
!vt....Version 1.0, NPA Mod for CQL3D

!BH100425: Adjusted.

!.......................................................................
!
!     This routine computes the CHARGE EXCHANGE (CX) contribution to
!     the NPA spectrum from H+ - H, or D+ -D collisions.
!
!     vel(1:nen_npa) are velocities (cms/sec) corresponding to the
!     the diagnostic energy grid.
!
!     pitch is local pitch angle w.r.t. B-field, towards the diagnostic.
!
!     thpol is the local poloidal angle (representing the point of
!     intersection of the viewing angle and the relevant flux
!     surface).
!
!     sigmanpa(1:nen_npa,) cross-section (cm**2) is input (except
!     sigmanpa(1:nen_npa,5) for e-recombination, which is calculated below.)
!
!     emitnpa(1:nen_npa) is neutral particle emittance
!     (#/(sec*cm**3*steradian*eV) headed towards the diagnostic
!     for the particular bin (pitch,rho,thpol) contribution,
!     where rho=rya(ibinn).
!
!.......................................................................


      real(c_double),dimension(nen_npa)::vel
      real(c_double),dimension(nen_npa,npaproc)::sigmanpa,emitnpa

      real(c_double),dimension(nen_npa)::EH,E_old   !debug purposes

!     Calculate the emittance array: #/(sec*vol*steradian*eV)
!     (Need to multiply by ergtkev/1e3 to convert from cgs (/erg) to /eV.
      call bcast(emitnpa,zero,nen_npa*npaproc)

      if (npa_process(1).eq.'cxh') then
      do j=1,nen_npa
         call tdfinterp(vel(j),pitch,rya(ibinn),thpol,fdist)
         emitnpa(j,1)= (ergtkev*1.d-3*fdist)/(vnorm3*fmass(1)) &
                       *(sigmanpa(j,1)*vel(j)**2)
      enddo
      endif  ! On npa_process(1)


      if (npaproc.ge.2 .and. npa_process(2).eq.'cxb4') then
      do j=1,nen_npa
         call tdfinterp(vel(j),pitch,rya(ibinn),thpol,fdist)
         emitnpa(j,2)= (ergtkev*1.d-3*fdist)/(vnorm3*fmass(1)) &
                       *(sigmanpa(j,2)*vel(j)**2)
      enddo
      endif  ! On npa_process(2)


      if (npaproc.ge.5 .and. npa_process(5).eq.'radrecom') then
!        Using Bethe-Salpeter recombination cross-section formula.
!        See web, H. Poth and A. Wolf, CERN-EP/82-189
!        (25 November 1982), 198301145.pdf, Eq. (1).
!        sigma_n=1.96*pi**2*alpha*lambda_e**2/(n*(1+n**2*Te/Ry)*Te/Ry)
!        See formula and NRL table for defn and values of constants.
!        Take Te=electron energy=(3/2) electron temperature
!        Relative electron velocity = sqrt(elect temp/me) = vth
!
!        110806: Updated to improved expression from Ian Hutchinson book,
!                formula 6.3.5  [Aaron Bader and BH].

      do j=1,nen_npa
!BH110725         t_elec=temp(kelecm,ibinn)*1.e3/13.6  ! normalized to electron
!BH110806         t_elec=1.5*temp(kelecm,ibinn)*1.e3/13.6  ! normalized to electron
                                              ! ionization energy.
!BH110806         sigmanpa(j,5)=2.105d-22*(1./((1.+t_elec)*t_elec)
!BH110806     +                          +1./(2.*(1.+4.*t_elec)*t_elec))
         zz=1 !for H/D
         chidt=ZZ*13.6d0/(1000.*temp(kelecm,ibinn))
         sigmanpa(j,5)=5.2d-14*zz*(chidt**1.5d0*log(1.d0/chidt)+ &
                          (chidt/4.d0)**1.5*log(4.d0/chidt)) &
                          /vth(kelecm,ibinn) !divide for sigma expression
         call tdfinterp(vel(j),pitch,rya(ibinn),thpol,fdist)

         emitnpa(j,5)= (ergtkev*1.d-3*fdist)/(vnorm3*fmass(1)) &
                       *(sigmanpa(j,5)*vth(kelecm,ibinn)*vel(j))

         zz=1 !for H/D/T
      enddo
      endif  ! On npa_process(5)

      return
      end
end module tdnpa_mod
