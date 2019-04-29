module ainspec_mod

  !---BEGIN USE

  !---END USE

!
!

contains

      subroutine ainspec
      use param_mod
      use comm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!.......................................................................
!     Determine
!     kelecg = the index of the general species electrons (0,if none).
!     kelecm = the background Maxwellian species electrons (0,if none).
!     kionn =  1st ionic species (be it general or background species).
!     kelec = kelecg (if kelecg>0) or  kelec = kelecm (if kelecg=0)
!     niong = number of general species ions
!     kiong(1:niong) = index of general ion specie (otherwise 0 for 1:ngen)
!     kiongg(k=1:ngen) = k if an ions species, otherwise 0.
!     nionm = number of ion background Maxwellian species.
!     kionm(1:nionm) = Maxwellian ion specie indices (otherwise 0).
!
!     There must be at least one electron and
!     one ionic species present.
!
!    From cqlinput_help:
!    In the following, the index (k) refers to the species.
!    If 1 .le. k .le. ngen, then k is a general (time-advanced) species.
!    k.gt.ngen means species k is background Maxwellian.
!    In k order, the general (non-Maxwellian time advanced species)
!    come before the background species. Electrons and
!    ions can be mixed, and any species can be represented
!    simultaneously as a general and a background species.
!    If there is more than one general species ions, they must
!    not be separated in k by the electron species
!    (i.e., put general electrons at front or back of
!    the general ion species).
!    If there is more than one Maxwellian ion, then they must
!    not be separated in k by the electron species.   If more than
!    one of the ion species has the same charge number bnumb(), (e.g.,
!    H+, D+, and/or T+) then they should be placed successively at the
!    beginning of the Maxwellian ion species.

!.......................................................................

!MPIINSERT_INCLUDE

      ntotal=ngen+nmax

      kelecg=0
!DIR$ NEXTSCALAR
      do 1000 k=1,ngen
        if (fmass(k) .lt. 1.d-27) then
          kelecg=k
          goto 1001
        endif
 1000 continue
 1001 continue
      kelecm=0
!DIR$ NEXTSCALAR
      do 1002 k=ngen+1,ntotal
        if (fmass(k) .lt. 1.d-27) then
          kelecm=k
          goto 1003
        endif
 1002 continue
 1003 continue
      if (kelecg.gt.0) then
        kelec=kelecg
      else
        kelec=kelecm
      endif
!BH180908      if (kelec.eq.0) call diagwrng(9)
      if (kelec.eq.0) then
!MPIINSERT_IF_RANK_EQ_0
         WRITE(*,*)
         WRITE(*,*) 'WARNING: Unphysical plasma, only one species.'
         WRITE(*,*)
!MPIINSERT_ENDIF_RANK
      endif


      kionn=0
!DIR$ NEXTSCALAR
      do 1005 k=1,ntotal
        if (fmass(k).gt.1.d-27) then
          kionn=k
          goto 1006
        endif
 1005 continue
 1006 continue
      do 1007 k=1,ngen
        kspeci(2,k)="general"
 1007 continue
      do 1008 k=ngen+1,ntotal
        kspeci(2,k)="maxwell"
 1008 continue

      niong=0
      do 1009  k=1,ngena
        kiong(k)=0
        kiongg(k)=0
 1009 continue
!     Electron mass is 9.1e-28 grams. Use this to distinguish ions.
      do 1010  k=1,ngen
        if (fmass(k).gt.1.d-27 .and. kspeci(2,k).eq."general") then
          niong=niong+1
          kiong(niong)=k
          kiongg(k)=k
        endif
 1010 continue
      nionm=0
      do 1011 k=1,nmaxa
        kionm(k)=0
 1011 continue
      do 1012 k=ngen+1,ntotal
        if (fmass(k).gt.1.d-27 .and. kspeci(2,k).eq."maxwell") then
          nionm=nionm+1
          kionm(nionm)=k
        endif
 1012 continue

!BHTemp      if (kionn.eq.0) call diagwrng(9)
      if (kionn.eq.0) then
!MPIINSERT_IF_RANK_EQ_0
         WRITE(*,*)
         WRITE(*,*) 'WARNING: Unphysical plasma, only one species.'
         WRITE(*,*)
!MPIINSERT_ENDIF_RANK
      endif

!MPIINSERT_IF_RANK_EQ_0
      WRITE(*,*)
      WRITE(*,*)'ainspec: ngen,nmax,ntotal', &
                          ngen,nmax,ntotal
      WRITE(*,*)'ainspec: kelecg,kelecm,kelec', &
                          kelecg,kelecm,kelec
      WRITE(*,*)'ainspec: niong,kiong(1:niong)',niong,kiong(1:niong)
      WRITE(*,*)'ainspec: nionm,kionm(1:nionm)',nionm,kionm(1:nionm)
      WRITE(*,*)'ainspec: nionm/kionm may be adjusted upwards when'
      WRITE(*,*)'ainspec:   iprozeff.ne."disabled"; see cqlinput_help'
      WRITE(*,*)
!MPIINSERT_ENDIF_RANK

      return
      end
end module ainspec_mod
