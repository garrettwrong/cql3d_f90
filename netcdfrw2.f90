!
!
module netcdfrw2_mod
      use param_mod, only : nmodsa
      use bcast_mod, only : bcast
      use bcast_mod, only : ibcast
      !use pack21_mod, only : pack21
      !use pack21_mod, only : unpack21
      external pack21    ! XXXX ugh, you should fix these.  YuP: don't know what exactly to fix
      external unpack21  ! XXXX ugh, you should fix these.
      !integer, private ::  numrec(nmodsa) !-YuP-> added: as a function of krf

      save

contains
!-----------------------------------------------------------------------


      integer function length_char(string)
!     Returns length of string, ignoring trailing blanks,
!     using the fortran intrinsic len().
      implicit none
      character*(*) string ! input
      integer i ! local
      do i=len(string),1,-1
         if(string(i:i) .ne. ' ') goto 20
      enddo
 20   length_char=i
      return
      end function length_char


      subroutine netcdfrw2(kopt)
      use cqlconf_mod, only : setup0
      use advnce_mod !here: in netcdfrw2().
      use bcast_mod, only : bcast
      use cqlcomm_mod
      use diagentr_mod, only : gfu
      use param_mod
      use prppr_mod, only : prppr
      use r8subs_mod, only : dcopy
      use zcunix_mod, only : terp1
      implicit none
      save

      integer :: kopt
      integer :: i, ii, ielem, istat, iyy
      integer :: jj
      integer :: k, kk, kkk, kbdim, kedim, kzdim
      integer :: l, ll, lrzadim, lsig
      integer :: nampfmax_dim, ngenadim, nonch_dim, nrho_dim, nrho_dim1
      integer :: numrec1, numrecsave, nwkpack

!MPIINSERT_INCLUDE

!     This subroutine is only called from MPI rank=0.

!     This subroutine uses netCDF-3 subroutines.

!     This subroutine creates/writes a netCDF file
!       of standard data from cql3d.  (Additional data is
!       written for RF into a seperate file, if urfmod.ne."disabled").
!     netCDF file ==> "mnemonic".nc
!     netCDF file ID ==> ncid
!     time-step counter ==> numrec1
!
!     While moments of distn f and the electric field are saved after
!     each time-step, the distribution function itself is saved
!     only after the last time-step [netcdfshort='disabled'],
!     for n.eq.nsave() [if netcdfshort='lngshrtf'], or
!     after every time step [netcdfshort='longer_f']
!
!     Action of the subroutine is controlled by kopt and n:
!       kopt=0, n=0, initialize "mnemonic".nc and write data
!       kopt=1, n.gt.0, write data
!

! --- include file for netCDF declarations
! --- (obtained from NetCDF distribution)
      include 'netcdf.inc'

!     Will crop distribution functions to used size,
!     that is, take off the extra row and column of
!     points associated with f(0,0,k,l).


! --- some stuff for netCDF file ---
      integer ncid,vid,istatus
      integer xdim,ydim,rdim,r0dim,r00dim,kdim,zdim,tdim,gdim !-YuP: gdim added
      integer tsavedim
      integer twodim,fourdim,fivedim,thirteendim,chardim,char64dim
      integer nendim,nvdim,nmodsdim,ntavgdim,mrfnp1_dim,mrfnp3_dim
      integer nt_deltadim,nr_deltadim,nz_deltadim
      integer start2(2), start3(3), start4(4)

      integer nv_fusdim, nrho_fusdim !!!
      integer start_fusn(3), count_fusn(3), fusn_dims(3) !!!
      ! other fusion diagn: for sigftt and other diagn.  !!!
      integer start_fus(3),  count_fus(3),  fus_dims(3)  !!!
      !------

      integer nen_npadim,nv_npadim,npaproc_dim
      integer nena_dim, nva_dim, nrho_npadim
      integer npaproc_dims(2)
      integer count_npaproc(2),start_npaproc(2)
      integer npaenn_dims(2)
      integer count_npaenn(2),start_npaenn(2)
!
      integer count_hibr(3),start_hibr(3)
      integer count_sorpw(2), start_sorpw(2)
!
      integer dims(4),count(4),start(4)
      integer dimsg(5),countg(5),startg(5)  !To accnt for ngen.gt.1 cases
      integer dimssave(4),countsave(4),startsave(4)
      integer dimsgsave(5),countgsave(5),startgsave(5)  !For ngen.gt.1 cases
      integer dimsf(5),countf(5),startf(5) ! YuP[2018-09-28]
      integer count1(4),start1(4)
      integer count_xr(3),start_xr(3)
      integer count_npa(3),start_npa(3)
      integer count_powurf(2),start_powurf(2)
      integer count_rfpwr(3),start_rfpwr(3)
      integer count_powrf(3),start_powrf(3)
      integer count_powrft(2),start_powrft(2)
      integer count_powers(4),start_powers(4)
      integer r_dims(2),r_count(2)
      integer rk_dims(3),start_rk(3),count_rk(3)
      integer r0k_dims(3),start_r0k(3),count_r0k(3)
      integer r0_dims(2),r0_count(2)
      integer r00_dims(2),r00_count(2)
      integer y_dims(2),y_count(2)
      integer species_dims(3),species_count(3)
      integer kspeci_dims(3),kspeci_count(3)
      integer tau_dims(2),tau_count(2)
      integer z_dims(3),z_count(3),z_count1(3)
      integer delta_dims(3),delta_count(3),delta_start(3)
!deltarhop      Remove the !deltarhop to include deltarhop in .nc file
!deltarhop      integer deltap_dims(3),deltap_count(3),deltap_start(3)
      integer xray_dims(3)
      integer npa_dims(3)
      integer rfpwr_dims(3)
      integer dims_powurf(2)
      integer dims_powrft(2),dims_powrf(3)
      integer powers_dims(4)
      integer currv_dims(3)
      integer currv_dimsg(4)
!
      integer hibr_dims(3)  !  JK
      integer sorpw_dims(2)
!
      integer start_elecfldn(3),count_elecfldn(3),elecfldn_dims(3)

      character*100 ltitle

      character*8 ndelta_op

      real*8,  allocatable ::  wkpack(:) ! local working array for  pack21
      integer, allocatable :: iwkpack(:) ! local working array for ipack21

!      character(len=8), dimension(npaproc) :: npa_proc  !automatic var

      data start/1,1,1,1/, start1/1,1,1,1/
      data startsave/1,1,1,1/
      data startg/1,1,1,1,1/, start_rk/1,1,1/, start_r0k/1,1,1/
      data startgsave/1,1,1,1,1/
      data start_xr/1,1,1/
      data start_npa/1,1,1/
      data start_npaproc/1,1/
      data start_npaenn/1,1/
      !  For fusion diagnostics:
      data start_fusn/1,1,1/, start_fus/1,1,1/
      !--
      data start_powurf/1,1/
      data start_rfpwr/1,1,1/,start_powers/1,1,1,1/
      data start_powrf/1,1,1/
      data start_powrft/1,1/
      data delta_start/1,1,1/
      data start_hibr/1,1,1/
      data start_sorpw/1,1/
      data start2/1,1/
      data start3/1,1,1/
      data start4/1,1,1,1/
      
      integer j ! local

      ! another way to do this, when the code is a mess.
      character(len=8) :: cqlpmod
      integer :: lrindx(0:lrorsa)
      integer :: lrz
      integer :: lrzmax
      character(len=256) :: mnemonic
      cqlpmod = setup0%cqlpmod
      lrindx = setup0%lrindx
      lrz = setup0%lrz
      lrzmax =  setup0%lrzmax
      mnemonic = setup0%mnemonic


!deltarhop      data deltap_start/1,1,1/
!
!MPIINSERT_IF_RANK_NE_0_RETURN
! This subroutine is only called from MPI rank=0.

       WRITE(*,*) 'inside netcdfrw2.f...'
!       WRITE(*,*) 'frmodp=',frmodp
!       WRITE(*,*) 'hibrzp(i,1,1),hibrzp(i,2,1),hibrzp(i,3,1)'
!       do i=1,22
!       do i=1,nconteqn
!         WRITE(*,'(i4,2x,0p9f9.4)') i, hibrzp(i,1,1),
!     >        hibrzp(i,2,1),hibrzp(i,3,1),
!     >        hibrzp(i,1,2),hibrzp(i,2,2),hibrzp(i,3,2),
!     >        hibrzp(i,1,3),hibrzp(i,2,3),hibrzp(i,3,3)
!       enddo
!

!     Maximum iy as function of radius:
      iyy=0  !-YuP-101215: Don't use iy=; it's in common /params/
             ! Don't let overwrite the cqlinput value!
      do l=1,lrz
         iyy=max(iyy,iy_(l))
      enddo
      if(iyy.gt.iy) stop 'netcdfrw2: iy_(l) should not exceed iy'

      if (.NOT.ALLOCATED(wkpack)) then ! allocate working array for pack21
         nwkpack=max(iyjx2, iy*lrors, ntotala*(lrza+1), &
                     lza*lrzmax, iy*lz, nena*nva, &
                     iy*lrzmax, lrza*npaproca) +10 ! +10 just in case
         allocate(wkpack(nwkpack),STAT=istat)
         call bcast(wkpack,zero,SIZE(wkpack))
      endif

      if (.NOT.ALLOCATED(iwkpack)) then ! allocate working array for ipack21
         nwkpack=max(iyjx2, iy*lrors, ntotala*(lrza+1), &
                     lza*lrzmax, iy*lz, nena*nva, &
                     iy*lrzmax, lrza*npaproca) +10 ! +10 just in case
         allocate(iwkpack(nwkpack),STAT=istat)
         call ibcast(iwkpack,0,SIZE(iwkpack))
      endif

!     Ensure that tem1-2 has sufficient dim for (pack21) equivalences:
!      if (iy*lrors.gt.iyjx2)
!     +    stop '1: in netcdfrw2: Need to set jx>lrza'
!      if (ntotala*(lrza+1).gt.iyjx2)
!     +    stop '2: in netcdfrw2: Need ntotala*(lrza+1)<(iy+2)*(jx+2)'
!      if (lza*lrzmax.gt.iyjx2)
!     +    stop '3: in netcdfrw2: Need lza*lrza<(iy+2)*(jx+2)'
!      if (iy*lz.gt.iyjx2)
!     +    stop '4: in netcdfrw2: Need iy*lz<(iy+2)*(jx+2)'
!      if (nena*nva.gt.iyjx2)
!     +    stop '5: in netcdfrw2: Need nena*nva<(iy+2)*(jx+2)'
!
!      if (iy*lrzmax.gt.iyjx2)
!     +    stop '6: in netcdfrw2: Need iy*lrza<(iy+2)*(jx+2)'
!      if (lrza*npaproca.gt.iyjx2)
!     +    stop '7: in netcdfrw2: Need lrza*npaproca<(iy+2)*(jx+2)'


!BH110320C-----------------------------------------------------------------------
!BH110320C     Only set up for cqlpmod.ne."enabled",ngen=1, for the time being.
!BH110320C-----------------------------------------------------------------------

!BH110320      if (cqlpmod.eq."enabled" .or. ngen.ne.1) then
!BH110320         WRITE(*,*) 'WARNING: netcdfrw2 subroutine not implemented'
!BH110320         return
!BH110320      endif

!cdir$ master

      count(1)=iy
      count(2)=jx
      count(3)=lrz
      count(4)=1

      countg(1)=iy
      countg(2)=jx
      countg(3)=1  !radial index
      countg(4)=1  !species index
      countg(5)=1  !record number

      count1(1)=iy
      count1(2)=jx
      count1(3)=1  !radial index
      count1(4)=1  !record number

      count_xr(1)=nen
      count_xr(2)=nv
      count_xr(3)=1

      count_npa(1)=nen_npa
      count_npa(2)=nv_npa
      count_npa(3)=1

      count_npaenn(1)=lrzmax
      count_npaenn(2)=npaproc

      count_npaproc(1)=8
      count_npaproc(2)=npaproc

      count_fus(1)=lrzmax
      count_fus(2)=4
      count_fus(3)=1

      count_powurf(1)=mrfn+1
      count_powurf(2)=1

      count_rfpwr(1)=lrz
      count_rfpwr(2)=mrfn+3
      count_rfpwr(3)=1

      count_powrf(1)=lrz
      count_powrf(2)=nmodsa
      count_powrf(3)=1

      count_powrft(1)=lrz
      count_powrft(2)=1

      count_powers(1)=lrz
      count_powers(2)=13
      count_powers(3)=ngen
      count_powers(4)=1

      r_count(1)=lrz
      r_count(2)=1

      count_rk(1)=lrz
      count_rk(2)=ngen
      count_rk(3)=1

      count_r0k(1)=lrzmax
      count_r0k(2)=ngen
      count_r0k(3)=1

      r0_count(1)=lrzmax
      r0_count(2)=1

      r00_count(1)=lrzmax+1
      r00_count(2)=1

      y_count(1)=iy
      y_count(2)=lrz

      species_count(1)=ntotal
      species_count(2)=lrzmax
      species_count(3)=1

      tau_count(1)=iy
      tau_count(2)=lrzmax

      z_count(1)=iy
      z_count(2)=lz
      z_count(3)=lrzmax

      z_count1(1)=iy
      z_count1(2)=lz
      z_count1(3)=1

      delta_count(1)=nr_delta
      delta_count(2)=nz_delta
      delta_count(3)=nt_delta

!deltarhop      deltap_count(1)=nt_delta
!deltarhop      deltap_count(2)=lz
!deltarhop      deltap_count(3)=lrzmax

      kspeci_count(1)=8
      kspeci_count(2)=2
      kspeci_count(3)=ntotal
!
      count_hibr(1)=kz
      count_hibr(2)=ke
      count_hibr(3)=kb
      count_sorpw(1)=ngena
      count_sorpw(2)=lrza
!
!.......................................................................
!cl    1. Initialize part, creating new netcdf file
!

! --- begin if ---
      WRITE(*,*)'netcdfrw2 initialize, kopt=',kopt,'n=',n
      if ((kopt.eq.0) .and. (n.eq.0)) then !endif at l 3282

!-----------------------------------------------------------------------
!
!cl     1.1 create netCDF file and define dimensions,variables
!          and attributes
!

!.......................................................................
!cl    1.1.1 create netCDF filename
!     CLOBber old file, if it exists.
!     istatus is 0, if no errors.
      write(t_,1000) mnemonic(1:length_char(mnemonic))
 1000 format(a,".nc")
      istatus = NF_CREATE(t_, NF_CLOBBER, ncid) !-YuP: NetCDF-f77
      call check_err(istatus)

!.......................................................................
!cl    1.1.2 define dimensions
!     p. 43 of netcdf-3 manual

      istatus= NF_DEF_DIM(ncid, 'xdim',     jx,       xdim)  !-YuP: NetCDF-f77
      istatus= NF_DEF_DIM(ncid, 'ydim',     iy,       ydim)
      istatus= NF_DEF_DIM(ncid, 'rdim',     lrz,      rdim)
      istatus= NF_DEF_DIM(ncid, 'r0dim',    lrzmax,   r0dim)
      istatus= NF_DEF_DIM(ncid, 'kzdim',     kz,      kzdim) ! JK: freya
      istatus= NF_DEF_DIM(ncid, 'kedim',     ke,      kedim)
      istatus= NF_DEF_DIM(ncid, 'kbdim',     kb,      kbdim)
      istatus= NF_DEF_DIM(ncid, 'ngenadim',  ngena,   ngenadim)
      istatus= NF_DEF_DIM(ncid, 'lrzadim',   lrza,    lrzadim)
!      WRITE(*,*) 'kzdim,kedim,kbdim = ',kzdim,kedim,kbdim
!      stop
      istatus= NF_DEF_DIM(ncid, 'twodim',       2,    twodim)
      istatus= NF_DEF_DIM(ncid, 'fourdim',      4,    fourdim)
      istatus= NF_DEF_DIM(ncid, 'fivedim',      5,    fivedim)
      istatus= NF_DEF_DIM(ncid, 'thirteendim', 13,    thirteendim)
      istatus= NF_DEF_DIM(ncid, 'r00dim',   lrzmax+1, r00dim)
      istatus= NF_DEF_DIM(ncid, 'zdim',     lz,       zdim)
      istatus= NF_DEF_DIM(ncid, 'nendim',   nen,      nendim)
      istatus= NF_DEF_DIM(ncid, 'nvdim',    nv,       nvdim)
      istatus= NF_DEF_DIM(ncid, 'nen_npadim',   nen_npa,  nen_npadim)
      istatus= NF_DEF_DIM(ncid, 'nv_npadim',    nv_npa,    nv_npadim)
      istatus= NF_DEF_DIM(ncid, 'npaproc_dim',  npaproc, npaproc_dim)
      istatus= NF_DEF_DIM(ncid, 'mrfnp1_dim', mrfn+1,   mrfnp1_dim)
      istatus= NF_DEF_DIM(ncid, 'mrfnp3_dim', mrfn+3,   mrfnp3_dim)
      istatus= NF_DEF_DIM(ncid, 'nmodsdim', nmodsa,    nmodsdim)
      istatus= NF_DEF_DIM(ncid, 'ntavgdim', ntavga,    ntavgdim)

      istatus= NF_DEF_DIM(ncid, 'gen_species_dim', ngen,   gdim)
      istatus= NF_DEF_DIM(ncid, 'species_dim',     ntotal, kdim)

      istatus= NF_DEF_DIM(ncid,'nt_deltadim', nt_delta, nt_deltadim)
      istatus= NF_DEF_DIM(ncid,'nr_deltadim', nr_delta, nr_deltadim)
      istatus= NF_DEF_DIM(ncid,'nz_deltadim', nz_delta, nz_deltadim)

      istatus= NF_DEF_DIM(ncid, 'chardim',      8,    chardim)
      istatus= NF_DEF_DIM(ncid, 'char64dim',   64,    char64dim)

      if(ampfmod.eq.'enabled')then
      start_elecfldn(1)=1 !0
      start_elecfldn(2)=1 !0
      start_elecfldn(3)=1
      count_elecfldn(1)=lrz+2 !Starts at l=0. Put boundary value as last
                              !radial entry.
      count_elecfldn(2)=nstop+1
      count_elecfldn(3)=nampfmax+1
      print*,'count_elecfldn(1:3)=',count_elecfldn
      istatus= NF_DEF_DIM(ncid,'nampfmax_dim',nampfmax+1,nampfmax_dim)
      istatus= NF_DEF_DIM(ncid,'nonch_dim',   nstop+1,   nonch_dim)
      istatus= NF_DEF_DIM(ncid,'nrho_dim',    lrz+1,     nrho_dim)
      istatus= NF_DEF_DIM(ncid,'nrho_dim1',    lrz+2,     nrho_dim1)
      endif

!     unlimited dimension for time, dimension name= 'tdim'
      !istatus= NF_DEF_DIM(ncid, 'tdim',NF_UNLIMITED,tdim)
      !YuP: why NF_UNLIMITED needed?
      istatus= NF_DEF_DIM(ncid, 'tdim',nstop+1,tdim) !YuP[2018-09-28]
      !changed to nstop+1, to accomodate dimsf() logic, see below.

      !YuP[2018-09-28] Added 'tsavedim' to save timet
      ! at selected t steps (corr. to nsave(i)=n steps)
      ! Can have only one NF_UNLIMITED ?
      !istatus=NF_DEF_DIM(ncid,'tsavedim',NF_UNLIMITED,tsavedim) ! failed
      ! nsavet is counted in tdchief
      istatus= NF_DEF_DIM(ncid,'tsavedim',nsavet+1,tsavedim)

!     define vector of dimensions, unlimited last
      dims(1)=ydim
      dims(2)=xdim
      dims(3)=rdim
      dims(4)=tdim

      dimsg(1)=ydim
      dimsg(2)=xdim
      dimsg(3)=rdim
      dimsg(4)=gdim
      dimsg(5)=tdim ! can be re-defined, see below

      !dimsf() are used for saving of f() only.
      dimsf(1)=ydim
      dimsf(2)=xdim
      dimsf(3)=rdim
      ! Now define the remaining dimsf(4:5)
      if (netcdfshort.eq.'enabled') then
         !Do nothing: no storage defined.
      elseif ( (netcdfshort.eq.'longer_f').or. &
               (netcdfshort.eq.'lngshrtf')     ) then ! define storage
         if (ngen.eq.1) then    !maintaining backwards compatability
            !to be set below: vid=ncvdef2(ncid,'f',NCDOUBLE,4,dimsf,istatus)
            dimsf(4)=tdim
            if(netcdfshort.eq.'lngshrtf') dimsf(4)=tsavedim
         else  !ngen.ge.2
            !to be set below: vid=ncvdef2(ncid,'f',NCDOUBLE,5,dimsg,istatus)
            !Additional dim included for ngen.gt.1 cases
            dimsg(4)=gdim ! to contain 'k' index
            dimsg(5)=tdim !
            if(netcdfshort.eq.'lngshrtf') dimsg(5)=tsavedim
         endif  !on ngen
      else     !disabled, Standard o/p: f at last time step
         if (ngen.eq.1) then    !maintaining backwards compatability
            !to be set below: vid=ncvdef2(ncid,'f',NCDOUBLE,3,dimsf,istatus)
         else  !ngen.ge.2
            !to be set below: vid=ncvdef2(ncid,'f',NCDOUBLE,4,dimsg,istatus)
            !Additional dim included for ngen.gt.1 cases
            dimsg(4)=gdim !! to contain 'k' index
         endif  !on ngen
      endif  !on netcdfshort : storage defined


      r_dims(1)=rdim
      r_dims(2)=tdim

      rk_dims(1)=rdim
      rk_dims(2)=gdim
      rk_dims(3)=tdim

      r0_dims(1)=r0dim
      r0_dims(2)=tdim

      r0k_dims(1)=r0dim
      r0k_dims(2)=gdim
      r0k_dims(3)=tdim

      r00_dims(1)=r00dim
      r00_dims(2)=tdim

      y_dims(1)=ydim
      y_dims(2)=rdim

      species_dims(1)=kdim
      species_dims(2)=r0dim
      species_dims(3)=tdim

      tau_dims(1)=ydim
      tau_dims(2)=r0dim

      z_dims(1)=ydim
      z_dims(2)=zdim
      z_dims(3)=r0dim

      if (softxry.ne."disabled" .and. softxry.ne."ncdf_all") then
         xray_dims(1)=nendim
         xray_dims(2)=nvdim
         xray_dims(3)=twodim    !First and last sets of spectra
      elseif ( softxry.eq."ncdf_all") then
         xray_dims(1)=nendim
         xray_dims(2)=nvdim
         xray_dims(3)=tdim      !All sets of spectra
      endif

      if (npa_diag.ne."disabled" .and. npa_diag.ne."ncdf_all") then
         npa_dims(1)=nen_npadim
         npa_dims(2)=nv_npadim
         npa_dims(3)=twodim    !First and last sets of spectra
      elseif ( npa_diag.eq."ncdf_all") then
         npa_dims(1)=nen_npadim
         npa_dims(2)=nv_npadim
         npa_dims(3)=tdim      !All sets of spectra
      endif

      if (npa_diag.ne."disabled") then
         npaproc_dims(1)=chardim
         npaproc_dims(2)=npaproc_dim
         npaenn_dims(1)=r0dim
         npaenn_dims(2)=npaproc_dim
      endif
!
      hibr_dims(2)=kedim
      hibr_dims(3)=kbdim
      sorpw_dims(1)=ngenadim
      sorpw_dims(2)=lrzadim
!
      fus_dims(1)=r0dim       !Radial profiles at last step
      fus_dims(2)=fourdim
      fus_dims(3)=tdim

      dims_powurf(1)=mrfnp1_dim
      dims_powurf(2)=tdim

      rfpwr_dims(1)=rdim
      rfpwr_dims(2)=mrfnp3_dim
      rfpwr_dims(3)=tdim

      dims_powrf(1)=rdim
      dims_powrf(2)=nmodsdim
      dims_powrf(3)=tdim

      dims_powrft(1)=rdim
      dims_powrft(2)=tdim

      powers_dims(1)=rdim
      powers_dims(2)=thirteendim
      powers_dims(3)=gdim
      powers_dims(4)=tdim

      currv_dims(1)=xdim
      currv_dims(2)=rdim
      currv_dims(3)=tdim

      currv_dimsg(1)=xdim
      currv_dimsg(2)=rdim
      currv_dimsg(3)=gdim
      currv_dimsg(4)=tdim

      delta_dims(1)=nr_deltadim
      delta_dims(2)=nz_deltadim
      delta_dims(3)=nt_deltadim

!deltarhop      deltap_dims(1)=nt_deltadim
!deltarhop      deltap_dims(2)=zdim
!deltarhop      deltap_dims(3)=r0dim

      kspeci_dims(1)=chardim
      kspeci_dims(2)=twodim
      kspeci_dims(3)=kdim

      if(ampfmod.eq.'enabled')then
      elecfldn_dims(1)= nrho_dim1 ! correspond to lrz+2 (0:lrz+1)
      elecfldn_dims(2)= nonch_dim ! corr. to nstop+1   (0:nstop)
      elecfldn_dims(3)= nampfmax_dim                 ! (0:nampfmax)
      endif

!.......................................................................
!cl    1.1.3 define variables

!     Note, the variable IDs (denoted below as vid) are
!     not saved here in this subroutine; rather, the IDs
!     are retrieved from the netCDF data file, as needed,
!     by calling the netCDF routine ncvid.

!     netcdf variable_define:
!     nf_def_var(netcdf_id,variable_name,variable_type,
!                number_of_dimensions,
!                vector_for_length_per_dimension,varid)
!     Note: Unlimited dimension must be last
!     Refer to p. 50, netcdf-3 manual

!     NF_DOUBLE for REAL*8:


!--------------------------
!     Time-independent data
!--------------------------

      ltitle='Main netcdf output from CQL3D version: '//version
      if( length_char(ltitle).gt.100 ) stop 'Adjust ltitle in netcdfrw2'
      call ncaptc2(ncid,NCGLOBAL,'title',NCCHAR,length_char(ltitle), &
           ltitle,istatus)

      vid=ncvdef0(ncid,'version',NCCHAR,1,char64dim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,20, &
                'CQL3D version number',istatus)

      vid=ncvdef0(ncid,'mnemonic',NCCHAR,1,char64dim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,23, &
                'Mnemonic run identifier',istatus)

      vid=ncvdef0(ncid,'ampfmod',NCCHAR,1,chardim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,27, &
                'Ampere-Faraday module switch',istatus)

      vid=ncvdef0(ncid,'urfmod',NCCHAR,1,chardim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,17, &
                'URF module switch',istatus)

      vid=ncvdef0(ncid,'rdcmod',NCCHAR,1,chardim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,45, &
                'RF Heating and CD from input diffusion coeffs',istatus)

      vid=ncvdef0(ncid,'frmod',NCCHAR,1,chardim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,17, &
                'NBI module switch',istatus)

      vid=ncvdef0(ncid,'beamplse',NCCHAR,1,chardim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,17, &
                'Pulsed NBI switch',istatus)

      vid=ncvdef0(ncid,'transp',NCCHAR,1,chardim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30, &
                'Radial transport module switch',istatus)

      vid=ncvdef0(ncid,'tavg',NCCHAR,1,chardim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,37, &
                'Indicates calc of time-avg distn favg',istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,40, &
                'in which case favg is o/p in place of f,',istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,35, &
                'except in netcdfshort=longer_f case',istatus)

!.......................................................................

      vid=ncvdef0(ncid,'f4d_out',NCCHAR,1,chardim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,44, &
                'Indicates output for 4D distn/ separate file',istatus)

!BH181025:  Not used
!      vid=ncvdef0(ncid,'f3d_out',NCCHAR,1,chardim,istatus)
!      call ncaptc2(ncid,vid,'long_name',NCCHAR,44, &
!                'Indicates output for 3D distn/ separate file',istatus)

      vid=ncvdef0(ncid,'netcdfshort',NCCHAR,1,chardim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,46, &
              'Indicates level of output of data to .nc files',istatus)

      vid=ncvdef0(ncid,'eqdskin',NCCHAR,1,char64dim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,39, &
                'Name of input eqdsk, for eqsource=eqdsk',istatus)

      vid=ncvdef0(ncid,'ngen',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,25, &
                 'Number of general species',istatus)
      call check_err(istatus)

      vid=ncvdef0(ncid,'ntotal',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,26, &
                 'Number of species: gen+max',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'kspeci',NCCHAR,3,kspeci_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,44, &
                'Name of species and spec as general or maxwl',istatus)

      vid=ncvdef0(ncid,'bnumb',NCDOUBLE,(1),kdim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,29, &
                 'atomic number or each species',istatus)
      call check_err(istatus)

      vid=ncvdef0(ncid,'fmass',NCDOUBLE,(1),kdim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,20, &
                 'mass of each species',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5, &
                           'grams',istatus)
      call check_err(istatus)

      vid=ncvdef0(ncid,'lrzmax',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,38, &
                'Number of radial bins(=r0dim, .ge.lrz)',istatus)
      call check_err(istatus)
      vid=ncvdef0(ncid,'radcoord',NCCHAR,1,chardim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,45, &
                'Radial coordinate is proportional to radcoord',istatus)

      vid=ncvdef0(ncid,'rya',NCDOUBLE,(1),r0dim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,37, &
                 'Normalized radial mesh at bin centers',istatus)

      vid=ncvdef0(ncid,'Rp',NCDOUBLE,1,r0dim, istatus)   ! rpcon(lrz) array
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,41, &
                 'Major radius of bin centers at outerboard',istatus)
      call check_err(istatus)

      vid=ncvdef0(ncid,'Rm',NCDOUBLE,1,r0dim,istatus)   ! rmcon(lrz) array
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,41, &
                 'Major radius of bin centers at innerboard',istatus)
      call check_err(istatus)

      vid=ncvdef0(ncid,'rhomax',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,31, &
                 'Generalized plasma minor radius',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                 'cms',istatus)

      vid=ncvdef0(ncid,'radmaj',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,20, &
                 'Nominal major radius',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                 'cms',istatus)

      vid=ncvdef0(ncid,'rpmconz',NCDOUBLE,(1),r00dim,istatus)
      call ncaptc2(ncid,vid,'long_name1',NCCHAR,43, &
                 'Major radius of bin boundaries on the outer',istatus)
      call ncaptc2(ncid,vid,'long_name2',NCCHAR,36, &
                 'equatorial plane, starting at radmaj',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                 'cms',istatus)

      vid=ncvdef0(ncid,'btor',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,29, &
                 'Nominal tor mag fld at radmaj',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5, &
                 'Gauss',istatus)

      vid=ncvdef0(ncid,'toteqd',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30, &
                 'Tor equilibrium plasma current',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,4, &
                 'Amps',istatus)

      vid=ncvdef0(ncid,'rgeomp',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,29, &
                 '0.5*(max-min) of major radius',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                 'cms',istatus)

      vid=ncvdef0(ncid,'r0geomp',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,29, &
                 '0.5*(max+min) of major radius',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                 'cms',istatus)

      vid=ncvdef0(ncid,'rmag',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,26, &
                 'Magnetic axis major radius',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                 'cms',istatus)

      vid=ncvdef0(ncid,'zmag',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,31, &
                 'Magnetic axis vertical position',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                 'cms',istatus)

      vid=ncvdef0(ncid,'eqsym',NCCHAR,1,chardim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,28, &
                'Indicator for symmetrization',istatus)

!.......................................................................
      vid=ncvdef0(ncid,'zshift',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,39, &
                 'Vertical shift of equilibrium per eqsym',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                 'cms',istatus)

      vid=ncvdef0(ncid,'eps0',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,32, &
                 'Inv aspect ratio = rgeomp/r0geomp',istatus)

      vid=ncvdef0(ncid,'elong',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,26, &
                 'Elongation = zgeomp/rgeomp',istatus)

      vid=ncvdef0(ncid,'zgeomp',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,26, &
                 'Approx half-height to LCFS',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                 'cms',istatus)

      vid=ncvdef0(ncid,'rgeom1',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,33, &
                 'Approx inner major radius to LCFS',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                 'cms',istatus)

      vid=ncvdef0(ncid,'rgeom2',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,33, &
                 'Approx outer major radius to LCFS',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                 'cms',istatus)

      vid=ncvdef0(ncid,'area',NCDOUBLE,(1),r0dim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30, &
                 'Cumulative area to bin centers',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5, &
                 'cms^2',istatus)

      vid=ncvdef0(ncid,'darea',NCDOUBLE,(1),r0dim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,35, &
                 'Incremental tor area of radial bins',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5, &
                 'cms^2',istatus)

      vid=ncvdef0(ncid,'vol',NCDOUBLE,(1),r0dim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,32, &
                 'Cumulative volume to bin centers',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5, &
                 'cms^3',istatus)

      vid=ncvdef0(ncid,'dvol',NCDOUBLE,(1),r0dim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,33, &
                 'Incremental volume to bin centers',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5, &
                 'cms^3',istatus)

      vid=ncvdef0(ncid,'equilpsi',NCDOUBLE,(1),r0dim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,37, &
                 'Poloidal flux function at bin centers',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                 'cgs',istatus)

      vid=ncvdef0(ncid,'psivalm',NCDOUBLE,(1),r0dim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,46, &
                 'Pol flux fctn at radial outer edge of each bin', &
                 istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                 'cgs',istatus)

      vid=ncvdef0(ncid,'dpsi',NCDOUBLE,(1),r0dim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,44, &
                 'Incremental pol flux function at bin centers',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                 'cgs',istatus)

      vid=ncvdef0(ncid,'psimag',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,34, &
                 'Pol flux function at magnetic axis',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                 'cgs',istatus)

      vid=ncvdef0(ncid,'psilim',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,36, &
                 'Pol flux function at plasma boundary',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                 'cgs',istatus)

      vid=ncvdef0(ncid,'h_r',NCDOUBLE,(1),r0dim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,28, &
                 'H*rho ONETWO geometry factor',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                 'cms',istatus)

      vid=ncvdef0(ncid,'qsafety',NCDOUBLE,(1),r0dim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,28, &
                 'Safety factor at bin centers',istatus)

!     Output current calculated from eqdsk:
      if (eqmod.eq."enabled") then
      vid=ncvdef0(ncid,'curreq',NCDOUBLE,(1),r0dim,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,35, &
              'Toroidal current density from EQDSK',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,9, &
              'Amps/cm**2',istatus)
      endif

      vid=ncvdef0(ncid,'lrz',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,41, &
                  'Number of FPd radial surface bins (=rdim)',istatus)
      call check_err(istatus)

      vid=ncvdef0(ncid,'lrindx',NCLONG,1,rdim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30, &
                 'Radial indices of FPd surfaces',istatus)

      vid=ncvdef0(ncid,'jx',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,35, &
                  'momentum-per-mass dimension (=xdim)',istatus)
      call check_err(istatus)

      vid=ncvdef0(ncid,'x',NCDOUBLE,(1),xdim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,28, &
                 'normalized momentum-per-mass',istatus)

      vid=ncvdef0(ncid,'enerkev',NCDOUBLE,(1),xdim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,31, &
                 'energy=restmkev*(gamma(x(:))-1)',istatus)

      vid=ncvdef0(ncid,'uoc',NCDOUBLE,(1),xdim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,34, &
                 'uoc=x(:)/cnorm=mom-per-mass/clight',istatus)

      vid=ncvdef0(ncid,'dx',NCDOUBLE,(1),xdim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,28, &
                 'dx centered on x-mesh points',istatus)

      vid=ncvdef0(ncid,'cint2',NCDOUBLE,(1),xdim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,7, &
                 'x**2*dx',istatus)

      vid=ncvdef0(ncid,'vnorm',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,33, &
                 'velocity (momentum-per-mass) norm',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7, &
                           'cms/sec',istatus)

      vid=ncvdef0(ncid,'enorm',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,20, &
                           'Energy normalization',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                           'keV',istatus)

      vid=ncvdef0(ncid,'iy',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,33, &
                  'max Pitch angle dimension (=ydim)',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'y',NCDOUBLE,2,y_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,11, &
                 'pitch angle',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7, &
                 'radians',istatus)

      vid=ncvdef2(ncid,'dy',NCDOUBLE,2,y_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,28, &
                 'dy centered on y-mesh points',istatus)

      vid=ncvdef2(ncid,'cynt2',NCDOUBLE,2,y_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,14, &
                 '2*pi*sin(y)*dy',istatus)

      vid=ncvdef0(ncid,'iy_',NCLONG,1,rdim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,41, &
                  'Pitch angle dimension at each r (le ydim)',istatus)
      call check_err(istatus)

      vid=ncvdef0(ncid,'itl',NCLONG,1,rdim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,32, &
                 'lower trapped-passing bndy index',istatus)

      vid=ncvdef0(ncid,'itu',NCLONG,1,rdim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,32, &
                 'upper trapped-passing bndy index',istatus)

      vid=ncvdef0(ncid,'lz',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,27, &
                  'dimension of z-grid along B',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'z',NCDOUBLE,2,z_dims(2),istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,16, &
                 'Distance along B',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                           'cms',istatus)

      vid=ncvdef2(ncid,'dz',NCDOUBLE,2,z_dims(2),istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,23, &
                 'dz centered on z-points',istatus)

      vid=ncvdef2(ncid,'solrz',NCDOUBLE,2,z_dims(2),istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,21, &
                 'Major radius of z pts',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                           'cms',istatus)

      vid=ncvdef2(ncid,'solzz',NCDOUBLE,2,z_dims(2),istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,15, &
                 'Height of z pts',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                           'cms',istatus)

      vid=ncvdef2(ncid,'pol',NCDOUBLE,2,z_dims(2),istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,39, &
                 'Poloidal angle, measured about mag axis',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7, &
                           'radians',istatus)

      vid=ncvdef2(ncid,'bbpsi',NCDOUBLE,2,z_dims(2),istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,12, &
                 'B(z)/B(z=0)',istatus)

      vid=ncvdef2(ncid,'imax',NCLONG,2,z_dims(2),istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,27, &
                 'Max i s.t. part passes z(l)',istatus)

      vid=ncvdef2(ncid,'lmax',NCLONG,2,tau_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,24, &
                 'Max l s.t. i passes z(l)',istatus)

      vid=ncvdef2(ncid,'zboun',NCDOUBLE,2,tau_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,20, &
                 'Bounce point z-value',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                           'cms',istatus)

      vid=ncvdef0(ncid,'zmaxpsi',NCDOUBLE,(1),r0dim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,15, &
                 'Integral dz/bbpsi',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                           'cms',istatus)

      vid=ncvdef2(ncid,'tau',NCDOUBLE,2,tau_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,23, &
                 'tau_bounce * abs(speed)',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                           'cms',istatus)

      vid=ncvdef2(ncid,'dtau',NCDOUBLE,3,z_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,24, &
                'abs(speed)*dtau in dz(l)',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                           'cms',istatus)

      vid=ncvdef0(ncid,'beampon',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,28, &
                'On time, per cylce, of beam',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,4, &
                           'secs',istatus)

      vid=ncvdef0(ncid,'beampoff',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,29, &
                'Off time, per cylce, of beam',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,4, &
                           'secs',istatus)

      vid=ncvdef0(ncid,'tavg1',NCDOUBLE,(1),ntavgdim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,31, &
                'tavg case:  interval start time',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,4, &
                           'secs',istatus)

      vid=ncvdef0(ncid,'tavg2',NCDOUBLE,(1),ntavgdim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30, &
                'tavg case:  interval stop time',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,4, &
                           'secs',istatus)

      vid=ncvdef0(ncid,'ndeltarho',NCCHAR,1,chardim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,41, &
                'Indicates first order-orbit shift options',istatus)

!.......................................................................

      ndelta_op="enabled"
      if ((ndeltarho.ne.'disabled'.or.lossmode(1).eq.'simplban').and. &
           ndelta_op.eq."enabled") then

         vid=ncvdef2(ncid,'deltarho',NCDOUBLE,3,z_dims,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,43, &
              'Orbit shift FROM BA FS vs theta0,z,rya,/|v|',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,26, &
              "Norm'd to rhomax*|x*vnorm|",istatus)

!deltarhop         vid=ncvdef(ncid,'deltarhop',NCDOUBLE,3,deltap_dims,istatus)
!deltarhop         call ncaptc2(ncid,vid,'long_name',NCCHAR,48,
!deltarhop     +       'Orbit shift FROM BA FS vs local theta,z,rya,/|v|',istatus)
!deltarhop         call ncaptc2(ncid,vid,'units',NCCHAR,26,
!deltarhop     +        "Norm'd to rhomax*|x*vnorm|",istatus)

         vid=ncvdef0(ncid,'r_delta',NCDOUBLE,(1),nr_deltadim,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,25, &
              'deltarz major radius mesh',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,3, &
              'cms',istatus)

         vid=ncvdef0(ncid,'z_delta',NCDOUBLE,(1),nz_deltadim,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,23, &
              'deltarz vertical Z-mesh',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,3, &
              'cms',istatus)

         vid=ncvdef0(ncid,'t_delta',NCDOUBLE,(1),nt_deltadim,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,30, &
              'deltarz local pitch angle mesh',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,7, &
              'radians',istatus)

         vid=ncvdef2(ncid,'deltarz',NCDOUBLE,3,delta_dims,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,46, &
              'Orbit shift TO BA FS vs local theta,R,Z,/|v|',istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,49, &
              'Thus, co-current ions at any R,Z have neg deltarz', &
              istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,54, &
              "Norm'd to max of coord specified by radcoord*|x*vnorm|", &
              istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,55, &
              'I.e., normed to 1.*velocity (cms/sec) for most radcoord', &
              istatus)

         vid=ncvdef2(ncid,'delta_bdb0',NCDOUBLE,2,delta_dims,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,54, &
              'Ratio of mag fld B at deltarz R,Z over B0 at midplane', &
              istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,8, &
              "unitless",istatus)

      endif  ! On deltarho
      if(ampfmod.eq.'enabled')then
         vid=ncvdef2(ncid,'elecfldn',NCDOUBLE,3,elecfldn_dims,istatus)
!         WRITE(*,*)'elecfldn_dims,ncid,vid,istatus=', &
!                    elecfldn_dims,ncid,vid,istatus
         call ncaptc2(ncid,vid,'long_name1',NCCHAR,53, &
         'Tor. Electric Field as func of iteration,time,radius,', &
         istatus)
         call ncaptc2(ncid,vid,'long_name3',NCCHAR,48, &
         'Bin centered, except last value the bndry value.',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,12, &
                               'statVolts/cm',istatus)
         ! elecfldn,delecfld0 are in cgs !
      endif


      vid=ncvdef0(ncid,'bthr',NCDOUBLE,(1),r0dim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,35, &
                 'Equil Pol B field at theta_pol=pi/2',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5, &
                           'Gauss',istatus)

      vid=ncvdef0(ncid,'btoru',NCDOUBLE,(1),r0dim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,35, &
                 'Equil Tor B field at theta_pol=pi/2',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5, &
                           'Gauss',istatus)

      vid=ncvdef0(ncid,'btor0',NCDOUBLE,(1),r0dim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,48, &
             'Tor mag fld strength at min |B| on flux surfaces',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5, &
                           'Gauss',istatus)

      vid=ncvdef0(ncid,'bmidplne',NCDOUBLE,(1),r0dim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,38, &
                 'Min mag fld |B| on a rad flux surfaces',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5, &
                           'Gauss',istatus)

      vid=ncvdef0(ncid,'efflag',NCCHAR,1,chardim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,41, &
                'Indicates elecfld is toroidal or parallel',istatus)

!     X-ray data:

      vid=ncvdef0(ncid,'softxry',NCCHAR,1,chardim,istatus)
      call ncaptc2(ncid,vid,'long_name0',NCCHAR,55, &
           'X-ray diagnostic is disabled,enabled, ncdf_all, or e-ion' , &
           istatus)
      call ncaptc2(ncid,vid,'long_name1',NCCHAR,36, &
            'Data for 1st and last step is output' ,istatus)
      call ncaptc2(ncid,vid,'long_name2',NCCHAR,35, &
            'or, data for all steps for ncdf_all' ,istatus)
      call check_err(istatus)

      if (softxry .ne. "disabled") then

         if (x_sxr(1).ne.zero  .or. z_sxr(1).ne.zero) then

      vid=ncvdef0(ncid,'x_sxr',NCDOUBLE,(1),nvdim,istatus)
            call ncaptc2(ncid,vid,'long_name',NCCHAR,42, &
                 'X-ray detector major radius (tor angle=0.)',istatus)
            call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                 'cms',istatus)

      vid=ncvdef0(ncid,'z_sxr',NCDOUBLE,(1),nvdim,istatus)
            call ncaptc2(ncid,vid,'long_name',NCCHAR,21, &
                 'X-ray detector height',istatus)
            call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                 'cms',istatus)

         else

      vid=ncvdef0(ncid,'rd',NCDOUBLE,(1),nvdim,istatus)
            call ncaptc2(ncid,vid,'long_name',NCCHAR,27, &
                 'X-ray detector minor radius ',istatus)
            call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                 'cms',istatus)

      vid=ncvdef0(ncid,'thetd',NCDOUBLE,(1),nvdim,istatus)
            call ncaptc2(ncid,vid,'long_name',NCCHAR,29, &
                 'X-ray detector poloidal angle',istatus)
            call ncaptc2(ncid,vid,'units',NCCHAR,7, &
                 'degrees',istatus)

         endif

      vid=ncvdef0(ncid,'nv',NCLONG,0,0,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,29, &
              'Number of X-ray viewing cords',istatus)
         call check_err(istatus)

      vid=ncvdef0(ncid,'nen',NCLONG,0,0,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,40, &
              'Number of equispaced energies in spectra',istatus)
         call check_err(istatus)

      vid=ncvdef0(ncid,'msxr',NCLONG,0,0,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,47, &
              'Highest order of Legendre poly in f expr for XR',istatus)
         call check_err(istatus)

      vid=ncvdef0(ncid,'enmin',NCDOUBLE,0,0,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,31, &
              'X-ray minimun energy in spectra',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,3, &
              'keV',istatus)

      vid=ncvdef0(ncid,'enmax',NCDOUBLE,0,0,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,31, &
              'X-ray maximum energy in spectra',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                 'keV',istatus)

      vid=ncvdef0(ncid,'en_',NCDOUBLE,(1),nendim,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,13, &
              'Photon Energy',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,3, &
              'keV',istatus)

      vid=ncvdef2(ncid,'eflux',NCDOUBLE,3,xray_dims,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,17, &
              'X-ray energy flux',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,22, &
              'ergs/cm**2/sec/ster/eV',istatus)

      vid=ncvdef2(ncid,'efluxt',NCDOUBLE,2,xray_dims(2),istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,28, &
              'Integrated X-ray energy flux',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,19, &
              'ergs/cm**2/sec/ster',istatus)

      endif ! On softxry .ne. "disabled"


!     NPA data:

      vid=ncvdef0(ncid,'npa_diag',NCCHAR,1,chardim,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name0',NCCHAR,47, &
            'NPA diagnostic is disabled,enabled, or ncdf_all' ,istatus)
      call ncaptc2(ncid,vid,'long_name1',NCCHAR,45, &
            'enabled: Data for 1st and last step is output' ,istatus)
      call ncaptc2(ncid,vid,'long_name2',NCCHAR,45, &
            'or, ncdf_all: data for all steps for ncdf_all' ,istatus)

      if (npa_diag .ne. "disabled") then

         if (x_npa(1).ne.zero  .or. z_npa(1).ne.zero) then

            vid=ncvdef0(ncid,'x_npa',NCDOUBLE,(1),nv_npadim,istatus)
            call check_err(istatus)
            call ncaptc2(ncid,vid,'long_name',NCCHAR,40, &
                 'NPA detector major radius (tor angle=0.)',istatus)
            call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                 'cms',istatus)

            vid=ncvdef0(ncid,'z_npa',NCDOUBLE,(1),nv_npadim,istatus)
            call check_err(istatus)
            call ncaptc2(ncid,vid,'long_name',NCCHAR,19, &
                 'NPA detector height',istatus)
            call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                 'cms',istatus)

         else ! rd_npa/thetd_npa input

            vid=ncvdef0(ncid,'rd_npa',NCDOUBLE,(1),nv_npadim,istatus)
            call check_err(istatus)
            call ncaptc2(ncid,vid,'long_name',NCCHAR,25, &
                 'NPA detector minor radius ',istatus)
            call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                 'cms',istatus)

            vid=ncvdef0(ncid,'thetd_npa',NCDOUBLE,(1),nv_npadim,istatus)
            call ncaptc2(ncid,vid,'long_name',NCCHAR,27, &
                 'NPA detector poloidal angle',istatus)
            call ncaptc2(ncid,vid,'units',NCCHAR,7, &
                 'degrees',istatus)
         endif

            vid=ncvdef0(ncid,'nv_npa',NCLONG,0,0,istatus)
            call ncaptc2(ncid,vid,'long_name',NCCHAR,27, &
              'Number of NPA viewing cords',istatus)
            call check_err(istatus)

            vid=ncvdef0(ncid,'nen_npa',NCLONG,0,0,istatus)
            call ncaptc2(ncid,vid,'long_name',NCCHAR,40, &
              'Number of equispaced energies in spectra',istatus)
            call check_err(istatus)

            vid=ncvdef0(ncid,'npaproc',NCLONG,0,0,istatus)
            call ncaptc2(ncid,vid,'long_name',NCCHAR,41, &
              'Maximum number of CX processes considered',istatus)
            call check_err(istatus)

            vid=ncvdef0(ncid,'enmin_npa',NCDOUBLE,0,0,istatus)
            call ncaptc2(ncid,vid,'long_name',NCCHAR,29, &
              'NPA minimun energy in spectra',istatus)
            call ncaptc2(ncid,vid,'units',NCCHAR,3, &
              'keV',istatus)

            vid=ncvdef0(ncid,'enmax_npa',NCDOUBLE,0,0,istatus)
            call ncaptc2(ncid,vid,'long_name',NCCHAR,29, &
              'NPA maximum energy in spectra',istatus)
            call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                 'keV',istatus)

         vid=ncvdef2(ncid,'npa_process',NCCHAR,2,npaproc_dims,istatus)
         call ncaptc2(ncid,vid,'long_name0',NCCHAR,54, &
            '(1)cxh,(2)cxb4,(3)cxhe,(4)cxc,(5)radrecom, or notset()', &
             istatus)
         call ncaptc2(ncid,vid,'long_name1',NCCHAR,40, &
            'Indicates which CX process(es) are included:' ,istatus)
         call check_err(istatus)

         vid=ncvdef0(ncid,'atten_npa',NCCHAR,1,chardim,istatus)
         call ncaptc2(ncid,vid,'long_name0',NCCHAR,42, &
            'enabled, normal calculation of attenuation' ,istatus)
         call ncaptc2(ncid,vid,'long_name1',NCCHAR,40, &
            'disabled, for numerical testing purposes' ,istatus)
         call check_err(istatus)

         vid=ncvdef0(ncid,'ipronn',NCCHAR,1,chardim,istatus)
         call ncaptc2(ncid,vid,'long_name0',NCCHAR,39, &
            'disabled, default, zero neutral density' ,istatus)
         call ncaptc2(ncid,vid,'long_name1',NCCHAR,51, &
         'exp, ennl exponential falloff with dist from radmin' ,istatus)
         call ncaptc2(ncid,vid,'long_name1',NCCHAR,46, &
         'spline, give neutral density profiles vs ryain' ,istatus)
         call check_err(istatus)

         vid=ncvdef0(ncid,'en_',NCDOUBLE,(1),nendim,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,17, &
              'CX Neutral Energy',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,3, &
              'keV',istatus)

         vid=ncvdef0(ncid,'ennscal',NCDOUBLE,(1),npaproc_dim,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,51, &
           'Scale factor for density assoc with each CX process' &
           ,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,49, &
           'ennscal has been used in calc of the enn profiles' &
           ,istatus)

      vid=ncvdef2(ncid,'enn',NCDOUBLE,2,npaenn_dims,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,40, &
              'Density associated with each npa_process',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,6, &
              '/cm**3',istatus)

      vid=ncvdef2(ncid,'eflux_npa',NCDOUBLE,3,npa_dims,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,12, &
              'Neutral flux',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,19, &
              '#/cm**2/sec/ster/eV',istatus)

         vid=ncvdef2(ncid,'efluxt',NCDOUBLE,2,npa_dims(2),istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,23, &
              'Integrated neutral flux',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,16, &
              '#/cm**2/sec/ster',istatus)

      endif  ! On npa_diag .ne. "disabled"


!     Fusion reaction rate data
      vid=ncvdef0(ncid,'sigmamod',NCCHAR,1,chardim,istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'long_name0',NCCHAR,46, &
            'fusion rates diagnostic is disabled or enabled', &
            istatus)
      call ncaptc2(ncid,vid,'long_name1',NCCHAR,45, &
            'Data for Rx 1:4 is output (see cqlinput_help)' ,istatus)
      call check_err(istatus)

      if (sigmamod .ne. "disabled") then ! define netcdf names

         vid=ncvdef0(ncid,'isigmas',NCLONG,1,fourdim,istatus)
         call check_err(istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,40, &
              'Indicators of the particular rates calcd',istatus)

         vid=ncvdef0(ncid,'mmsv',NCLONG,0,0,istatus)
         call check_err(istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,44, &
              'Max order Legendre expansion of Rx cross-sec',istatus)

         vid=ncvdef0(ncid,'isigsgv1',NCLONG,0,0,istatus)
         call check_err(istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,47, &
              'Indicator of Rx due to gen distn with self, =0 ',istatus)

         vid=ncvdef0(ncid,'isigsgv2',NCLONG,0,0,istatus)
         call check_err(istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,40, &
              'Indicator of Bgrnd Max not included, =0 ',istatus)
         call check_err(istatus)

      endif  !  On sigmamod

!--------------------------
!     Time-dependent data
!--------------------------

      vid=ncvdef0(ncid,'time',NCDOUBLE,(1),tdim,istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7, &
                           'seconds',istatus)

      !YuP[2018-09-28], BH181112 added for 'lngshrtf' option,
      !for saving f() distr.func. at selected (nsave()) t steps only.
      vid=ncvdef0(ncid,'nsave',NCLONG,1,tsavedim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,41, &
                 'Selected time steps, n.eq.nsave(1:nsavet)',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,4,'none',istatus)
      
      if((netcdfshort.eq.'lngshrtf').and.isave.ne.0) then  !isave/nsavet set in tdchief
        !YuP[2019-06-08] Added  if()
        ! Define 'tsave' only for this case, otherwise 'tsave' will be empty
        vid=ncvdef0(ncid,'tsave',NCDOUBLE,(1),tsavedim,istatus)
        call ncaptc2(ncid,vid,'long_name',NCCHAR,41, &
                 'Times selected using n.eq.nsave(1:nsavet)',istatus)
        call ncaptc2(ncid,vid,'units',NCCHAR,7,'seconds',istatus)
      endif

      vid=ncvdef2(ncid,'den_e',NCDOUBLE,2,r00_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,16, &
                 'Electron density',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,6, &
                 '/cm**3',istatus)

      vid=ncvdef2(ncid,'density',NCDOUBLE,3,species_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,34, &
                 'Densities, general and Maxwellians',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,6, &
                 '/cm**3',istatus)

      vid=ncvdef2(ncid,'zeff',NCDOUBLE,2,r0_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,4, &
                 'Zeff',istatus)

      vid=ncvdef2(ncid,'temp',NCDOUBLE,3,species_dims,istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                 'keV',istatus)

      vid=ncvdef2(ncid,'energy',NCDOUBLE,3,species_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,23, &
                 'FSA Energy per particle',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                 'keV',istatus)

      if (ngen.eq.1) then
         vid=ncvdef2(ncid,'wpar',NCDOUBLE,2,r_dims,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,32, &
              'FSA Parallel Energy per particle',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,3, &
              'keV',istatus)

         vid=ncvdef2(ncid,'wperp',NCDOUBLE,2,r_dims,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,37, &
              'FSA Perpendicular Energy per particle',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,3, &
              'keV',istatus)
      else  !  ngen.ge.2
         vid=ncvdef2(ncid,'wpar',NCDOUBLE,3,rk_dims,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,32, &
              'FSA Parallel Energy per particle',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,3, &
              'keV',istatus)

         vid=ncvdef2(ncid,'wperp',NCDOUBLE,3,rk_dims,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,37, &
              'FSA Perpendicular Energy per particle',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,3, &
              'keV',istatus)
      endif  !  on ngen

      vid=ncvdef2(ncid,'elecfld',NCDOUBLE,2,r00_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,23, &
                 'Parallel Electric Field',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,8, &
                 'Volts/cm',istatus)

      vid=ncvdef2(ncid,'edreicer',NCDOUBLE,2,r0_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,45, &
                'E_D Dreicer elec fld, e.g., Kulsrud PRL(1973)',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,8, &
                 'Volts/cm',istatus)

      vid=ncvdef2(ncid,'runaway_rate',NCDOUBLE,2,r_dims,istatus)
      call ncaptc2(ncid,vid,'long_name0',NCCHAR,45, &
                'Runaway rate, determined from e flux off grid',istatus)
      call ncaptc2(ncid,vid,'long_name1',NCCHAR,39, &
                'Runaway rate = 1/n * dn/dt / nu_Kulsrud',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,8, &
                 'Unitless',istatus)

      vid=ncvdef2(ncid,'denra',NCDOUBLE,2,r_dims,istatus)
      call ncaptc2(ncid,vid,'long_name0',NCCHAR,31, &
                'Runaway FSA density above ucrit',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,6, &
                 '/cm**3',istatus)

      vid=ncvdef2(ncid,'curra',NCDOUBLE,2,r_dims,istatus)
      call ncaptc2(ncid,vid,'long_name0',NCCHAR,44, &
                'Runaway FSA parallel cur density above ucrit',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,10, &
                 'Amps/cm**2',istatus)

      vid=ncvdef2(ncid,'ucrit',NCDOUBLE,2,r_dims,istatus)
      call ncaptc2(ncid,vid,'long_name0',NCCHAR,38, &
                'Critical momentum per mass for runaway',istatus)
!BH071013      call ncaptc2(ncid,vid,'units',NCCHAR,6, &
!BH071013                 'cm/sec',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,19, &
                 'Normalized to vnorm',istatus)

!     Knockon electron data
      vid=ncvdef0(ncid,'knockon',NCCHAR,1,chardim,istatus)
      call ncaptc2(ncid,vid,'long_name0',NCCHAR,47, &
            'Knockon src of high energy elec is en-/disabled' ,istatus)

      if (knockon.ne."disabled") then
      vid=ncvdef2(ncid,'eoe0',NCDOUBLE,2,r_dims,istatus)
         call ncaptc2(ncid,vid,'long_name0',NCCHAR,39, &
              'Elecfld/Critical knockon electric field',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,8, &
              'Unitless',istatus)

      vid=ncvdef2(ncid,'srckotot',NCDOUBLE,2,r_dims,istatus)
         call ncaptc2(ncid,vid,'long_name0',NCCHAR,31, &
              'FSA Knockon source density rate',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,11, &
              '#/cm**3*sec',istatus)

      vid=ncvdef2(ncid,'denfl',NCDOUBLE,2,r_dims,istatus)
         call ncaptc2(ncid,vid,'long_name0',NCCHAR,34, &
              'FSA Elec Density from KO Reduced Distn',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,7, &
              '#/cm**3',istatus)
      endif  ! on knockon
!     In this section: we assume not both rdcmod and urfmod.ne.disabled
!     mrfn=number of modes =1, in subroutine rdc_multi.
      if (rdcmod.ne."disabled") then
         if (rdcmod.eq."enabled" .and. urfmod.ne."disabled") then
            WRITE(*,*)
            WRITE(*,*)'Warning: netcdfrw2 not set up for both'
            WRITE(*,*)'Warning: rdcmod and urfmod.ne.disabled.'
            WRITE(*,*)'Warning: urfmod output will overwrite rdcmod o/p'
            WRITE(*,*)
         endif
         vid=ncvdef2(ncid,'rfpwr',NCDOUBLE,3,rfpwr_dims,istatus)
         call ncaptc2(ncid,vid,'long_name0',NCCHAR,40, &
              'RF power densities (sorpw_rf(*,1:mrfn)):',istatus)
         call ncaptc2(ncid,vid,'long_name1',NCCHAR,41, &
              'Radially integrated: sorpw_rfi(*,mrfn+3)=',istatus)
         call ncaptc2(ncid,vid,'long_name2',NCCHAR,40, &
              'Summed rf+nbi pwr den: sorpwt(*,mrfn+2)=',istatus)
         call ncaptc2(ncid,vid,'long_name3',NCCHAR,39, &
              'Radially integrated: sorpwti(*,mrfn+3)=',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,37, &
              'Watts/cm**3, except Watts for sorpwti',istatus)

      endif  !On rdcmod

!
!... new Freya stuff
!
        WRITE(*,*) frmodp
      if (frmodp.eq."enabled") then

           vid=ncvdef2(ncid,'hibrz',NCDOUBLE,3,hibr_dims,istatus)
           call ncaptc2(ncid,vid,'long_name0',NCCHAR,29, &
               'Normalized hot ion birth rate',istatus)

           vid=ncvdef2(ncid,'sorpw_nbi',NCDOUBLE,2,sorpw_dims,istatus)
           call ncaptc2(ncid,vid,'long_name0',NCCHAR,28, &
               'NBI+FUS Source Power Density',istatus)
           call ncaptc2(ncid,vid,'units',NCCHAR,10, &
               'Watts/cm^3',istatus)

           vid=ncvdef2(ncid,'sorpw_nbii',NCDOUBLE,2,sorpw_dims,istatus)
           call ncaptc2(ncid,vid,'long_name0',NCCHAR,38, &
               'Radially Integrated Power from NBI+FUS',istatus)
           call ncaptc2(ncid,vid,'units',NCCHAR,5, &
               'Watts',istatus)
          !YuP[06-2016] Now the corresponding arrays contain
          ! a COMBINED NBI+FUSproduct source.
          ! But the netcdf names like 'sorpw_nbi' are not changed
          ! because Python plotting script must be changed in such case.

      endif
!

      if (urfmod.ne."disabled") then
         vid=ncvdef0(ncid,'mrfn',NCLONG,0,0,istatus) ! YuP92017]added
         call ncaptc2(ncid,vid,'long_name',NCCHAR,59, &
         'number of rf modes (sum over all wave types and all nharms)', &
         istatus)

         vid=ncvdef2(ncid,'powurf',NCDOUBLE,2,dims_powurf,istatus)
         call ncaptc2(ncid,vid,'long_name0',NCCHAR,33, &
              'URF power in each mode, and total',istatus)
         call ncaptc2(ncid,vid,'long_name1',NCCHAR,41, &
              'powurf(1:mrfn)=power for individual modes',istatus)
         call ncaptc2(ncid,vid,'long_name2',NCCHAR,38, &
              'powurf(mrfn+1)=power summed over modes',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,5, &
              'Watts',istatus)

         vid=ncvdef2(ncid,'rfpwr',NCDOUBLE,3,rfpwr_dims,istatus)
         call ncaptc2(ncid,vid,'long_name0',NCCHAR,37, &
              'RF power densities (powrf(*,1:mrfn)):',istatus)
         call ncaptc2(ncid,vid,'long_name1',NCCHAR,45, &
              'rfpwr(*,1:mrfn)=pwr den from individual modes',istatus)
         call ncaptc2(ncid,vid,'long_name2',NCCHAR,24, &
              'rfpwr(*,mrfn+1)=powrft()',istatus)
         call ncaptc2(ncid,vid,'long_name3',NCCHAR,40, &
              'Summed rf+nbi pwr den: sorpwt(*,mrfn+2)=',istatus)
         call ncaptc2(ncid,vid,'long_name4',NCCHAR,39, &
              'Radially integrated: sorpwti(*,mrfn+3)=',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,37, &
              'Watts/cm**3, except Watts for sorpwti',istatus)

         vid=ncvdef2(ncid,'powrf',NCDOUBLE,3,dims_powrf,istatus)
         call ncaptc2(ncid,vid,'long_name0',NCCHAR,47, &
              'RF power densities due to mode (or harmonic for',istatus)
         call ncaptc2(ncid,vid,'long_name1',NCCHAR,39, &
              'nharms.gt.1 cases [powrf(lrza,nmodsa)])',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,11, &
              'Watts/cm**3',istatus)
         vid=ncvdef0(ncid,'nrfspecies',NCLONG,1,dims_powrf(2),istatus) !nmodsa
         !YuP[11-2017] added
         call ncaptc2(ncid,vid,'long_name',NCCHAR,60, &
         'nrfspecies(nmodsa)= general species index for each wave type', &
         istatus)
         call check_err(istatus)

         vid=ncvdef2(ncid,'powrfl',NCDOUBLE,3,dims_powrf,istatus)
         call ncaptc2(ncid,vid,'long_name0',NCCHAR,46, &
              'RF power densities due to salphal(1:nmodsdim)',istatus)
         call ncaptc2(ncid,vid,'long_name1',NCCHAR,40, &
              '(For multi-harmonic or multi-mode cases)',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,11, &
              'Watts/cm**3',istatus)

         vid=ncvdef2(ncid,'powurfl',NCDOUBLE,2,dims_powrf(2),istatus)
         call ncaptc2(ncid,vid,'long_name0',NCCHAR,46, &
              'Tot RF power absrbd due to salphal(1:nmodsdim)',istatus)
         call ncaptc2(ncid,vid,'long_name1',NCCHAR,40, &
              '(For multi-harmonic or multi-mode cases)',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,5, &
              'Watts',istatus)

         vid=ncvdef2(ncid,'powrfc',NCDOUBLE,3,dims_powrf,istatus)
         call ncaptc2(ncid,vid,'long_name0',NCCHAR,50, &
           'Coll RF power densities due to salphac(1:nmodsdim)',istatus)
         call ncaptc2(ncid,vid,'long_name1',NCCHAR,40, &
              '(For multi-harmonic or multi-mode cases)',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,11, &
              'Watts/cm**3',istatus)

         vid=ncvdef2(ncid,'powurfc',NCDOUBLE,2,dims_powrf(2),istatus)
         call ncaptc2(ncid,vid,'long_name0',NCCHAR,49, &
           'Tot Coll RF pwr absrbd due to salphac(1:nmodsdim)',istatus)
         call ncaptc2(ncid,vid,'long_name1',NCCHAR,40, &
              '(For multi-harmonic or multi-mode cases)',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,5, &
              'Watts',istatus)

         vid=ncvdef2(ncid,'powrft',NCDOUBLE,2,dims_powrft,istatus)
         call ncaptc2(ncid,vid,'long_name0',NCCHAR,50, &
           'RF power densities summed over modes or harmonics,',istatus)
         call ncaptc2(ncid,vid,'long_name1',NCCHAR,43, &
              'due to urf, collisional and add. linear abs',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,11, &
              'Watts/cm**3',istatus)

      endif  !On urfmod


      vid=ncvdef2(ncid,'curtor',NCDOUBLE,2,r0_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,36, &
                 'Toroidal current density at min B pt',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,10, &
                 'Amps/cm**2',istatus)

      vid=ncvdef2(ncid,'ccurtor',NCDOUBLE,2,r0_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,32, &
                 'Area Integrated toroidal current',istatus)
      call ncaptc2(ncid,vid,'long_name1',NCCHAR,39, &
                 'accounting for pol variation of tor cur',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,4, &
                 'Amps',istatus)

      vid=ncvdef2(ncid,'curpol',NCDOUBLE,2,r0_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,36, &
                 'Poloidal current density at min B pt',istatus)
       call ncaptc2(ncid,vid,'units',NCCHAR,10, &
                 'Amps/cm**2',istatus)

      vid=ncvdef2(ncid,'ccurpol',NCDOUBLE,2,r0_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,44, &
                 'Integrated poloidal current density at min B',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,4, &
                 'Amps',istatus)


      if (kelecg.ne.0) then ! e as gen.species
         vid=ncvdef2(ncid,'currm_e',NCDOUBLE,2,r0_dims,istatus)
      else ! ion as general species
         vid=ncvdef2(ncid,'currm_i',NCDOUBLE,2,r0_dims,istatus)
      endif

      call ncaptc2(ncid,vid,'long_name0',NCCHAR,38, &
                 'Parallel elec current density at min B',istatus)
      call ncaptc2(ncid,vid,'long_name1',NCCHAR,43, &
                 'Electrons, or first gen species if kelecg=0',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,10, &
                 'Amps/cm**2',istatus)

      if (ngen.eq.1) then
         vid=ncvdef2(ncid,'curr',NCDOUBLE,2,r0_dims,istatus)
         call ncaptc2(ncid,vid,'long_name0',NCCHAR,20, &
              'FSA Parallel current',istatus)
         call ncaptc2(ncid,vid,'long_name1',NCCHAR,44, &
              'i.e., Par curr per poloidal area between FSs',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,10, &
              'Amps/cm**2',istatus)
         !YuP[07-31-2014] Added:
         vid=ncvdef2(ncid,'energym',NCDOUBLE,2,r_dims,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,40, &
                 'Energy per particle, from f0 at midplane',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                 'keV',istatus)
      else  !  ngen.ge.2
         vid=ncvdef2(ncid,'curr',NCDOUBLE,3,r0k_dims,istatus)
         call ncaptc2(ncid,vid,'long_name0',NCCHAR,20, &
              'FSA Parallel current',istatus)
         call ncaptc2(ncid,vid,'long_name1',NCCHAR,44, &
              'i.e., Par curr per poloidal area between FSs',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,10, &
              'Amps/cm**2',istatus)
         !YuP[07-31-2014] Added:
         vid=ncvdef2(ncid,'energym',NCDOUBLE,3,r0k_dims,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,40, &
                 'Energy per particle, from f0 at midplane',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                 'keV',istatus)

      endif !ngen

      vid=ncvdef2(ncid,'restp',NCDOUBLE,2,r0_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,19, &
                 '<E_phi/R>/<j_phi/R>',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,12, &
                 'cgs, seconds',istatus)

      vid=ncvdef2(ncid,'restnp',NCDOUBLE,2,r0_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,45, &
                'neoclassical resist <E_parall*B>/<j_parall*B>',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,12, &
                 'cgs, seconds',istatus)

      vid=ncvdef2(ncid,'sptzrp',NCDOUBLE,2,r0_dims,istatus)
      call ncaptc2(ncid,vid,'long_name0',NCCHAR,41, &
                 'Spitzer resistivity, incl Zeff dependence',istatus)
      call ncaptc2(ncid,vid,'long_name1',NCCHAR,41, &
                 'Eq. 4.2-77 of ONETWO Manual, Eq. 5.66 H&H',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,12, &
                 'cgs, seconds',istatus)

      vid=ncvdef2(ncid,'rovsc',NCDOUBLE,2,r0_dims,istatus)
      call ncaptc2(ncid,vid,'long_name0',NCCHAR,31, &
                 'Connor resistivity over Spitzer',istatus)
      call ncaptc2(ncid,vid,'long_name1',NCCHAR,42, &
                 'J.W. Connor et al, Nucl Fus 13, 211 (1973)',istatus)

      vid=ncvdef2(ncid,'rovsc_hi',NCDOUBLE,2,r0_dims,istatus)
      call ncaptc2(ncid,vid,'long_name0',NCCHAR,39, &
                 'Connor resistivity over Spitzer, hi eps',istatus)
      call ncaptc2(ncid,vid,'long_name1',NCCHAR,42, &
                 'J.W. Connor et al, Nucl Fus 13, 211 (1973)',istatus)

      vid=ncvdef2(ncid,'zreskim',NCDOUBLE,2,r0_dims,istatus)
      call ncaptc2(ncid,vid,'long_name0',NCCHAR,44, &
                 'Hirshman/Sigmar/Kim resistivity over Spitzer',istatus)
      call ncaptc2(ncid,vid,'long_name1',NCCHAR,43, &
                 'Reference is Hir/Sig (~1988) and Kim theses',istatus)

      vid=ncvdef2(ncid,'taueeh',NCDOUBLE,2,r0_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,36, &
                 'Hinton-Hazeltine(Eq 5.4)-ONETWO-taue',istatus)

      vid=ncvdef2(ncid,'nuestar',NCDOUBLE,2,r0_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,25, &
                 'ONETWO Manual, Eq. 4.2-30',istatus)

      vid=ncvdef2(ncid,'powers',NCDOUBLE,4,powers_dims,istatus)
      call ncaptc2(ncid,vid,'long_name0',NCCHAR,55, &
          'Component by component FSA powers to gen species k vs t', &
           istatus)
      call ncaptc2(ncid,vid,'long_name1',NCCHAR,53, &
        'powers(*,1,k,t)=due to collisions with Maxw electrons',istatus)
      call ncaptc2(ncid,vid,'long_name2',NCCHAR,48, &
         'powers(*,2,k,t)=due to collisions with Maxw ions',istatus)
      call ncaptc2(ncid,vid,'long_name3',NCCHAR,25, &
         'powers(*,3,k,t)=Ohmic E.v',istatus)
      call ncaptc2(ncid,vid,'long_name4',NCCHAR,51, &
         'powers(*,4,k,t)=due to collisions with general spec',istatus)
      call ncaptc2(ncid,vid,'long_name5',NCCHAR,24, &
         'powers(*,5,k,t)=RF power',istatus)
      call ncaptc2(ncid,vid,'long_name6',NCCHAR,35, &
         'powers(*,6,k,t)=Ion particle source',istatus)
      call ncaptc2(ncid,vid,'long_name7',NCCHAR,34, &
         'powers(*,7,k,t)=losses by lossmode',istatus)
      call ncaptc2(ncid,vid,'long_name8',NCCHAR,33, &
         'powers(*,8,k,t)=losses by torloss',istatus)
      call ncaptc2(ncid,vid,'long_name9',NCCHAR,30, &
         'powers(*,9,k,t)=Runaway losses',istatus)
      call ncaptc2(ncid,vid,'long_name10',NCCHAR,45, &
         'powers(*,10,k,t)=Synchrotron radiation losses',istatus)
      call ncaptc2(ncid,vid,'long_name11',NCCHAR,39, &
         'powers(*,11,k,t)=Setting neg. j to zero',istatus)
      call ncaptc2(ncid,vid,'long_name12',NCCHAR,40, &
         'powers(*,12,k,t)=Phenomenological losses',istatus)
      call ncaptc2(ncid,vid,'long_name13',NCCHAR,22, &
         'powers(*,13,k,t)=Total',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,11, &
         'Watts/cm**3',istatus)

      vid=ncvdef2(ncid,'powers_int',NCDOUBLE,3,powers_dims(2),istatus)
      call ncaptc2(ncid,vid,'long_name0',NCCHAR,58, &
          'Vol int of FSA powers, respectively, to gen species k vs t', &
           istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5, &
         'Watts',istatus)

      if (sigmamod.eq.'enabled') then ! define netcdf names
      vid=ncvdef2(ncid,'sigftt',NCDOUBLE,2,fus_dims(2),istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,29, &
              'Total fusion Rx rates vs time',istatus)
      endif


      if (netcdfshort.eq.'enabled') then
!        Do nothing: no storage defined.
      elseif ( (netcdfshort.eq.'longer_f').or. &
               (netcdfshort.eq.'lngshrtf')     ) then ! define storage

         if (ngen.eq.1) then    !maintaining backwards compatability
            vid=ncvdef2(ncid,'f',NCDOUBLE,4,dimsf(1:4),istatus) ! ngen=1
            ! here dimsf={ydim,xdim,rdim,tdim(or tsavedim)}
            call ncaptc2(ncid,vid,'long_name',NCCHAR,21, &
                 'Distribution function',istatus)
            call ncaptc2(ncid,vid,'units',NCCHAR,28, &
                 'vnorm**3/(cm**3*(cm/sec)**3)',istatus)
         else  !ngen.ge.2
            vid=ncvdef2(ncid,'f',NCDOUBLE,5,dimsg(1:5),istatus) ! ngen>1
                !Additional dim included for ngen.gt.1 cases
            ! here dimsg={ydim,xdim,rdim,gdim,tdim(or tsavedim)}
            call ncaptc2(ncid,vid,'long_name',NCCHAR,21, &
                 'Distribution function',istatus)
            call ncaptc2(ncid,vid,'units',NCCHAR,28, &
                 'vnorm**3/(cm**3*(cm/sec)**3)',istatus)
            call ncaptc2(ncid,vid,'comment',NCCHAR,44, &
                 'Additional dimension added for multi-species',istatus)
         endif  !on ngen

      else     !disabled, Standard o/p: f at last time step

         if (ngen.eq.1) then    !maintaining backwards compatability
            vid=ncvdef2(ncid,'f',NCDOUBLE,3,dimsf(1:3),istatus)
            ! here dimsf={ydim,xdim,rdim}
            call ncaptc2(ncid,vid,'long_name',NCCHAR,21, &
                 'Distribution function',istatus)
            call ncaptc2(ncid,vid,'units',NCCHAR,28, &
                 'vnorm**3/(cm**3*(cm/sec)**3)',istatus)
            call ncaptc2(ncid,vid,'comment',NCCHAR,39, &
                 'Facility set up only for single species',istatus)
         else  !ngen.ge.2
            vid=ncvdef2(ncid,'f',NCDOUBLE,4,dimsg(1:4),istatus)
                !Additional dim included for ngen.gt.1 cases
            ! here dimsg={ydim,xdim,rdim,gdim}
            call ncaptc2(ncid,vid,'long_name',NCCHAR,21, &
                 'Distribution function',istatus)
            call ncaptc2(ncid,vid,'units',NCCHAR,28, &
                 'vnorm**3/(cm**3*(cm/sec)**3)',istatus)
            call ncaptc2(ncid,vid,'comment',NCCHAR,44, &
                 'Additional dimension added for multi-species',istatus)
         endif  !on ngen

      endif  !on netcdfshort : storage defined


!BH100912:  Option to output specific curr and rf pwr at each time step
      if (netcdfshort.eq.'long_jp') then

      if (ngen.eq.1) then
      vid=ncvdef2(ncid,'currv',NCDOUBLE,3,currv_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,44, &
                 'Specific Current Density j_u(u) at each step',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,43, &
                 'Amps/cm^2 (int:0,1 over dx =current density',istatus)
      call ncaptc2(ncid,vid,'comment',NCCHAR,34, &
                 'Facility set up for single species',istatus)

      vid=ncvdef2(ncid,'pwrrf',NCDOUBLE,3,currv_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,47, &
              'Specific RF Power Density pwrrf(u) at each step',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,44, &
                'W/cm^3 (int:0,1 over dx =RF power density',istatus)
      call ncaptc2(ncid,vid,'comment',NCCHAR,34, &
                 'Facility set up for single species',istatus)
      else  !ngen.ge.2
      vid=ncvdef2(ncid,'currv',NCDOUBLE,4,currv_dimsg,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,44, &
                 'Specific Current Density j_u(u) at each step',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,43, &
                 'Amps/cm^2 (int:0,1 over dx =current density',istatus)
      call ncaptc2(ncid,vid,'comment',NCCHAR,23, &
                 'Setup for multi-species',istatus)

      vid=ncvdef2(ncid,'pwrrf',NCDOUBLE,4,currv_dimsg,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,47, &
              'Specific RF Power Density pwrrf(u) at each step',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,41, &
                'W/cm^3 (int:0,1 over dx =RF power density',istatus)
      call ncaptc2(ncid,vid,'comment',NCCHAR,23, &
                 'Setup for multi-species',istatus)
      endif  !on ngen

      else  ! Setup for last time step only

      if (ngen.eq.1) then
      vid=ncvdef2(ncid,'currv',NCDOUBLE,2,currv_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,31, &
                 'Specific Current Density j_u(u)',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,44, &
                 'Amps/cm^2 (int:0,1 over dx =current density)',istatus)
      call ncaptc2(ncid,vid,'comment',NCCHAR,34, &
                 'Facility set up for single species',istatus)

      vid=ncvdef2(ncid,'pwrrf',NCDOUBLE,2,currv_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,34, &
                 'Specific RF Power Density pwrrf(u)',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,42, &
                'W/cm^3 (int:0,1 over dx =RF power density)',istatus)
      call ncaptc2(ncid,vid,'comment',NCCHAR,34, &
                 'Facility set up for single species',istatus)
      else  !ngen.ge.2
      vid=ncvdef2(ncid,'currv',NCDOUBLE,3,currv_dimsg,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,31, &
                 'Specific Current Density j_u(u)',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,44, &
                 'Amps/cm^2 (int:0,1 over dx =current density)',istatus)
      call ncaptc2(ncid,vid,'comment',NCCHAR,23, &
                 'Setup for multi-species',istatus)

      vid=ncvdef2(ncid,'pwrrf',NCDOUBLE,3,currv_dimsg,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,34, &
                 'Specific RF Power Density pwrrf(u)',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,42, &
                'W/cm^3 (int:0,1 over dx =RF power density)',istatus)
      call ncaptc2(ncid,vid,'comment',NCCHAR,23, &
                 'Setup for multi-species',istatus)
      endif  !on ngen


      endif  !  On netcdfshort.eq.'long_jp'


!--------------------------
!     Last time step data
!--------------------------

!BH100912
!$$$      vid=ncvdef2(ncid,'currv',NCDOUBLE,2,currv_dims,istatus)
!$$$      call ncaptc2(ncid,vid,'long_name',NCCHAR,31,
!$$$     +           'Specific Current Density j_u(u)',istatus)
!$$$      call ncaptc2(ncid,vid,'units',NCCHAR,44,
!$$$     +           'StatA/cm^2 (int:0,1 over dx =current density',istatus)
!$$$      call ncaptc2(ncid,vid,'comment',NCCHAR,39,
!$$$     +           'Facility set up only for single species',istatus)
!$$$
!$$$      vid=ncvdef2(ncid,'pwrrf',NCDOUBLE,2,currv_dims,istatus)
!$$$      call ncaptc2(ncid,vid,'long_name',NCCHAR,34,
!$$$     +           'Specific RF Power Density pwrrf(u)',istatus)
!$$$      call ncaptc2(ncid,vid,'units',NCCHAR,45,
!$$$     +          'StatA/cm^2 (int:0,1 over dx =RF power density',istatus)
!$$$      call ncaptc2(ncid,vid,'comment',NCCHAR,39,
!$$$     +           'Facility set up only for single species',istatus)


      if (sigmamod.eq.'enabled') then  ! define netcdf names
         vid=ncvdef0(ncid,'fuspwrvt',NCDOUBLE,(1),fus_dims(2),istatus)
         call ncaptc2(ncid,vid,'long_name0',NCCHAR,38, &
              'Total fusion power, for four reactions',istatus)
         call ncaptc2(ncid,vid,'long_name1',NCCHAR,18, &
              'At final time step',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,5, &
              'Watts',istatus)

         vid=ncvdef2(ncid,'fuspwrv',NCDOUBLE,2,fus_dims(1),istatus)
         call ncaptc2(ncid,vid,'long_name0',NCCHAR,46, &
              'Fusion power versus radius, for four reactions',istatus)
         call ncaptc2(ncid,vid,'long_name1',NCCHAR,18, &
              'At final time step',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,11, &
              'Watts/cm**3',istatus)
      endif


!.......................................................................
!cl    1.1.4 end the define-mode and start the data-mode
!     p. 51-2 of manual

      istatus= NF_ENDDEF(ncid) !-YuP: NetCDF-f77
      call check_err(istatus)

!-----------------------------------------------------------------------
!cl    1.2 Write data:   First call to netcdfrw2
!

! --- set the time-step counter ==> numrec1
      numrec1=1
      numrecsave=1   ! for (netcdfshort.eq.'lngshrtf')
      start(4)=numrec1 !=1 here;  used for 'time', etc.
      start1(4)=numrec1
      startg(5)=numrec1 ! here numrec1=1 =numrecsave=1
      start_rk(3)=numrec1
      start_r0k(3)=numrec1
      start_powurf(2)=numrec1
      start_rfpwr(3)=numrec1
      start_powrf(3)=numrec1
      start_powrft(2)=numrec1
      start_powers(4)=numrec1
      start_fus(3)=numrec1
      start_xr(3)=numrec1
      start_npa(3)=numrec1

! --- initialize data file ---
!     First get variable_id:
!             nf_inq_varid(netcdf_id,variable_name,integer_info)
!     Then write data with nc variable_put:
!             nf_put_var... ()
!
!     function nf_put_var...(ncid,variable_id,index,val)
!     p. 52,, 54-65 of netcdf-3 manual
!

      istatus= NF_INQ_VARID(ncid,'version',vid)  !-YuP: NetCDF-f77 get vid
      ll=length_char(version)
      call ncvptc0(ncid,vid,1,ll,version,ll,istatus)

      istatus= NF_INQ_VARID(ncid,'mnemonic',vid)  !-YuP: NetCDF-f77 get vid
      ll=length_char(mnemonic)
      call ncvptc0(ncid,vid,1,ll,mnemonic,ll,istatus)

      istatus= NF_INQ_VARID(ncid,'ampfmod',vid)
      ll=length_char(frmodp)
      call ncvptc0(ncid,vid,1,ll,ampfmod,ll,istatus)

      istatus= NF_INQ_VARID(ncid,'urfmod',vid)
      ll=length_char(urfmod)
      call ncvptc0(ncid,vid,1,ll,urfmod,ll,istatus)

      istatus= NF_INQ_VARID(ncid,'rdcmod',vid)
      ll=length_char(rdcmod)
      call ncvptc0(ncid,vid,1,ll,rdcmod,ll,istatus)

      istatus= NF_INQ_VARID(ncid,'frmod',vid)
      ll=length_char(frmodp)
      call ncvptc0(ncid,vid,1,ll,frmodp,ll,istatus)

      istatus= NF_INQ_VARID(ncid,'beamplse',vid)
      ll=length_char(beamplsep)
      call ncvptc0(ncid,vid,1,ll,beamplsep,ll,istatus)

      istatus= NF_INQ_VARID(ncid,'transp',vid)
      ll=length_char(transp)
      call ncvptc0(ncid,vid,1,ll,transp,ll,istatus)

      istatus= NF_INQ_VARID(ncid,'tavg',vid)
      ll=length_char(tavg)
      call ncvptc0(ncid,vid,1,ll,tavg,ll,istatus)

      istatus= NF_INQ_VARID(ncid,'f4d_out',vid)
      ll=length_char(f4d_out)
      call ncvptc0(ncid,vid,1,ll,f4d_out,ll,istatus)

      istatus= NF_INQ_VARID(ncid,'netcdfshort',vid)
      ll=length_char(netcdfshort)
      call ncvptc0(ncid,vid,1,ll,netcdfshort,ll,istatus)

      istatus= NF_INQ_VARID(ncid,'eqdskin',vid)  !-YuP: NetCDF-f77 get vid
      ll=length_char(eqdskin)
      call ncvptc0(ncid,vid,1,ll,eqdskin,ll,istatus)

      istatus= NF_INQ_VARID(ncid,'ngen',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_int0(ncid,vid,1,1,(ngen),istatus)

      istatus= NF_INQ_VARID(ncid,'ntotal',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_int0(ncid,vid,1,1,(ntotal),istatus)

      istatus= NF_INQ_VARID(ncid,'kspeci',vid)  !-YuP: NetCDF-f77 get vid
      !call ncvptc2(ncid,vid,start(1:3),kspeci_count(1:3),kspeci(1:2,1:ntotal),8,istatus)
      !??? istatus= NF_PUT_VAR_TEXT(ncid,vid,kspeci) !YuP: does not work, either
      !Note: start(1:3)= {1,1,1}
      !kspeci_count(1:3)={chardim,twodim,kdim} <--> {8,2,ntotal}
      !kspeci(1:2,1:ntotala) is character(len=8)

      istatus= NF_INQ_VARID(ncid,'bnumb',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,(1),ntotal,bnumb,istatus)

      istatus= NF_INQ_VARID(ncid,'fmass',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,(1),ntotal,fmass,istatus)

      istatus= NF_INQ_VARID(ncid,'lrzmax',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_int0(ncid,vid,1,1,(lrzmax),istatus)

      istatus= NF_INQ_VARID(ncid,'radcoord',vid)  !-YuP: NetCDF-f77 get vid
      call ncvptc0(ncid,vid,1,8,radcoord,8,istatus)

      istatus= NF_INQ_VARID(ncid,'rya',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,(1),lrzmax,rya(1),istatus)

      istatus= NF_INQ_VARID(ncid,'Rp',vid)  ! rpcon(1:lrz) array
      call ncvpt_doubl0(ncid,vid,(1),lrzmax,rpcon,istatus)
      call check_err(istatus)

      istatus= NF_INQ_VARID(ncid,'Rm',vid)  ! rmcon(1:lrz) array
      call ncvpt_doubl0(ncid,vid,(1),lrzmax,rmcon,istatus)
      call check_err(istatus)


      istatus= NF_INQ_VARID(ncid,'rhomax',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,1,1,(rhomax),istatus)

      istatus= NF_INQ_VARID(ncid,'radmaj',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,1,1,(radmaj),istatus)

      istatus= NF_INQ_VARID(ncid,'rpmconz',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,1,lrzmax+1,rpmconz(0),istatus)

      istatus= NF_INQ_VARID(ncid,'btor',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,1,1,(btor),istatus)

      istatus= NF_INQ_VARID(ncid,'toteqd',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,(1),1,toteqd/3.e9,istatus)

      istatus= NF_INQ_VARID(ncid,'rgeomp',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,1,1,(rgeomp),istatus)

      istatus= NF_INQ_VARID(ncid,'r0geomp',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,1,1,(r0geomp),istatus)

      istatus= NF_INQ_VARID(ncid,'rmag',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,1,1,(rmag),istatus)

      istatus= NF_INQ_VARID(ncid,'zmag',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,1,1,(zmag),istatus)

      istatus= NF_INQ_VARID(ncid,'eqsym',vid)  !-YuP: NetCDF-f77 get vid
      call ncvptc0(ncid,vid,1,8,eqsym,8,istatus)

      istatus= NF_INQ_VARID(ncid,'zshift',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,1,1,(zshift),istatus)

      istatus= NF_INQ_VARID(ncid,'eps0',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,1,1,(eps0),istatus)

      istatus= NF_INQ_VARID(ncid,'elong',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,(1),1,zgeomp/rgeomp,istatus)

      istatus= NF_INQ_VARID(ncid,'area',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,(1),lrzmax,area(1),istatus)

      istatus= NF_INQ_VARID(ncid,'darea',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,(1),lrzmax,darea(1),istatus)

      istatus= NF_INQ_VARID(ncid,'vol',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,(1),lrzmax,vol(1),istatus)

      istatus= NF_INQ_VARID(ncid,'dvol',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,(1),lrzmax,dvol(1),istatus)

      istatus= NF_INQ_VARID(ncid,'equilpsi',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,(1),lrzmax,equilpsi(1),istatus)

      istatus= NF_INQ_VARID(ncid,'psivalm',vid)
      call ncvpt_doubl0(ncid,vid,(1),lrzmax,psivalm(1),istatus)

      istatus= NF_INQ_VARID(ncid,'psimag',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,(1),1,(psimag),istatus)

      istatus= NF_INQ_VARID(ncid,'psilim',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,(1),1,(psilim),istatus)

      istatus= NF_INQ_VARID(ncid,'dpsi',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,(1),lrzmax,dpsi(1),istatus)

      istatus= NF_INQ_VARID(ncid,'h_r',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,(1),lrzmax,h_r(1),istatus)

      istatus= NF_INQ_VARID(ncid,'qsafety',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,(1),lrzmax,qsafety(1),istatus)

      if (eqmod.eq."enabled") then
         istatus= NF_INQ_VARID(ncid,'curreq',vid)  !-YuP: NetCDF-f77 get vid
         call ncvpt_doubl0(ncid,vid,(1),lrzmax,curreq(1),istatus)
      endif

      istatus= NF_INQ_VARID(ncid,'lrz',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_int0(ncid,vid,1,1,(lrz),istatus)

      istatus= NF_INQ_VARID(ncid,'lrindx',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_int0(ncid,vid,1,lrz,setup0%lrindx(1),istatus)

      istatus= NF_INQ_VARID(ncid,'jx',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_int0(ncid,vid,1,1,(jx),istatus)

      istatus= NF_INQ_VARID(ncid,'x',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,(1),jx,x,istatus)

      istatus= NF_INQ_VARID(ncid,'enerkev',vid)
      call ncvpt_doubl0(ncid,vid,(1),jx,enerkev,istatus)

      istatus= NF_INQ_VARID(ncid,'uoc',vid)
      call ncvpt_doubl0(ncid,vid,(1),jx,uoc,istatus)

      istatus= NF_INQ_VARID(ncid,'dx',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,(1),jx,dx,istatus)

      istatus= NF_INQ_VARID(ncid,'cint2',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,(1),jx,cint2,istatus)

      istatus= NF_INQ_VARID(ncid,'vnorm',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,(1),1,(vnorm),istatus)

      istatus= NF_INQ_VARID(ncid,'enorm',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,(1),1,(enorm),istatus)

      istatus= NF_INQ_VARID(ncid,'iy',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_int0(ncid,vid,1,1,(iy),istatus)

      call pack21(y,1,iy,1,lrors,wkpack,iy,lrors)
      istatus= NF_INQ_VARID(ncid,'y',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start,y_count,wkpack,istatus)

      call pack21(dy,1,iy,1,lrors,wkpack,iy,lrors)
      istatus= NF_INQ_VARID(ncid,'dy',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start,y_count,wkpack,istatus)

      call pack21(cynt2,1,iy,1,lrors,wkpack,iy,lrors)
      istatus= NF_INQ_VARID(ncid,'cynt2',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start,y_count,wkpack,istatus)

      istatus= NF_INQ_VARID(ncid,'iy_',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_int0(ncid,vid,1,lrz,iy_,istatus)

      istatus= NF_INQ_VARID(ncid,'itl',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_int0(ncid,vid,1,lrz,itl_,istatus)

      istatus= NF_INQ_VARID(ncid,'itu',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_int0(ncid,vid,1,lrz,itu_,istatus)

      istatus= NF_INQ_VARID(ncid,'lz',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_int0(ncid,vid,1,1,(lz),istatus)

      call pack21(z,1,lza,1,lrzmax,wkpack,lz,lrzmax)
      istatus= NF_INQ_VARID(ncid,'z',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start,z_count(2),wkpack,istatus)

      call pack21(dz,1,lza,1,lrzmax,wkpack,lz,lrzmax)
      istatus= NF_INQ_VARID(ncid,'dz',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start,z_count(2),wkpack,istatus)

      call pack21(solrz,1,lza,1,lrzmax,wkpack,lz,lrzmax)
      istatus= NF_INQ_VARID(ncid,'solrz',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start,z_count(2),wkpack,istatus)

      call pack21(solzz,1,lza,1,lrzmax,wkpack,lz,lrzmax)
      istatus= NF_INQ_VARID(ncid,'solzz',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start,z_count(2),wkpack,istatus)

      call pack21(pol,1,lza,1,lrzmax,wkpack,lz,lrzmax)
      istatus= NF_INQ_VARID(ncid,'pol',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start,z_count(2),wkpack,istatus)

      call pack21(bbpsi,1,lza,1,lrzmax,wkpack,lz,lrzmax)
      istatus= NF_INQ_VARID(ncid,'bbpsi',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start,z_count(2),wkpack,istatus)

      call ipack21(imax,1,lza,1,lrzmax,iwkpack,lz,lrzmax) !YuP[2019-06-19] was wkpack. BUG?
      istatus= NF_INQ_VARID(ncid,'imax',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_int2(ncid,vid,start,z_count(2),iwkpack,istatus)

      call ipack21(lmax,1,iy,1,lrzmax,item1,iy,lrzmax)
      istatus= NF_INQ_VARID(ncid,'lmax',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_int2(ncid,vid,start,tau_count,item1,istatus)

      call pack21(zboun,1,iy,1,lrzmax,wkpack,iy,lrzmax)
      istatus= NF_INQ_VARID(ncid,'zboun',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start,tau_count,wkpack,istatus)

      istatus= NF_INQ_VARID(ncid,'zmaxpsi',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start,tau_count(2),zmaxpsi(1),istatus)

      call pack21(tau,1,iy,1,lrzmax,wkpack,iy,lrzmax)
      istatus= NF_INQ_VARID(ncid,'tau',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start,tau_count,wkpack,istatus)

      istatus= NF_INQ_VARID(ncid,'dtau',vid)  !-YuP: NetCDF-f77 get vid
      do ll=1,lrzmax
!BH         do l=1,lz
!BH            do i=1,iy
!BH               temp1(i,l)=dtau(i,l,ll)
!BH            enddo
!BH         enddo
!BH         call pack21(temp1,0,iyp1,0,jxp1,wkpack,iy,lz)
!BH  Actually, with dynamic dimensioning to exact size, can
!BH  leave out following step:
!BH         call pack21(dtau(1,1,ll),1,iy,1,lz,wkpack,iy,lz)
         start1(3)=ll
!BH         call ncvpt_doubl2(ncid,vid,start1,z_count1,wkpack,istatus)
         call ncvpt_doubl2(ncid,vid,start1(1:3),z_count1(1:3),dtau(1:iy,1:lz,ll), &
              istatus)
      enddo

      istatus= NF_INQ_VARID(ncid,'beampon',vid)
      call ncvpt_doubl0(ncid,vid,(1),1,(beamponp),istatus)

      istatus= NF_INQ_VARID(ncid,'beampoff',vid)
      call ncvpt_doubl0(ncid,vid,(1),1,(beampoffp),istatus)

      istatus= NF_INQ_VARID(ncid,'tavg1',vid)
      call ncvpt_doubl0(ncid,vid,(1),ntavga,(tavg1),istatus)

      istatus= NF_INQ_VARID(ncid,'tavg2',vid)
      call ncvpt_doubl0(ncid,vid,(1),ntavga,(tavg2),istatus)

      istatus= NF_INQ_VARID(ncid,'ndeltarho',vid)
      call ncvptc0(ncid,vid,1,8,ndeltarho,8,istatus)

      if ((ndeltarho.ne.'disabled'.or.lossmode(1).eq.'simplban').and. &
           ndelta_op.eq."enabled") then
!BH         vid=ncvid(ncid,'deltarho',istatus)
         istatus= NF_INQ_VARID(ncid,'deltarho',vid)
!$$$         do ll=1,lrzmax
!$$$            do l=1,lz
!$$$               do i=1,iy
!$$$                  temp1(i,l)=deltarho(i,l,ll)
!$$$               enddo
!$$$            enddo
!$$$            call pack21(temp1,0,iyp1a,0,jxp1a,tem2,iy,lz)
!$$$            start1(3)=ll
!$$$            call ncvpt(ncid,vid,start1,z_count1,tem2,istatus)
!$$$         enddo
!  Don't need above temp1 flail, since deltarho dynamically
!  dimensioned to size deltarho(iy,lz,lrzmax)
!BH         call ncvpt(ncid,vid,start,z_count,deltarho,istatus)
         call ncvpt_doubl2(ncid,vid,start,z_count,deltarho,istatus)


!deltarhop         WRITE(*,*)'netcdfrw2:deltap_start,deltap_count=',
!deltarhop     +       deltap_start(1:3),deltap_count(1:3)
!deltarhop         vid=ncvid(ncid,'deltarhop',istatus)
!deltarhop         call ncvpt(ncid,vid,deltap_start,deltap_count,deltarhop,
!deltarhop     +        istatus)

!BH         vid=ncvid(ncid,'r_delta',istatus)
!BH         call ncvpt(ncid,vid,1,nr_delta,r_delta(1),istatus)
         istatus= NF_INQ_VARID(ncid,'r_delta',vid)
         call ncvpt_doubl0(ncid,vid,1,nr_delta,r_delta(1:nr_delta),istatus)

!BH         vid=ncvid(ncid,'z_delta',istatus)
!BH         call ncvpt(ncid,vid,1,nz_delta,z_delta(1),istatus)
         istatus= NF_INQ_VARID(ncid,'z_delta',vid)
         call ncvpt_doubl0(ncid,vid,1,nz_delta,z_delta(1:nz_delta),istatus)

!BH         vid=ncvid(ncid,'t_delta',istatus)
!BH         call ncvpt(ncid,vid,1,nt_delta,t_delta(1),istatus)
         istatus= NF_INQ_VARID(ncid,'t_delta',vid)
         call ncvpt_doubl0(ncid,vid,1,nt_delta,t_delta(1:nt_delta),istatus)

!BH         vid=ncvid(ncid,'deltarz',istatus)
!BH         call ncvpt_doubl2(ncid,vid,delta_start,delta_count,deltarz,istatus)
         istatus= NF_INQ_VARID(ncid,'deltarz',vid)
         call ncvpt_doubl2(ncid,vid,delta_start,delta_count,deltarz, &
              istatus)
        WRITE(*,*)'netcdfrw2:start,z_count1=', &
             start(1:3),z_count1(1:3)
        WRITE(*,*)'netcdfrw2:delta_start,delta_count=', &
             delta_start(1:3),delta_count(1:3)

!BH         vid=ncvid_doubl2(ncid,'delta_bdb0',istatus)
!BH         call ncvpt(ncid,vid,delta_start,delta_count,delta_bdb0,
!BH     +        istatus)
         istatus= NF_INQ_VARID(ncid,'delta_bdb0',vid)
         call ncvpt_doubl2(ncid,vid,delta_start,delta_count,delta_bdb0, &
              istatus)

      endif  ! On ndeltarho

      istatus= NF_INQ_VARID(ncid,'bthr',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,(1),lrzmax,bthr(1),istatus)

      istatus= NF_INQ_VARID(ncid,'btoru',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,(1),lrzmax,btoru(1),istatus)

      istatus= NF_INQ_VARID(ncid,'btor0',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,(1),lrzmax,btor0(1),istatus)

      istatus= NF_INQ_VARID(ncid,'bmidplne',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,(1),lrzmax,bmidplne(1),istatus)

      istatus= NF_INQ_VARID(ncid,'softxry',vid)  !-YuP: NetCDF-f77 get vid
      call ncvptc0(ncid,vid,1,8,softxry,8,istatus)

      if (softxry .ne. "disabled") then

         if (x_sxr(1).ne.zero  .or. z_sxr(1).ne.zero) then

            istatus= NF_INQ_VARID(ncid,'x_sxr',vid)  !-YuP: NetCDF-f77 get vid
            call ncvpt_doubl0(ncid,vid,(1),nv,x_sxr(1),istatus)

            istatus= NF_INQ_VARID(ncid,'z_sxr',vid)  !-YuP: NetCDF-f77 get vid
            call ncvpt_doubl0(ncid,vid,(1),nv,z_sxr(1),istatus)

         else

            istatus= NF_INQ_VARID(ncid,'rd',vid)  !-YuP: NetCDF-f77 get vid
            call ncvpt_doubl0(ncid,vid,(1),nv,rd(1),istatus)

            istatus= NF_INQ_VARID(ncid,'thetd',vid)  !-YuP: NetCDF-f77 get vid
            call ncvpt_doubl0(ncid,vid,(1),nv,thetd(1),istatus)

         endif

         istatus= NF_INQ_VARID(ncid,'nv',vid)  !-YuP: NetCDF-f77 get vid
         call ncvpt_int0(ncid,vid,1,1,(nv),istatus)

         istatus= NF_INQ_VARID(ncid,'nen',vid)  !-YuP: NetCDF-f77 get vid
         call ncvpt_int0(ncid,vid,1,1,(nen),istatus)

         istatus= NF_INQ_VARID(ncid,'msxr',vid)  !-YuP: NetCDF-f77 get vid
         call ncvpt_int0(ncid,vid,1,1,(msxr),istatus)

         istatus= NF_INQ_VARID(ncid,'enmin',vid)  !-YuP: NetCDF-f77 get vid
         call ncvpt_doubl0(ncid,vid,(1),1,(enmin),istatus)

         istatus= NF_INQ_VARID(ncid,'enmax',vid)  !-YuP: NetCDF-f77 get vid
         call ncvpt_doubl0(ncid,vid,(1),1,(enmax),istatus)

         istatus= NF_INQ_VARID(ncid,'en_',vid)  !-YuP: NetCDF-f77 get vid
         call ncvpt_doubl0(ncid,vid,(1),nen,en_(1),istatus)

         istatus= NF_INQ_VARID(ncid,'eflux',vid)  !-YuP: NetCDF-f77 get vid
         call pack21(eflux,1,nena,1,nva,wkpack,nen,nv)
         call ncvpt_doubl2(ncid,vid,start_xr,count_xr,wkpack,istatus)

         istatus= NF_INQ_VARID(ncid,'efluxt',vid)  !-YuP: NetCDF-f77 get vid
         call ncvpt_doubl2(ncid,vid,start_xr(2),count_xr(2),efluxt(1), &
              istatus)

      endif  ! On softxry .ne. "disabled"

      istatus= NF_INQ_VARID(ncid,'npa_diag',vid)
      call ncvptc0(ncid,vid,1,8,npa_diag,8,istatus)

      if (npa_diag .ne. "disabled") then

         if (x_npa(1).ne.zero  .or. z_npa(1).ne.zero) then

            istatus= NF_INQ_VARID(ncid,'x_npa',vid)
            call ncvpt_doubl0(ncid,vid,(1),nv_npa,x_npa(1),istatus)

            istatus= NF_INQ_VARID(ncid,'z_npa',vid)
            call ncvpt_doubl0(ncid,vid,(1),nv_npa,z_npa(1),istatus)

         else

            istatus= NF_INQ_VARID(ncid,'rd_npa',vid)
            call ncvpt_doubl0(ncid,vid,(1),nv_npa,rd_npa(1),istatus)

            istatus= NF_INQ_VARID(ncid,'thetd_npa',vid)
            call ncvpt_doubl0(ncid,vid,(1),nv_npa,thetd_npa(1),istatus)

         endif

         istatus= NF_INQ_VARID(ncid,'nv_npa',vid)
         call ncvpt_int0(ncid,vid,1,1,(nv_npa),istatus)

         istatus= NF_INQ_VARID(ncid,'nen_npa',vid)
         call ncvpt_int0(ncid,vid,1,1,(nen_npa),istatus)

         istatus= NF_INQ_VARID(ncid,'npaproc',vid)
         call ncvpt_int0(ncid,vid,1,1,(npaproc),istatus)

         istatus= NF_INQ_VARID(ncid,'enmin_npa',vid)
         call ncvpt_doubl0(ncid,vid,(1),1,(enmin_npa),istatus)

         istatus= NF_INQ_VARID(ncid,'enmax_npa',vid)
         call ncvpt_doubl0(ncid,vid,(1),1,(enmax_npa),istatus)

         istatus= NF_INQ_VARID(ncid,'npa_process',vid)
!         do ii=1,npaproc
!            npa_proc(ii)=npa_process(ii)
!         enddo
         !call ncvptc2(ncid,vid,start_npaproc(1:2),count_npaproc(1:2), &
         !            npa_process,8*npaproc,istatus)
         !istatus = NF_PUT_VAR_TEXT(ncid, vid, npa_process) !YuP: does not work either
         !Note: start_npaproc={1,1}
         !count_npaproc={8,npaproc}
         !npa_process(1:npaproca) is character(len=8)

         istatus= NF_INQ_VARID(ncid,'atten_npa',vid)
         call ncvptc0(ncid,vid,1,8,atten_npa,8,istatus)

         istatus= NF_INQ_VARID(ncid,'ipronn',vid)
         call ncvptc0(ncid,vid,1,8,ipronn,8,istatus)

         istatus= NF_INQ_VARID(ncid,'ennscal',vid)
         call ncvpt_doubl0(ncid,vid,(1),npaproc,ennscal(1),istatus)

         istatus= NF_INQ_VARID(ncid,'en_',vid)
         call ncvpt_doubl0(ncid,vid,(1),nen_npa,en_(1),istatus)

         istatus= NF_INQ_VARID(ncid,'enn',vid)
         call pack21(enn,1,lrza,1,npaproca,wkpack,lrzmax,npaproc)
         call ncvpt_doubl2(ncid,vid,start_npaenn,count_npaenn,wkpack, &
              istatus)

         istatus= NF_INQ_VARID(ncid,'eflux_npa',vid)
         call pack21(eflux,1,nena,1,nva,wkpack,nen_npa,nv_npa)
         call ncvpt_doubl2(ncid,vid,start_npa,count_npa,wkpack,istatus)

         istatus= NF_INQ_VARID(ncid,'efluxt',vid)
         call ncvpt_doubl2(ncid,vid,start_npa(2),count_npa(2),efluxt(1), &
              istatus)

      endif  ! On npa_diag .ne. "disabled"


      istatus= NF_INQ_VARID(ncid,'sigmamod',vid)  !-YuP: NetCDF-f77 get vid
      call ncvptc0(ncid,vid,1,8,sigmamod,8,istatus)
      call check_err(istatus)

      call check_err(istatus)

      if (sigmamod .eq. "enabled") then ! n=0: save values

         istatus= NF_INQ_VARID(ncid,'isigmas',vid)  !-YuP: NetCDF-f77 get vid
         call ncvpt_int0(ncid,vid,1,4,isigmas,istatus)

         istatus= NF_INQ_VARID(ncid,'isigsgv1',vid)  !-YuP: NetCDF-f77 get vid
         call ncvpt_int0(ncid,vid,1,1,(isigsgv1),istatus)

         istatus= NF_INQ_VARID(ncid,'isigsgv2',vid)  !-YuP: NetCDF-f77 get vid
         call ncvpt_int0(ncid,vid,1,1,(isigsgv2),istatus)

         istatus= NF_INQ_VARID(ncid,'mmsv',vid)  !-YuP: NetCDF-f77 get vid
         call ncvpt_int0(ncid,vid,1,1,(mmsv),istatus)

      endif  !  On sigmamod



!.......................................................................
!     Time-Dependent data (numrec1=1)
!     (Continue for numrec1.gt.1 with additional write below.)
!.......................................................................
!     Time-Dependent data (numrec1=1)

      istatus= NF_INQ_VARID(ncid,'time',vid)
      call ncvpt_doubl0(ncid,vid,start(4),1,(timet),istatus)

      !YuP[2018-09-28] added for 'lngshrtf' option,
      !for saving f() distr.func. at selected t steps only.
      if((netcdfshort.eq.'lngshrtf').and.isave.ne.0) then  !isave/nsavet set in tdchief
        !YuP[2019-06-08] Added (netcdfshort.eq.'lngshrtf') in if()
         istatus= NF_INQ_VARID(ncid,'nsave',vid) !here n=0
         call ncvpt_int0(ncid,vid,isave,1,nsave,istatus)
         istatus= NF_INQ_VARID(ncid,'tsave',vid) !here n=0
         call ncvpt_doubl0(ncid,vid,startg(5),1,(timet),istatus) !here startg(5)=1
      endif

      do ll=0,lrzmax
         tr(ll)=reden(kelec,ll)
      enddo
      istatus= NF_INQ_VARID(ncid,'den_e',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start(3),r00_count,tr,istatus)

      call pack21(reden,1,ntotala,0,lrza,wkpack,ntotal,lrzmax)
      istatus= NF_INQ_VARID(ncid,'density',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start(2),species_count,wkpack,istatus)

      istatus= NF_INQ_VARID(ncid,'zeff',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start(3),r0_count(1),zeff,istatus)

      call pack21(temp,1,ntotala,0,lrza,wkpack,ntotal,lrzmax)
      istatus= NF_INQ_VARID(ncid,'temp',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start(2),species_count,wkpack,istatus)

      call pack21(energy,1,ntotala,1,lrza,wkpack,ntotal,lrzmax)
      istatus= NF_INQ_VARID(ncid,'energy',vid) !<..>_FSA
      call ncvpt_doubl2(ncid,vid,start(2),species_count,wkpack,istatus)

      if (ngen.eq.1) then
      call pack21(wpar,1,ngena,1,lrza,wkpack,ngen,lrz)
      istatus= NF_INQ_VARID(ncid,'wpar',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start(3),r_count(1),wkpack,istatus)

      call pack21(wperp,1,ngena,1,lrza,wkpack,ngen,lrz)
      istatus= NF_INQ_VARID(ncid,'wperp',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start(3),r_count(1),wkpack,istatus)

      else  !  ngen.ge.2

      do ll=1,lrz
         do k=1,ngen
            tem1(ll+(k-1)*lrz)=wpar(k,lrindx(ll))
         enddo
      enddo
      istatus= NF_INQ_VARID(ncid,'wpar',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start_rk,count_rk,tem1,istatus)

      do ll=1,lrz
         do k=1,ngen
            tem1(ll+(k-1)*lrz)=wperp(k,lrindx(ll))
         enddo
      enddo
      istatus= NF_INQ_VARID(ncid,'wperp',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start_rk,count_rk,tem1,istatus)

      endif !  on ngen

      istatus= NF_INQ_VARID(ncid,'elecfld',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start(3),r00_count,elecfld,istatus)

      istatus= NF_INQ_VARID(ncid,'edreicer',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start(3),r0_count(1),elecr,istatus)

      call bcast(tr(1:lrzmax),zero,lrzmax)
      do ll=1,lrz
         tr(ll)=vfluxz(lrindx(ll))
      enddo
      istatus= NF_INQ_VARID(ncid,'runaway_rate',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start(3),r_count(1),tr(1:r_count(1)),istatus)

      call bcast(tr(1:lrzmax),zero,lrzmax)
      do ll=1,lrz
         tr(ll)=denra(1,ll) !YuP[2018-09-24]
      enddo
      istatus= NF_INQ_VARID(ncid,'denra',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start(3),r_count(1),tr(1:r_count(1)),istatus)

      call bcast(tr(1:lrzmax),zero,lrzmax)
      do ll=1,lrz
         tr(ll)=curra(1,ll)/3.e9  !Scaling from statA/cm**2 ==> A/cm**2
      enddo
      istatus= NF_INQ_VARID(ncid,'curra',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start(3),r_count(1),tr(1:r_count(1)),istatus)

      tr=zero !YuP[2019-06-08]was call bcast(tr(1),zero,lrzmax)
      do ll=1,lrz
         tr(ll)=ucrit(1,ll) !YuP[2018-09-24]
      enddo
      istatus= NF_INQ_VARID(ncid,'ucrit',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start(3),r_count(1),tr(1:r_count(1)),istatus)

      istatus= NF_INQ_VARID(ncid,'knockon',vid)  !-YuP: NetCDF-f77 get vid
      call ncvptc0(ncid,vid,1,8,knockon,8,istatus)

      if (knockon.ne."disabled") then

      call bcast(tr(1:lrzmax),zero,lrzmax)
      do ll=1,lrz
         tr(ll)=eoe0(1,ll) !YuP[2018-09-24]
      enddo
         istatus= NF_INQ_VARID(ncid,'eoe0',vid)  !-YuP: NetCDF-f77 get vid
         call ncvpt_doubl2(ncid,vid,start(3),r_count(1),tr(1:r_count(1)),istatus)

         istatus= NF_INQ_VARID(ncid,'srckotot',vid)  !-YuP: NetCDF-f77 get vid
         call ncvpt_doubl2(ncid,vid,start(3),r_count(1),srckotot,istatus)

         istatus= NF_INQ_VARID(ncid,'denfl',vid)  !-YuP: NetCDF-f77 get vid
         call ncvpt_doubl2(ncid,vid,start(3),r_count(1),denfl,istatus)

      endif                     ! on knockon
!
!...new Freya stuff
!
      if (frmodp.eq."enabled") then
           istatus= NF_INQ_VARID(ncid,'hibrz',vid)
           call ncvpt_doubl2(ncid,vid,start_hibr,count_hibr, &
                       hibrzp,istatus)

           istatus= NF_INQ_VARID(ncid,'sorpw_nbi',vid)
           call ncvpt_doubl2(ncid,vid,start_sorpw,count_sorpw, &
                       sorpw_nbi,istatus)

           istatus= NF_INQ_VARID(ncid,'sorpw_nbii',vid)
           call ncvpt_doubl2(ncid,vid,start_sorpw,count_sorpw, &
                       sorpw_nbii,istatus)

      endif
!


      if (rdcmod.ne."disabled") then
         if ((mrfn+3)*lrz .gt. iyjx2) stop &
           'netcdfrw2:  Need  (mrfn+3)*lrz<(iy+2)*(jx+2)'
         call bcast(tem1,zero,(mrfn+3)*lrz)
         do ll=1,lrz
            do kk=1,mrfn
               kkk=kk
               tem1(ll+(kk-1)*lrz)=sorpw_rf(kk,lrindx(ll))
            enddo
            tem1(ll+kkk*lrz)=sorpw_rfi(1,lrindx(ll))
            tem1(ll+(kkk+1)*lrz)=sorpwt(lrindx(ll))
            tem1(ll+(kkk+2)*lrz)=sorpwti(lrindx(ll))
         enddo
         istatus= NF_INQ_VARID(ncid,'rfpwr',vid) !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start_rfpwr,count_rfpwr,tem1,istatus)

      endif  !On rdcmod


      if (urfmod.ne."disabled") then
      if ((mrfn+3)*lrz .gt. iyjx2) stop &
           'netcdfrw2:  Need  (mrfn+3)*lrz<(iy+2)*(jx+2)'
      istatus= NF_INQ_VARID(ncid,'mrfn',vid)  !-YuP[2017]added
      call ncvpt_int0(ncid,vid,1,1,(mrfn),istatus)
      !number of rf modes (sum over all wave types and all nharms)
      call bcast(tem1,zero,(mrfn+3)*lrz)
!BH120223:  Removed erroneous divide by dvol(ll) from powrf/powrft
      do ll=1,lrz
         do kk=1,mrfn
            kkk=kk
            tem1(ll+(kk-1)*lrz)=powrf(lrindx(ll),kk)
         enddo
         tem1(ll+kkk*lrz)=powrft(lrindx(ll))
         tem1(ll+(kkk+1)*lrz)=sorpwt(lrindx(ll))
         tem1(ll+(kkk+2)*lrz)=sorpwti(lrindx(ll))
      enddo
      istatus= NF_INQ_VARID(ncid,'rfpwr',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start_rfpwr,count_rfpwr,tem1,istatus)

      call bcast(tem1,zero,mrfn+1)
      do kk=1,mrfn
         tem1(kk)=powurf(kk)
      enddo
      tem1(mrfn+1)=powurf(0)
      istatus= NF_INQ_VARID(ncid,'powurf',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start_powurf,count_powurf,tem1,istatus)

!BH120202:  Should only store mrfn*lrz, at most: NEEDS ADJUSTMENT

      call bcast(tem1,zero,nmodsa*lrz)
      do ll=1,lrz
         do kk=1,nmodsa
            tem1(ll+(kk-1)*lrz)=powrfl(lrindx(ll),kk)
         enddo
      enddo

      istatus= NF_INQ_VARID(ncid,'powrfl',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start_powrf,count_powrf,tem1,istatus)

      istatus= NF_INQ_VARID(ncid,'powurfl',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start_powrf(2),count_powrf(2), &
                 powurfl(1),istatus)


!BH120202:  Should only store mrfn*lrz, at most: NEEDS ADJUSTMENT
      call bcast(tem1,zero,nmodsa*lrz)
      do ll=1,lrz
         do kk=1,nmodsa
            tem1(ll+(kk-1)*lrz)=powrf(lrindx(ll),kk)
         enddo
      enddo
      istatus= NF_INQ_VARID(ncid,'powrf',vid)
      call ncvpt_doubl2(ncid,vid,start_powrf,count_powrf,tem1,istatus)

!BH120202:  Should only store mrfn*lrz, at most: NEEDS ADJUSTMENT
      call bcast(tem1,zero,nmodsa*lrz)
      do ll=1,lrz
         do kk=1,nmodsa
            tem1(ll+(kk-1)*lrz)=powrfc(lrindx(ll),kk)
         enddo
      enddo
      istatus= NF_INQ_VARID(ncid,'powrfc',vid)
      call ncvpt_doubl2(ncid,vid,start_powrf,count_powrf,tem1,istatus)

      istatus= NF_INQ_VARID(ncid,'powurfc',vid)
      call ncvpt_doubl2(ncid,vid,start_powrf(2),count_powrf(2), &
                 powurfc(1),istatus)

      call bcast(tem1,zero,lrz)
      do ll=1,lrz
            tem1(ll)=powrft(lrindx(ll))
      enddo
      istatus= NF_INQ_VARID(ncid,'powrft',vid)
      call ncvpt_doubl2(ncid,vid,start_powrft,count_powrft,tem1,istatus)

      istatus= NF_INQ_VARID(ncid,'nrfspecies',vid)  !YuP[11-2017]
      call ncvpt_int0(ncid,vid,1,count_powrf(2),nrfspecies,istatus)
      endif  !On urfmod

      istatus= NF_INQ_VARID(ncid,'curtor',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start(3),r0_count(1),curtor,istatus)

      istatus= NF_INQ_VARID(ncid,'ccurtor',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start(3),r0_count(1),ccurtor(1),istatus)

      istatus= NF_INQ_VARID(ncid,'curpol',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start(3),r0_count(1),curpol,istatus)

      istatus= NF_INQ_VARID(ncid,'ccurpol',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start(3),r0_count(1),ccurpol(1),istatus)

      if (ngen.eq.1) then

         k=1
         do ll=1,lrzmax
            tr(ll)=curr(k,ll)/3.e9  !Scaling statA/cm**2 ==> A/cm**2
         enddo
         istatus= NF_INQ_VARID(ncid,'curr',vid)
         call ncvpt_doubl2(ncid,vid,start(3),r0_count(1),tr(1),istatus)
         !YuP[07-31-2014] Added:
         do ll=1,lrzmax
            tr(ll)=energym(k,ll)
         enddo
         istatus= NF_INQ_VARID(ncid,'energym',vid)
         call ncvpt_doubl2(ncid,vid,start(3),r0_count(1),tr(1),istatus)

      else  !  ngen.ge.2
         do k=1,ngen
         do ll=1,lrzmax
               tem1(ll+(k-1)*lrzmax)=curr(k,ll)/3.e9
            ! 1/3e9 is for statA/cm**2 ==> A/cm**2
         enddo
         enddo
         istatus= NF_INQ_VARID(ncid,'curr',vid)
         call ncvpt_doubl2(ncid,vid,start_r0k,count_r0k,tem1,istatus)
         !YuP[07-2017] added
         do k=1,ngen
         do ll=1,lrzmax
            tem1(ll+(k-1)*lrzmax)=energym(k,ll) !YuP[07-31-2014] Added
         enddo
         enddo
         istatus= NF_INQ_VARID(ncid,'energym',vid)
         call ncvpt_doubl2(ncid,vid,start_r0k,count_r0k,tem1,istatus)

      endif ! ngen>1


      istatus= NF_INQ_VARID(ncid,'efflag',vid)  !-YuP: NetCDF-f77 get vid
      call ncvptc0(ncid,vid,1,8,efflag,8,istatus)

      kk=1
      if (kelecg.ne.0) kk=kelecg
      do ll=1,lrzmax
         tr(ll)=currm(kk,ll)/3.e9  !Scaling statA/cm**2 ==> A/cm**2
      enddo
      if (kelecg.ne.0) then ! e as gen.species
         istatus= NF_INQ_VARID(ncid,'currm_e',vid)
      else ! ion as general species
         istatus= NF_INQ_VARID(ncid,'currm_i',vid)
      endif
      call ncvpt_doubl2(ncid,vid,start(3),r0_count(1),tr(1),istatus)

      do ll=1,lrzmax
         tr(ll)=restp(nch(ll),ll)
      enddo
      istatus= NF_INQ_VARID(ncid,'restp',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start(3),r0_count(1),tr(1),istatus)

      do ll=1,lrzmax
         tr(ll)=restnp(nch(ll),ll)
      enddo
      istatus= NF_INQ_VARID(ncid,'restnp',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start(3),r0_count(1),tr(1),istatus)

      do ll=1,lrzmax
         tr(ll)=sptzrp(nch(ll),ll)
      enddo
      istatus= NF_INQ_VARID(ncid,'sptzrp',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start(3),r0_count(1),tr(1),istatus)

      istatus= NF_INQ_VARID(ncid,'rovsc',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start(3),r0_count(1),rovsc,istatus)

      istatus= NF_INQ_VARID(ncid,'rovsc_hi',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start(3),r0_count(1),rovsc,istatus)

      istatus= NF_INQ_VARID(ncid,'zreskim',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start(3),r0_count(1),zreskim,istatus)

      istatus= NF_INQ_VARID(ncid,'taueeh',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start(3),r0_count(1),taueeh,istatus)

      istatus= NF_INQ_VARID(ncid,'nuestar',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start(3),r0_count(1),starnue,istatus)

!BH110320c     Following only set up for ngen=1
!BH110320      if (ngen.ge.2)
!BH110320     1     WRITE(*,*)'netcdfrw2: Tot pwrs only set up for ngen=1'
      if (13*lrz*ngen .gt. iyjx2) stop &
              'netcdfrw2: Need 13*lrz*ngen<(iy+2)*(jx+2)'
      do k=1,ngen
         kkk=(k-1)*13*lrz
         call bcast(tem1,zero,13*lrz*ngen)
         do ll=1,lrz
            kk=1
            tem1(ll+(kk-1)*lrz+kkk)=entr(k,-1,ll)
            kk=2
            tem1(ll+(kk-1)*lrz+kkk)=entr(k,0,ll)
            kk=3
            tem1(ll+(kk-1)*lrz+kkk)=entr(k,2,ll)
            kk=4
            tem1(ll+(kk-1)*lrz+kkk)=entr(k,1,ll)
            kk=5
            tem1(ll+(kk-1)*lrz+kkk)=entr(k,3,ll)
            kk=6
            tem1(ll+(kk-1)*lrz+kkk)=entr(k,5,ll)
            kk=7
            tem1(ll+(kk-1)*lrz+kkk)=entr(k,6,ll)
            kk=8
            tem1(ll+(kk-1)*lrz+kkk)=entr(k,7,ll)
            kk=9
            tem1(ll+(kk-1)*lrz+kkk)=entr(k,8,ll)
            kk=10
            tem1(ll+(kk-1)*lrz+kkk)=entr(k,11,ll)
            kk=11
            tem1(ll+(kk-1)*lrz+kkk)=entr(k,10,ll)
            kk=12
            tem1(ll+(kk-1)*lrz+kkk)=entr(k,12,ll)
            kk=13
            tem1(ll+(kk-1)*lrz+kkk)=entr(k,4,ll)
         enddo
      enddo
      istatus= NF_INQ_VARID(ncid,'powers',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start_powers,count_powers,tem1,istatus)

!     Following only set up for ngen=1
!BH110320      if (ngen.ge.2)
!BH110320     1     WRITE(*,*)'netcdfrw2: Tot pwrs only set up for ngen=1'
      if (13*ngen .gt.iyjx2) stop &
           'netcdfrw2: Need 13*ngen<(iy+2)*(jx+2)'
      do k=1,ngen
            kkk=(k-1)*13
            call bcast(tem1,zero,13*ngen)
            kk=1
            tem1(kk+kkk)=entrintr(k,-1)
            kk=2
            tem1(kk+kkk)=entrintr(k,0)
            kk=3
            tem1(kk+kkk)=entrintr(k,2)
            kk=4
            tem1(kk+kkk)=entrintr(k,1)
            kk=5
            tem1(kk+kkk)=entrintr(k,3)
            kk=6
            tem1(kk+kkk)=entrintr(k,5)
            kk=7
            tem1(kk+kkk)=entrintr(k,6)
            kk=8
            tem1(kk+kkk)=entrintr(k,7)
            kk=9
            tem1(kk+kkk)=entrintr(k,8)
            kk=10
            tem1(kk+kkk)=entrintr(k,11)
            kk=11
            tem1(kk+kkk)=entrintr(k,10)
            kk=12
            tem1(kk+kkk)=entrintr(k,12)
            kk=13
            tem1(kk+kkk)=entrintr(k,4)
      enddo
      istatus= NF_INQ_VARID(ncid,'powers_int',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start_powers(2),count_powers(2), &
                        tem1,istatus)


      if (sigmamod .eq. "enabled") then ! n>0: save sigftt
         istatus= NF_INQ_VARID(ncid,'sigftt',vid)  !-YuP: NetCDF-f77 get vid
         do lsig=1,4
            tem2(lsig)=sigftt(nch(1),lsig)
         enddo
         call ncvpt_doubl2(ncid,vid,start_fus(2),count_fus(2), &
                           tem2,istatus)
      endif

      if ( netcdfshort.eq.'long_jp' ) then


         if (ngen.eq.1) then

            k=1
            do ll=1,lrz
               do j=1,jx
                  i=j+(lrindx(ll)-1)*jx
                  tem1(i)=currv(j,k,lrindx(ll))/3.e9
               enddo
            enddo
            istatus= NF_INQ_VARID(ncid,'currv',vid)  !-YuP: NetCDF-f77 get vid
            call ncvpt_doubl2(ncid,vid,start(2),count(2),tem1,istatus)

            do ll=1,lrz
               do j=1,jx
                  i=j+(lrindx(ll)-1)*jx
                  tem1(i)=pwrrf(j,k,lrindx(ll))
               enddo
            enddo
            istatus= NF_INQ_VARID(ncid,'pwrrf',vid) !-YuP: NetCDF-f77 get vid
            call ncvpt_doubl2(ncid,vid,start(2),count(2),tem1,istatus)

         else  !  ngen.ge.2

         do k=1,ngen
            do ll=1,lrz
               do j=1,jx
                  i=j+(lrindx(ll)-1)*jx
                  tem1(i)=currv(j,k,lrindx(ll))/3.e9
               enddo
            enddo
            startg(4)=k
            countg(3)=lrz
            istatus= NF_INQ_VARID(ncid,'currv',vid)  !-YuP: NetCDF-f77 get vid
            call ncvpt_doubl2(ncid,vid,startg(2),countg(2),tem1,istatus)

            do ll=1,lrz
               do j=1,jx
                  i=j+(lrindx(ll)-1)*jx
                  tem1(i)=pwrrf(j,k,lrindx(ll))
               enddo
            enddo
            startg(4)=k
            countg(3)=lrz
            istatus= NF_INQ_VARID(ncid,'pwrrf',vid) !-YuP: NetCDF-f77 get vid
            call ncvpt_doubl2(ncid,vid,startg(2),countg(2),tem1,istatus)
         enddo  !  on k=1,ngen

         endif  !  on ngen
         startg(4)=1 ! restored
         countg(3)=1 ! restored

      endif  ! on netcdfshort.eq.'long_jp'


!YuP      if ( netcdfshort.eq.'longer_f' ) then  !endif at line 3090
      if ( (netcdfshort.eq.'longer_f').or. &
           (netcdfshort.eq.'lngshrtf')    ) then
         !---> Here n=0. Should we record f() even if isave=0?
         !---> If-yes, do not use ".and.(isave.ne.0)" as in n>0 case.
         !write(*,*)'netcdfrw2: netcdfshort,n,isave=',netcdfshort,n,isave
         istatus= NF_INQ_VARID(ncid,'f',vid)  !-YuP: NetCDF-f77 get vid
         if (ngen.eq.1) then

            if(netcdfshort.eq.'lngshrtf') then
              start1(4)=numrecsave !YuP[2018-09-28]netcdfshort.eq.'lngshrtf'
            else ! 'longer_f'
              start1(4)=numrec1 ! usual (saved at every t step) Here numrec1=1
            endif
            do ll=1,lrz
               do j=1,jx
                  do i=1,iy
                     temp1(i,j)=f(i,j,1,lrindx(ll))
                  enddo
               enddo
               call pack21(temp1,0,iyp1,0,jxp1,wkpack,iy,jx)
               start1(3)=ll
               call ncvpt_doubl2(ncid,vid,start1,count1,wkpack,istatus)
               !So, here start1={1,   1,   ll,   numrec1(or numrecsave)}
               !         count1={iy,  jx,  1,    1}
               !         dimsf= {ydim,xdim,rdim, tdim(or tsavedim)}
            enddo
            start1(4)=numrec1 ! restore (=1 here, anyway)

         else  !ngen.ge.2

            if(netcdfshort.eq.'lngshrtf') then
              startg(5)=numrecsave !YuP[2018-09-28]netcdfshort.eq.'lngshrtf'
            else ! 'longer_f'
              startg(5)=numrec1 ! usual (saved at every t step) here numrec1=1
            endif
            do k=1,ngen
            do ll=1,lrz
               do j=1,jx
                  do i=1,iy
                     temp1(i,j)=f(i,j,k,lrindx(ll))
                  enddo
               enddo
               call pack21(temp1,0,iyp1,0,jxp1,wkpack,iy,jx)
               startg(3)=ll
               startg(4)=k
               call ncvpt_doubl2(ncid,vid,startg,countg,wkpack,istatus)
               !So, here startg={1,   1,   ll,  k,   numrec1(or numrecsave)}
               !         countg={iy,  jx,  1,   1,   1}
               !         dimsg= {ydim,xdim,rdim,gdim,tdim(or tsavedim)}
            enddo  !  On ll
            enddo  !  On k
            startg(3)=1
            startg(4)=1
            startg(5)=numrec1 ! restore (=1 here, anyway)

         endif  ! on ngen


      endif  ! on netcdfshort.eq.'longer_f'




! --- output large arrays only if nstop=0
      if (nstop.eq.0) then

!$$$         istatus= NF_INQ_VARID(ncid,'currv',vid)  !-YuP: NetCDF-f77 get vid
!$$$         k=1
!$$$         do ll=1,lrz
!$$$            do j=1,jx
!$$$               i=j+(lrindx(ll)-1)*jx
!$$$               tem1(i)=currv(j,k,lrindx(ll))/3.e9
!$$$            enddo
!$$$         enddo
!$$$         call ncvpt_doubl2(ncid,vid,start(2),count(2),tem1,istatus)
!$$$
!$$$         istatus= NF_INQ_VARID(ncid,'pwrrf',vid)  !-YuP: NetCDF-f77 get vid
!$$$         k=1
!$$$         do ll=1,lrz
!$$$            do j=1,jx
!$$$               i=j+(lrindx(ll)-1)*jx
!$$$               tem1(i)=pwrrf(j,k,lrindx(ll))
!$$$            enddo
!$$$         enddo
!$$$         call ncvpt_doubl2(ncid,vid,start(2),count(2),tem1,istatus)

         if ( (netcdfshort.ne.'enabled')  .and. &
              (netcdfshort.ne.'longer_f') .and. &
              (netcdfshort.ne.'lngshrtf')       ) then ! here n=0 (nstop=0)

            istatus= NF_INQ_VARID(ncid,'f',vid)  !-YuP: NetCDF-f77 get vid
!BH011221: This storage is set up for constant iy as function of radius.
!BH011221: Needs generalizing. Should we store in reg array to max(iy_)?
!BH011221: For now, simply stop.
            do ll=1,lrz
               if (iy_(ll).ne.iymax) &
                    stop 'netcdfrw2: Cant handle iy.ne.iymax'
            enddo

            if (ngen.eq.1) then

               k=1 !
               if (tavg.ne."disabled") then
                 PRINT *,'netcdfrw2 WARNING: favg is saved (NOT f)'
                 PRINT *,'netcdfrw2 WARN:  but if time(nstop)<tavg1(1)'
                 PRINT *,'netcdfrw2 WARNING: then favg = f.'
               endif
               do ll=1,lrz
                  if (tavg.eq."disabled") then
                     do j=1,jx
                     do i=1,iy
                        temp1(i,j)=f(i,j,k,lrindx(ll))
                        if(gone(i,j,k,lrindx(ll)).lt.-0.1) then
                           temp1(i,j)=em90
                        endif
                     enddo
                     enddo
                  else  !On tavg = enabled
                     ! Bob, here is the nstop=0 part,
                     ! so how can we have favg?
                     do j=1,jx
                     do i=1,iy
                        temp1(i,j)=favg(i,j,k,lrindx(ll))
                     !if(gone(i,j,k,lrindx(ll)).lt.-0.1) temp1(i,j)=em90
                     enddo
                     enddo
                     WRITE(*,'(a,i6,2e13.4)') &
                     'nstop0.netcdfrw2/tavg=en ll,sumij(favg),sumij(f)', &
                       ll,sum(temp1),sum(f(:,:,k,lrindx(ll)))
                  endif  !On tavg
!                 temp1 dimensnd 0:iyp1,0,jxp1. Pack in to (1:iy,1:jx)
                  call pack21(temp1,0,iyp1,0,jxp1,wkpack,iy,jx)
                  start1(3)=ll
                  call ncvpt_doubl2(ncid,vid,start1(1:3),count1(1:3), &
                                    wkpack,istatus)
               !So, here start1={1,   1,   ll}
               !         count1={iy,  jx,  1}
               !         dimsf= {ydim,xdim,rdim}
               enddo  !On ll
               start1(3)=1 !restore

            else  !  ngen.ge.2

               do k=1,ngen
               if (tavg.eq."disabled") then
               do ll=1,lrz
                  do j=1,jx
                     do i=1,iy
                        temp1(i,j)=f(i,j,k,lrindx(ll))
                     enddo
                  enddo
                  call pack21(temp1,0,iyp1,0,jxp1,wkpack,iy,jx)
                  startg(3)=ll
                  startg(4)=k
                 call ncvpt_doubl2(ncid,vid,startg(1:4),countg(1:4), &
                                   wkpack,istatus)
               !So, here startg={1,   1,   ll,  k}
               !         countg={iy,  jx,  1,   1}
               !         dimsg= {ydim,xdim,rdim,gdim}
               enddo  ! on ll
               else  ! On tavg
               do ll=1,lrz
                  do j=1,jx
                     do i=1,iy
                        temp1(i,j)=favg(i,j,k,lrindx(ll))
                     enddo
                  enddo
                  call pack21(temp1,0,iyp1,0,jxp1,wkpack,iy,jx)
                  startg(3)=ll
                  startg(4)=k
                  call ncvpt_doubl2(ncid,vid,startg(1:4),countg(1:4), &
                                    wkpack,istatus)
               !So, here startg={1,   1,   ll,  k}
               !         countg={iy,  jx,  1,   1}
               !         dimsg= {ydim,xdim,rdim,gdim}
               enddo  ! on ll
               endif  !On tavg
               enddo  ! on k
               startg(3)=1
               startg(4)=1

            endif  !  on ngen

         endif  !  on netcdfshort
      endif  ! on nstop.eq.0

! --- endif ((kopt.eq.0) .and. (n.eq.0)) ---
      endif







!$$$!-----------------------------------------------------------------------
!$$$!
!$$$!cl    2. Restart from previous run
!            (THIS COMMENTED OUT SECTION IS PART OF FPET
!             ROUTINES USED AS A TEMPLATE FOR THIS CQL3D
!             NETCDF SUBROUTINE SET.  BobH, 2000)
!$$$!        The following data is NOT read in:
!$$$!          vnorm,x,y,z,dvol,flux,press,hflux,pdens,flux_a,flux_b,
!$$$!          hflux_a,hflux_b,pdens0
!$$$!
!$$$
!$$$! --- begin if ---
!$$$      if ((kopt.ne.0) .and. (n.eq.0)) then
!$$$
!$$$!.......................................................................
!$$$!cl    2.1 Open previous netCDF file
!$$$
!$$$      ncid = ncopn(filename,NCWRITE,istatus)
!$$$
!$$$!.......................................................................
!$$$!cl    2.2 read in dimension IDs and sizes
!$$$
!$$$      xdim = ncdid(ncid,'xdim',istatus)
!$$$      ydim = ncdid(ncid,'ydim',istatus)
!$$$      rdim = ncdid(ncid,'zdim',istatus)
!$$$      tdim = ncdid(ncid,'taudim',istatus)
!$$$
!$$$! --- inquire about dimension sizes ---
!$$$!     ncdinq(netCDF_id, dimension_id_from_ncdid, returned_dim_name,
!$$$!     returned_dim_size)
!$$$!     Note: for unlimited dimension, returned_dim_size=current maximum
!$$$!     which is the same as the maximum record number
!$$$
!$$$      call ncdinq(ncid,ydim,name,iyp,istatus)
!$$$      call ncdinq(ncid,xdim,name,jxp,istatus)
!$$$      call ncdinq(ncid,rdim,name,kzp,istatus)
!$$$
!$$$! --- stop if dimension sizes don't match ---
!$$$      if ((iyp.ne.iy) .or. (jxp.ne.jx) .or. (kzp.ne.kz)) then
!$$$         write(6,*) "set the following parameter values in RUN_FPET"
!$$$         write(6,*) "  iy = ",iyp,"  jx = ",jxp, "  kz = ",kzp
!$$$         stop "non-matching dimensions in NCDFWRITE"
!$$$      endif
!$$$
!$$$! --- set the time-step counter ==> numrec1
!$$$      call ncdinq(ncid,tdim,name,numrec1,istatus)
!$$$      start(4)=numrec1
!$$$
!$$$!.......................................................................
!$$$!cl    2.3 Read datas
!$$$!
!$$$!     Here we read in only what's needed to re-start
!$$$!     the run from the previous time-step.
!$$$!
!$$$!     ncvgt(netCDF_id, variable_id_from_ncvid, vector_of_starting_index,
!$$$!     vector_of_lengths, returned_variable, integer_info)
!$$$!
!$$$      vid = ncvid(ncid,'tau',istatus)
!$$$      call ncvgt(ncid,vid,start(4),1,tau,istatus)
!$$$
!$$$      vid = ncvid(ncid,'dtau',istatus)
!$$$      call ncvgt(ncid,vid,1,1,dtau,istatus)
!$$$
!$$$      vid = ncvid(ncid,'elecfld',istatus)
!$$$      call ncvgt(ncid,vid,start(3),count(3),zarray1d(1),istatus)
!$$$      do k=start(3),start(3)+count(3)-1
!$$$        elecfld(k)=zarray1d(k)
!$$$      enddo
!$$$      .   .   .   .  .   .   .    .    .     .
!$$$
!$$$! --- endif ((kopt.ne.0) .and. (n.eq.0)) ---
!$$$      endif

!-----------------------------------------------------------------------









!
!cl    3. Periodic save at each time-step (numrec1.gt.1)
!

! --- begin if ---
      if (n.gt.0) then  !endif at line 3376

!.......................................................................
!cl    3.1 set the time-step counter ==> numrec1

      ! This is for saving data at EACH t step:
      numrec1=numrec1+1
      start(4)=numrec1
      start1(4)=numrec1
      startg(5)=numrec1
      start_powurf(2)=numrec1
      start_rfpwr(3)=numrec1
      start_powrf(3)=numrec1
      start_powrft(2)=numrec1
      start_powers(4)=numrec1
      start_fus(3)=numrec1
      start_xr(3)=numrec1
      start_npa(3)=numrec1
      start_rk(3)=numrec1
      start_r0k(3)=numrec1

      if((netcdfshort.eq.'lngshrtf').and.(isave.ne.0))then
      !YuP[2018-09-28] Saving distr.func. at selected t steps only.
      !Increment index only for steps specified by nsave() array:
      numrecsave=numrecsave+1 ! for (netcdfshort.eq.'lngshrtf')
      startgsave(4)=numrecsave ! for (netcdfshort.eq.'lngshrtf')
      startg(5)=numrecsave ! for (netcdfshort.eq.'lngshrtf')
      !YuP[2018-09-28] added for 'lngshrtf' option,
      !for saving f() distr.func. at selected t steps only.
      istatus= NF_INQ_VARID(ncid,'tsave',vid) ! here n>0
      call ncvpt_doubl0(ncid,vid,startgsave(4),1,(timet),istatus) !can use startg(5) here
      !Note: for netcdfshort.ne.'lngshrtf' startgsave(4) remains =1
      !and 'tsave' was recorded for n=0 only
      endif

!.......................................................................
!cl    3.2 Variables saved at each time-step (numrec1.gt.1)

      istatus= NF_INQ_VARID(ncid,'time',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,start(4),1,(timet),istatus)


      do ll=0,lrzmax
         tr(ll)=reden(kelec,ll)
      enddo
      istatus= NF_INQ_VARID(ncid,'den_e',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start(3),r00_count(1),tr,istatus)

      call pack21(reden,1,ntotala,0,lrza,wkpack,ntotal,lrzmax)
      istatus= NF_INQ_VARID(ncid,'density',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start(2),species_count,wkpack,istatus)

      istatus= NF_INQ_VARID(ncid,'zeff',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start(3),r0_count(1),zeff,istatus)

      istatus= NF_INQ_VARID(ncid,'temp',vid)  !-YuP: NetCDF-f77 get vid
      call pack21(temp,1,ntotala,0,lrza,wkpack,ntotal,lrzmax)
      call ncvpt_doubl2(ncid,vid,start(2),species_count,wkpack,istatus)

      istatus= NF_INQ_VARID(ncid,'energy',vid)
      call pack21(energy,1,ntotala,1,lrza,wkpack,ntotal,lrzmax)
      call ncvpt_doubl2(ncid,vid,start(2),species_count,wkpack,istatus)


      if (ngen.eq.1) then
      call pack21(wpar,1,ngena,1,lrza,wkpack,ngen,lrz)
      istatus= NF_INQ_VARID(ncid,'wpar',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start(3),r_count(1),wkpack,istatus)

      call pack21(wperp,1,ngena,1,lrza,wkpack,ngen,lrz)
      istatus= NF_INQ_VARID(ncid,'wperp',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start(3),r_count(1),wkpack,istatus)

      else  !  ngen.ge.2

      do ll=1,lrz
         do k=1,ngen
            tem1(ll+(k-1)*lrz)=wpar(k,lrindx(ll))
         enddo
      enddo
      istatus= NF_INQ_VARID(ncid,'wpar',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start_rk,count_rk,tem1,istatus)

      do ll=1,lrz
         do k=1,ngen
            tem1(ll+(k-1)*lrz)=wperp(k,lrindx(ll))
         enddo
      enddo
      istatus= NF_INQ_VARID(ncid,'wperp',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start_rk,count_rk,tem1,istatus)

      endif ! ngen

      istatus= NF_INQ_VARID(ncid,'elecfld',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,start(3),r00_count(1),elecfld,istatus)

      istatus= NF_INQ_VARID(ncid,'edreicer',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,start(3),r0_count(1),elecr,istatus)

      call bcast(tr(1:lrzmax),zero,lrzmax)
      do ll=1,lrz
         tr(ll)=vfluxz(lrindx(ll))
      enddo
      istatus= NF_INQ_VARID(ncid,'runaway_rate',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,start(3),r_count(1),tr(1:r_count(1)),istatus)

      call bcast(tr(1:lrzmax),zero,lrzmax)
      do ll=1,lrz
         tr(ll)=denra(1,ll) !YuP[2018-09-24]
      enddo
      istatus= NF_INQ_VARID(ncid,'denra',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,start(3),r_count(1),tr(1:r_count(1)),istatus)

      call bcast(tr(1:lrzmax),zero,lrzmax)
      do ll=1,lrz
         tr(ll)=curra(1,ll)/3.e9
      enddo
      istatus= NF_INQ_VARID(ncid,'curra',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,start(3),r_count(1),tr(1:r_count(1)),istatus)

      call bcast(tr(1:lrzmax),zero,lrzmax)
      do ll=1,lrz
         tr(ll)=ucrit(1,ll) !YuP[2018-09-24]
      enddo
      istatus= NF_INQ_VARID(ncid,'ucrit',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,start(3),r_count(1),tr(1:r_count(1)),istatus)

      if (knockon.ne."disabled") then

      call bcast(tr(1:lrzmax),zero,lrzmax)
      do ll=1,lrz
         tr(ll)=eoe0(1,ll) !YuP[2018-09-24]
      enddo
      istatus= NF_INQ_VARID(ncid,'eoe0',vid)  !-YuP: NetCDF-f77 get vid
         call ncvpt_doubl0(ncid,vid,start(3),r_count(1),tr(1:r_count(1)),istatus)

      istatus= NF_INQ_VARID(ncid,'srckotot',vid)  !-YuP: NetCDF-f77 get vid
         call ncvpt_doubl0(ncid,vid,start(3),r_count(1),srckotot,istatus)

      istatus= NF_INQ_VARID(ncid,'denfl',vid)  !-YuP: NetCDF-f77 get vid
         call ncvpt_doubl0(ncid,vid,start(3),r_count(1),denfl,istatus)

      endif                     ! on knockon

!
!...new Freya stuff
!
      if (frmodp.eq."enabled") then
         istatus= NF_INQ_VARID(ncid,'hibrz',vid)
         call ncvpt_doubl2(ncid,vid,start_hibr,count_hibr, &
                       hibrzp,istatus)

         istatus= NF_INQ_VARID(ncid,'sorpw_nbi',vid)
         call ncvpt_doubl2(ncid,vid,start_sorpw,count_sorpw, &
                       sorpw_nbi,istatus)

         istatus= NF_INQ_VARID(ncid,'sorpw_nbii',vid)
         call ncvpt_doubl2(ncid,vid,start_sorpw,count_sorpw, &
                           sorpw_nbii,istatus)
      endif
!

      if (rdcmod.ne."disabled") then

         if ((mrfn+3)*lrz .gt. iyjx2) stop &
              'netcdfrw2:  Need  (mrfn+3)*lrz<(iy+2)*(jx+2)'
         call bcast(tem1,zero,(mrfn+3)*lrz)
!BH120223:  Removed erroneous divide by dvol(ll) from powrf/powrft
         do ll=1,lrz
            do kk=1,mrfn
               kkk=kk
               tem1(ll+(kk-1)*lrz)=sorpw_rf(kk,lrindx(ll))
            enddo
            tem1(ll+kkk*lrz)=sorpw_rfi(1,lrindx(ll)) !For 1 gen species
            tem1(ll+(kkk+1)*lrz)=sorpwt(lrindx(ll))
            tem1(ll+(kkk+2)*lrz)=sorpwti(lrindx(ll))
         enddo
!     -YuP:      vid=ncvid(ncid,'rfpwr',istatus)
        istatus= NF_INQ_VARID(ncid,'rfpwr',vid) !-YuP: NetCDF-f77 get vid
        call ncvpt_doubl2(ncid,vid,start_rfpwr,count_rfpwr,tem1,istatus)

      endif  !On rdcmod


      if (urfmod.ne."disabled") then

      if ((mrfn+3)*lrz .gt. iyjx2) stop &
           'netcdfrw2:  Need  (mrfn+3)*lrz<(iy+2)*(jx+2)'
      call bcast(tem1,zero,(mrfn+3)*lrz)
!BH120223:  Removed erroneous divide by dvol(ll) from powrf/powrft
      do ll=1,lrz
         do kk=1,mrfn
            kkk=kk
            tem1(ll+(kk-1)*lrz)=powrf(lrindx(ll),kk)
         enddo
         tem1(ll+kkk*lrz)=powrft(lrindx(ll))
         tem1(ll+(kkk+1)*lrz)=sorpwt(lrindx(ll))
         tem1(ll+(kkk+2)*lrz)=sorpwti(lrindx(ll))
      enddo
      istatus= NF_INQ_VARID(ncid,'rfpwr',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start_rfpwr,count_rfpwr,tem1,istatus)

      call bcast(tem1,zero,mrfn+1)
      do kk=1,mrfn
         tem1(kk)=powurf(kk)
      enddo
      tem1(mrfn+1)=powurf(0)
      istatus= NF_INQ_VARID(ncid,'powurf',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start_powurf,count_powurf,tem1,istatus)

!BH120202:  Should only store mrfn*lrz, at most: NEEDS ADJUSTMENT
      call bcast(tem1,zero,nmodsa*lrz)
      do ll=1,lrz
         do kk=1,nmodsa
            tem1(ll+(kk-1)*lrz)=powrfl(lrindx(ll),kk)
         enddo
      enddo
      istatus= NF_INQ_VARID(ncid,'powrfl',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start_powrf,count_powrf,tem1,istatus)

      istatus= NF_INQ_VARID(ncid,'powurfl',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start_powrf(2),count_powrf(2), &
                 powurfl(1),istatus)

!BH120202:  Should only store mrfn*lrz, at most: NEEDS ADJUSTMENT
      call bcast(tem1,zero,nmodsa*lrz)
      do ll=1,lrz
         do kk=1,nmodsa
            tem1(ll+(kk-1)*lrz)=powrf(lrindx(ll),kk)
         enddo
      enddo
      istatus= NF_INQ_VARID(ncid,'powrf',vid)
      call ncvpt_doubl2(ncid,vid,start_powrf,count_powrf,tem1,istatus)


!BH120202:  Should only store mrfn*lrz, at most: NEEDS ADJUSTMENT
      call bcast(tem1,zero,nmodsa*lrz)
      do ll=1,lrz
         do kk=1,nmodsa
            tem1(ll+(kk-1)*lrz)=powrfc(lrindx(ll),kk)
         enddo
      enddo

      istatus= NF_INQ_VARID(ncid,'powrfc',vid)
      call ncvpt_doubl2(ncid,vid,start_powrf,count_powrf,tem1,istatus)

      istatus= NF_INQ_VARID(ncid,'powurfc',vid)
      call ncvpt_doubl2(ncid,vid,start_powrf(2),count_powrf(2), &
                 powurfc(1),istatus)

      call bcast(tem1,zero,lrz)
      do ll=1,lrz
            tem1(ll)=powrft(lrindx(ll))
      enddo
      istatus= NF_INQ_VARID(ncid,'powrft',vid)
      call ncvpt_doubl2(ncid,vid,start_powrft,count_powrft,tem1,istatus)

      endif ! urfmod.ne."disabled"


      istatus= NF_INQ_VARID(ncid,'curtor',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,start(3),r0_count(1),curtor,istatus)

      istatus= NF_INQ_VARID(ncid,'ccurtor',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,start(3),r0_count(1),ccurtor(1),istatus)

      istatus= NF_INQ_VARID(ncid,'curpol',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,start(3),r0_count(1),curpol,istatus)

      istatus= NF_INQ_VARID(ncid,'ccurpol',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,start(3),r0_count(1),ccurpol(1),istatus)

      if (ngen.eq.1) then
         k=1
         do ll=1,lrzmax
            tr(ll)=curr(k,ll)/3.e9
         enddo
         istatus= NF_INQ_VARID(ncid,'curr',vid)
         call ncvpt_doubl0(ncid,vid,start(3),r0_count(1),tr(1:),istatus)
         !YuP[07-31-2014] Added:
         do ll=1,lrzmax
            tr(ll)=energym(k,ll)
         enddo
         istatus= NF_INQ_VARID(ncid,'energym',vid)
         call ncvpt_doubl0(ncid,vid,start(3),r0_count(1),tr(1:),istatus)
      else  !  ngen.ge.2
         do k=1,ngen
         do ll=1,lrzmax
               tem1(ll+(k-1)*lrzmax)=curr(k,ll)/3.e9
         enddo
         enddo
         istatus= NF_INQ_VARID(ncid,'curr',vid)
         call ncvpt_doubl2(ncid,vid,start_r0k,count_r0k,tem1,istatus)
         !YuP[07-2017] added
         !do k=1,ngen
         !do ll=1,lrzmax
         !      tem1(ll+(k-1)*lrzmax)=den_fsa(k,ll)
         !enddo
         !enddo
         !istatus= NF_INQ_VARID(ncid,'den_fsa',vid)
         !call ncvpt_doubl2(ncid,vid,start_r0k,count_r0k,tem1,istatus)
         !YuP[07-2017] added
         do k=1,ngen
         do ll=1,lrzmax
               tem1(ll+(k-1)*lrzmax)=energym(k,ll)
         enddo
         enddo
         istatus= NF_INQ_VARID(ncid,'energym',vid)
         call ncvpt_doubl2(ncid,vid,start_r0k,count_r0k,tem1,istatus)
      endif ! ngen

      kk=1
      if (kelecg.ne.0) kk=kelecg
      do ll=1,lrzmax
         tr(ll)=currm(kk,ll)/3.e9
      enddo
      if (kelecg.ne.0) then ! e as gen.species
         istatus= NF_INQ_VARID(ncid,'currm_e',vid)
      else ! ion as general species
         istatus= NF_INQ_VARID(ncid,'currm_i',vid)
      endif
      call ncvpt_doubl0(ncid,vid,start(3),r0_count(1),tr(1:),istatus)

      do ll=1,lrzmax
         tr(ll)=restp(nch(ll),ll)
      enddo
      istatus= NF_INQ_VARID(ncid,'restp',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,start(3),r0_count(1),tr(1:),istatus)

      do ll=1,lrzmax
         tr(ll)=restnp(nch(ll),ll)
      enddo
      istatus= NF_INQ_VARID(ncid,'restnp',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,start(3),r0_count(1),tr(1:),istatus)

      do ll=1,lrzmax
         tr(ll)=sptzrp(nch(ll),ll)
      enddo
      istatus= NF_INQ_VARID(ncid,'sptzrp',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,start(3),r0_count(1),tr(1:),istatus)

      istatus= NF_INQ_VARID(ncid,'rovsc',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,start(3),r0_count(1),rovsc,istatus)

      istatus= NF_INQ_VARID(ncid,'rovsc_hi',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,start(3),r0_count(1),rovsc_hi,istatus)

      istatus= NF_INQ_VARID(ncid,'zreskim',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,start(3),r0_count(1),zreskim,istatus)

      istatus= NF_INQ_VARID(ncid,'taueeh',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,start(3),r0_count(1),taueeh,istatus)

      istatus= NF_INQ_VARID(ncid,'nuestar',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,start(3),r0_count(1),starnue,istatus)

!BH110320c     Following only set up for ngen=1
!BH110320      if (ngen.ge.2)
!BH110320     1     WRITE(*,*)'netcdfrw2: Tot pwrs only set up for ngen=1'
      if (13*lrz .gt. iyjx2) stop 'netcdfrw2: Need 13*lrz<(iy+2)*(jx+2)'
      do k=1,ngen
         kkk=(k-1)*13*lrz
         call bcast(tem1,zero,13*lrz)
         do ll=1,lrz
            kk=1
            tem1(ll+(kk-1)*lrz+kkk)=entr(k,-1,ll)
            kk=2
            tem1(ll+(kk-1)*lrz+kkk)=entr(k,0,ll)
            kk=3
            tem1(ll+(kk-1)*lrz+kkk)=entr(k,2,ll)
            kk=4
            tem1(ll+(kk-1)*lrz+kkk)=entr(k,1,ll)
            kk=5
            tem1(ll+(kk-1)*lrz+kkk)=entr(k,3,ll)
            kk=6
            tem1(ll+(kk-1)*lrz+kkk)=entr(k,5,ll)
            kk=7
            tem1(ll+(kk-1)*lrz+kkk)=entr(k,6,ll)
            kk=8
            tem1(ll+(kk-1)*lrz+kkk)=entr(k,7,ll)
            kk=9
            tem1(ll+(kk-1)*lrz+kkk)=entr(k,8,ll)
            kk=10
            tem1(ll+(kk-1)*lrz+kkk)=entr(k,11,ll)
            kk=11
            tem1(ll+(kk-1)*lrz+kkk)=entr(k,10,ll)
            kk=12
            tem1(ll+(kk-1)*lrz+kkk)=entr(k,12,ll)
            kk=13
            tem1(ll+(kk-1)*lrz+kkk)=entr(k,4,ll)
         enddo
      enddo
      istatus= NF_INQ_VARID(ncid,'powers',vid)  !-YuP: NetCDF-f77 get vid

      call ncvpt_doubl2(ncid,vid,start_powers,count_powers,tem1,istatus)

!BH110320c     Following only set up for ngen=1
!BH110320      if (ngen.ge.2)
!BH110320    1     WRITE(*,*)'netcdfrw2: Tot pwrs only set up for ngen=1'
      if (lrz .gt. iyjx2) stop 'netcdfrw2: Need lrz<(iy+2)*(jx+2)'
      do k=1,ngen
            kkk=(k-1)*13
            call bcast(tem1,zero,13)
            kk=1
            tem1(kk+kkk)=entrintr(k,-1)
            kk=2
            tem1(kk+kkk)=entrintr(k,0)
            kk=3
            tem1(kk+kkk)=entrintr(k,2)
            kk=4
            tem1(kk+kkk)=entrintr(k,1)
            kk=5
            tem1(kk+kkk)=entrintr(k,3)
            kk=6
            tem1(kk+kkk)=entrintr(k,5)
            kk=7
            tem1(kk+kkk)=entrintr(k,6)
            kk=8
            tem1(kk+kkk)=entrintr(k,7)
            kk=9
            tem1(kk+kkk)=entrintr(k,8)
            kk=10
            tem1(kk+kkk)=entrintr(k,11)
            kk=11
            tem1(kk+kkk)=entrintr(k,10)
            kk=12
            tem1(kk+kkk)=entrintr(k,12)
            kk=13
            tem1(kk+kkk)=entrintr(k,4)
      enddo  !  on k=1,ngen
      istatus= NF_INQ_VARID(ncid,'powers_int',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl0(ncid,vid,start_powers(2),count_powers(2), &
                        tem1,istatus)

      if (sigmamod .eq. "enabled") then ! n>0: save sigftt
         istatus= NF_INQ_VARID(ncid,'sigftt',vid)  !-YuP: NetCDF-f77 get vid
         do lsig=1,4
            tem2(lsig)=sigftt(nch(1),lsig)
         enddo
         call ncvpt_doubl0(ncid,vid,start_fus(2),count_fus(2), &
                           tem2,istatus)
         call check_err(istatus)
      endif ! sigmamod


!      if (npa_diag .eq. "ncdf_all") then   NEED at least do first and last
      if (npa_diag .ne. "disabled" .and. numrec1.ge.2) then

         istatus= NF_INQ_VARID(ncid,'eflux_npa',vid)
         call pack21(eflux,1,nena,1,nva,wkpack,nen_npa,nv_npa)
         call ncvpt_doubl2(ncid,vid,start_npa,count_npa,wkpack,istatus)

         istatus= NF_INQ_VARID(ncid,'efluxt',vid)  !-YuP: NetCDF-f77 get vid
         call ncvpt_doubl2(ncid,vid,start_npa(2),count_npa(2),efluxt(1), &
              istatus)

      endif

      if (softxry .eq. "ncdf_all") then

         istatus= NF_INQ_VARID(ncid,'eflux',vid)  !-YuP: NetCDF-f77 get vid
         call pack21(eflux,1,nena,1,nva,wkpack,nen,nv)
         call ncvpt_doubl2(ncid,vid,start_xr,count_xr,wkpack,istatus)

         istatus= NF_INQ_VARID(ncid,'efluxt',vid)  !-YuP: NetCDF-f77 get vid
         call ncvpt_doubl0(ncid,vid,start_xr(2),count_xr(2),efluxt(1), &
              istatus)

      endif


      if (npa_diag .eq. "ncdf_all") then

         istatus= NF_INQ_VARID(ncid,'eflux_npa',vid)
         call pack21(eflux,1,nena,1,nva,wkpack,nen_npa,nv_npa)
         call ncvpt_doubl2(ncid,vid,start_npa,count_npa,wkpack,istatus)

         istatus= NF_INQ_VARID(ncid,'efluxt',vid)
         call ncvpt_doubl0(ncid,vid,start_npa(2),count_npa(2),efluxt(1), &
              istatus)

      endif  ! On npa_diag .eq. "ncdf_all"


!YuP      if ( netcdfshort.eq.'longer_f' ) then
      if ( (netcdfshort.eq.'longer_f').or. &
           ((netcdfshort.eq.'lngshrtf').and.(isave.ne.0))  ) then

         !YuP[2018-09-27] Added option to save f() only at nsave() steps,
         !rather than at every step.
         !write(*,*)'netcdfrw2: netcdfshort,n,isave=',netcdfshort,n,isave
         ! Here: n>0
         istatus= NF_INQ_VARID(ncid,'f',vid)  !-YuP: NetCDF-f77 get vid

         if (ngen.eq.1) then

           if(netcdfshort.eq.'lngshrtf') then
             start1(4)=numrecsave !YuP[2018-09-28]netcdfshort.eq.'lngshrtf'
           else ! 'longer_f'
             start1(4)=numrec1 ! usual (saved at every t step)
           endif
           do ll=1,lrz
              do j=1,jx
                do i=1,iy
                  temp1(i,j)=f(i,j,1,lrindx(ll))
                enddo
              enddo
              call pack21(temp1,0,iyp1,0,jxp1,wkpack,iy,jx)
              start1(3)=ll
              call ncvpt_doubl2(ncid,vid,start1(1:4),count1(1:4), &
                                wkpack,istatus)
               !So, here start1={1,   1,   ll,  numrec1(or numrecsave)}
               !         count1={iy,  jx,  1,   1}
               !         dimsf= {ydim,xdim,rdim,tdim(or tsavedim)}
           enddo
           start1(4)=numrec1 ! restore

         else ! ngen>1

           if(netcdfshort.eq.'lngshrtf') then
             startg(5)=numrecsave !YuP[2018-09-28]netcdfshort.eq.'lngshrtf'
           else ! 'longer_f'
             startg(5)=numrec1 ! usual (saved at every t step)
           endif
           do k=1,ngen
            do ll=1,lrz
               do j=1,jx
                  do i=1,iy
                     temp1(i,j)=f(i,j,k,lrindx(ll))
                  enddo
               enddo
               call pack21(temp1,0,iyp1,0,jxp1,wkpack,iy,jx)
               startg(3)=ll
               startg(4)=k
               call ncvpt_doubl2(ncid,vid,startg,countg,wkpack,istatus)
               !So, here startg={1,   1,   ll,  k,   numrec1(or numrecsave)}
               !         countg={iy,  jx,  1,   1,   1}
               !         dimsg= {ydim,xdim,rdim,gdim,tdim(or tsavedim)}
            enddo  !  On ll
           enddo  !  On k
           startg(3)=1
           startg(4)=1
           startg(5)=numrec1 ! restore

         endif ! ngen

      endif ! netcdfshort

!BH131030:  Could replace above k=1 output of distn by following
!BH131030:  output of f for all k, but for now leave it out as
!BH131030:  the netcdf file many easily get unwieldly large.
!YuP[2018-09-28] BUT WE MUST DO IT because 'f' was defined
!YuP             with either dims() or dimg(), depending on ngen.


      if (netcdfshort.ne.'long_jp') then

         if (ngen.eq.1) then

         k=1
         do ll=1,lrz
            do j=1,jx
               i=j+(lrindx(ll)-1)*jx
               tem1(i)=currv(j,k,lrindx(ll))/3.e9
            enddo
         enddo
         istatus= NF_INQ_VARID(ncid,'currv',vid)  !-YuP: NetCDF-f77 get vid
         call ncvpt_doubl0(ncid,vid,start(2),count(2),tem1,istatus)

         do ll=1,lrz
            do j=1,jx
               i=j+(lrindx(ll)-1)*jx
               tem1(i)=pwrrf(j,k,lrindx(ll))
            enddo
         enddo
         istatus= NF_INQ_VARID(ncid,'pwrrf',vid)  !-YuP: NetCDF-f77 get vid
         call ncvpt_doubl0(ncid,vid,start(2),count(2),tem1,istatus)

         else  !  ngen.ge.2

         do k=1,ngen
         do ll=1,lrz
            do j=1,jx
               i=j+(lrindx(ll)-1)*jx
               tem1(i)=currv(j,k,lrindx(ll))/3.e9
            enddo
         enddo
         startg(4)=k
         countg(3)=lrz
         istatus= NF_INQ_VARID(ncid,'currv',vid)  !-YuP: NetCDF-f77 get vid
         call ncvpt_doubl0(ncid,vid,startg(2),countg(2),tem1,istatus)

         do ll=1,lrz
            do j=1,jx
               i=j+(lrindx(ll)-1)*jx
               tem1(i)=pwrrf(j,k,lrindx(ll))
            enddo
         enddo
         startg(4)=k
         countg(3)=lrz
         istatus= NF_INQ_VARID(ncid,'pwrrf',vid)  !-YuP: NetCDF-f77 get vid
         call ncvpt_doubl0(ncid,vid,startg(2),countg(2),tem1,istatus)
         enddo  !  on k=1,ngen

         endif  !  on ngen
         startg(4)=1 ! restored
         countg(3)=1 ! restored

      endif ! 'long_jp'



! --- endif (n.gt.0) ---     if at l 3378
      endif








!.......................................................................
!cl    3.3 Variables saved only at last time-step

! --- begin if ---
      if (n.eq.nstop) then  !endif at line 4139

! --- full update of data file (last time-step only) ---

      if (sigmamod.eq.'enabled') then ! n=nstop: save fuspwrvt,flux_neutron_f,...
         istatus= NF_INQ_VARID(ncid,'fuspwrvt',vid)  !-YuP: NetCDF-f77 get vid
         call ncvpt_doubl0(ncid,vid,(1),4,fuspwrvt,istatus)
         call check_err(istatus)
         istatus= NF_INQ_VARID(ncid,'fuspwrv',vid)  !-YuP: NetCDF-f77 get vid
         do ii=1,4
            do ll=1,lrzmax
               ielem=ll+(ii-1)*lrzmax
               tem2(ielem)=fuspwrv(ii,ll)
            enddo
         enddo
         call ncvpt_doubl2(ncid,vid,start_fus,count_fus,tem2,istatus)
         call check_err(istatus)
      endif

      if(ampfmod.eq.'enabled')then !elecfldn(0:,0:,0:) saved at the last t step
         istatus= NF_INQ_VARID(ncid,'elecfldn',vid)
         write(*,*)'netcdfrw2: n=1,it=0,elecfldn(:,n,0)*300=', &
                                    elecfldn(0:lrz+1,1,0)*300
         write(*,*)'netcdfrw2: n=1,it=1,elecfldn(:,n,1)*300=', &
                                    elecfldn(0:lrz+1,1,1)*300
         call ncvpt_doubl2(ncid,vid,start_elecfldn,count_elecfldn, &
              elecfldn(0:,0:,0:),istatus)
      endif

      if ( (netcdfshort.ne.'enabled')  .and. &
           (netcdfshort.ne.'longer_f') .and. &
           (netcdfshort.ne.'lngshrtf')        ) then ! here: n=nstop

         istatus= NF_INQ_VARID(ncid,'f',vid)  !-YuP: NetCDF-f77 get vid
!BH011221: This storage is set up for constant iy as function of radius.
!BH011221: Needs generalizing. Should we store in reg array to max(iy_)?
!BH011221: For now, simply stop.
         do ll=1,lrz
            if (iy_(ll).ne.iymax) &
                 stop 'netcdfrw2: Cant handle iy.ne.iymax'
         enddo

         if (ngen.eq.1) then

            k=1 !
            if (tavg.ne."disabled") then
              PRINT *,'netcdfrw2 WARNING: favg is saved (NOT f)'
              PRINT *,'netcdfrw2 WARNING: but if time(nstop) < tavg1(1)'
              PRINT *,'netcdfrw2 WARNING: then favg == f.'
            endif
            WRITE(*,*) &
              'netcdfrw2[netcdfshort=disabled]:Write f into mnemonic.nc'
            WRITE(*,*)'netcdfrw2_4060: For checkup SUM(f),SUM(gone)=', &
                       SUM(f),SUM(gone)
            WRITE(*,*)'netcdfrw2_4060: For checkup MIN(f),MAX(f)=', &
                       MINVAL(f),MAXVAL(f)
            do ll=1,lrz
               if (tavg.eq."disabled") then
                  do j=1,jx
                  do i=1,iy
                     temp1(i,j)=f(i,j,k,lrindx(ll))
                     if(gone(i,j,k,lrindx(ll)).lt.-0.1) then
                        temp1(i,j)=em90
                     endif
                  enddo
                  enddo
               else  !On tavg = enabled
                  do j=1,jx
                  do i=1,iy
                     temp1(i,j)=favg(i,j,k,lrindx(ll))
                     !if(gone(i,j,k,lrindx(ll)).lt.-0.1) temp1(i,j)=em90
                     !Note: because of gone, the saved favg can be
                     !smaller than the original favg (removed loss cone).
                     !Comment the above if(gone...)
                     !if you want to leave favg intact.
                  enddo
                  enddo
                  WRITE(*,'(a,i6,2e13.4)') &
                  'nstop.netcdfrw2/tavg=en. ll, sumij(favg), sumij(f)', &
                   ll,sum(temp1),sum(f(:,:,k,lrindx(ll)))
               endif  !On tavg
               call pack21(temp1,0,iyp1,0,jxp1,wkpack,iy,jx)
               start1(3)=ll
               call ncvpt_doubl2(ncid,vid,start1(1:3),count1(1:3), &
                                 wkpack,istatus)
               !So, here start1={1,   1,   ll}
               !         count1={iy,  jx,  1}
               !         dimsf= {ydim,xdim,rdim}
            enddo ! ll
            start1(3)=1 ! restore

         else  !  ngen.ge.1

            do k=1,ngen
            do ll=1,lrz
               if (tavg.eq."disabled") then
               do j=1,jx
                  do i=1,iy
                     temp1(i,j)=f(i,j,k,lrindx(ll))
                  enddo
               enddo
               else  !On tavg
               do j=1,jx
                  do i=1,iy
                     temp1(i,j)=favg(i,j,k,lrindx(ll))
                  enddo
               enddo
               endif  !On tavg
               call pack21(temp1,0,iyp1,0,jxp1,wkpack,iy,jx)
               startg(3)=ll
               startg(4)=k
               call ncvpt_doubl2(ncid,vid,startg(1:4),countg(1:4), &
                                 wkpack,istatus)
               !So, here startg={1,   1,   ll,  k}
               !         countg={iy,  jx,  1,   1}
               !         dimsg= {ydim,xdim,rdim,gdim}
            enddo ! on ll
            enddo ! on k
            startg(3)=1
            startg(4)=1

         endif !  on ngen

      endif


      if (softxry .ne. "disabled" .and. softxry.ne."ncdf_all") then

         start_xr(3)=2
         istatus= NF_INQ_VARID(ncid,'eflux',vid)  !-YuP: NetCDF-f77 get vid
         call pack21(eflux,1,nena,1,nva,wkpack,nen,nv)
         call ncvpt_doubl2(ncid,vid,start_xr,count_xr,wkpack,istatus)

         istatus= NF_INQ_VARID(ncid,'efluxt',vid)  !-YuP: NetCDF-f77 get vid
         call ncvpt_doubl0(ncid,vid,start_xr(2),count_xr(2),efluxt(1), &
              istatus)

      endif

! --- endif (n.eq.nstop) ---
      endif  !if at l 4016

!-----------------------------------------------------------------------
!
!cl    4. Close netCDF file
!
      if (n.eq.nstop) then
         istatus = NF_CLOSE(ncid) !-YuP: NetCDF-f77
         call check_err(istatus)
      endif

      return
      end subroutine netcdfrw2
!
!

      subroutine check_err(iret)
      implicit none
      integer iret
      include 'netcdf.inc'
      if (iret .ne. NF_NOERR) then
         PRINT *, 'netCDF error:', iret
         PRINT *, NF_STRERROR(iret) ! print error explanation
         call abort
      endif
      end subroutine check_err

!
!
      subroutine netcdfmain
      use param_mod
      use cqlcomm_mod
      implicit none
      integer igrid ! local

!.......................................................................
!     Controls some netcdf output (could make it all).
!     Presently, just output of velocity-space flux data.
!.......................................................................
      save
!MPIINSERT_INCLUDE

!MPIINSERT_IF_RANK_NE_0_RETURN
 ! save data on mpirank.eq.0 only

      if (n.ne.nstop) return
      if ((netcdfvecs.eq."irzplt" .and. mplot(l_).eq."enabled") &
           .or. (netcdfvecs.eq."all")) then

         if (netcdfvecal.ne."disabled") then
            igrid=0
            if(netcdfvecal.eq."x-theta") igrid=1
            call netcdfvec(4,igrid)
         endif
         if (netcdfvecrf.ne."disabled") then
            igrid=0
            if(netcdfvecrf.eq."x-theta") igrid=1
            call netcdfvec(3,igrid)
         endif
         if (netcdfvece.ne."disabled") then
            igrid=0
            if(netcdfvece.eq."x-theta") igrid=1
            call netcdfvec(2,igrid)
         endif
         if (netcdfvecc.ne."disabled") then
            igrid=0
            if(netcdfvecc.eq."x-theta") igrid=1
            call netcdfvec(1,igrid)
         endif

      endif

      return
      end subroutine netcdfmain
!
!
      subroutine netcdfvec(lefct,igrid)
      use advnce_mod, only : hfu, hfi
      use bcast_mod, only : bcast
      use coefefad_mod, only : coefefad
      use coeffpad_mod, only : coeffpad
      use coefmidt_mod, only : coefmidt
      use coefmidv_mod, only : coefmidv
      use coefrfad_mod, only : coefrfad
      use coefstup_mod, only : coefstup
      use cqlcomm_mod
      use diagentr_mod, only : gfi
      use diagentr_mod, only : gfu
      use param_mod
      use prppr_mod, only : prppr
      use r8subs_mod, only : dcopy

      implicit none

!...................................................................
!
!     This subroutine writes out flux data in netCDF format.
!     It is called if one of netcdfvecXX = .ne."disabled".
!     igrid=1:  indicates netcdfvecXX="u-theta", then write the
!       x,theta components of the fluxes on the code x,theta grid
!       (same as distribution functions).
!     igrid=0: indicates netcdfvecXX = .ne."disabled" and .ne."u-theta",
!       then write velocity-space fluxes suitable for
!       vector plot in (x-par,x-perp) coordinates, on an xpar,xperp
!       grid.
!
!...................................................................

      save

      include 'netcdf.inc'

      integer ncid,vid,istatus
      integer ipxydim,jpxydim
      integer dims(4),count(4),start(4)
      integer chardim,char64dim
      integer xdim,ydim,rdim,gdim
      integer y_dims(2),y_count(2),y_start(2)

      character*8 target
      character*18 fluxcmpt(4)

      integer lfirst(4)
      data lfirst/4*0/
      data fluxcmpt/"collisional flux  ","electric fld flux ", &
                    "rf diffusion flux ","sum of fluxes     "/

      real*8, allocatable :: wkpack(:) ! local working array for pack21
      
      integer :: lefct, igrid ! input
      integer :: nwkpack, istat, ll, llcount, i, j, k ! local
      integer :: n_netcdf ! local
      real*8 xll,xlu,xlp,xpl,xpu,xmaxq,pltlimm,pltlimmm ! local scalars, input for other functions
      real*8 xrf ! local scalar, input for call coefrfad()
      real*8 rya_netcdf(lrza) ! local, and input for ncvpt_doubl2
      integer :: ll_netcdf(lrza),itl_netcdf(lrza),itu_netcdf(lrza) ! local, and input for ncvpt_int2


!     another way to do this, when the code is a mess.
      character(len=8) :: cqlpmod
      integer :: lrindx(0:lrorsa)
      integer :: lrz
      integer :: lrzmax
      character(len=256) :: mnemonic
      cqlpmod = setup0%cqlpmod
      lrindx = setup0%lrindx
      lrz = setup0%lrz
      lrzmax =  setup0%lrzmax
      mnemonic = setup0%mnemonic


      if ((lefct.lt.1) .or. (lefct.gt.4)) stop 'pltvec:lefct'

!...................................................................
!     Set up netCDF output on last time step, if
!     this is first call for a given value of lefct.
!     It is assumed that lefct is an integer.ge.1, and that
!     there are lrz (or nirzplt, if netcdfvecs="irzplt")
!     calls to this routine (corresponding
!     to each of the FP'd flux surfaces for each
!     value of lefct.
!
!     lefct=1  ==> collisional flux
!     lefct=2  ==> electric fld flux
!     lefct=3  ==> rf diffusion flux
!     lefct=4  ==> sum of fluxes
!
!     This routine can readily be adapted for
!     cqlpmod="enabled".
!...................................................................


      if (cqlpmod.eq."enabled") then
         stop "netcdfrw2: Need to adapt sub for cqlpmod.eq.enabled"
      endif

      if (.NOT.ALLOCATED(wkpack)) then ! allocate working array for pack21
         nwkpack=max(iyjx2, jpxy*ipxy) + 10   ! +10 just in case
         allocate(wkpack(nwkpack),STAT=istat)
         call bcast(wkpack,zero,SIZE(wkpack))
      endif

!     lfirst(lefct) is an indicator of whether this is the first
!     call for given lefct.
      if (netcdfnm.ne."disabled" .and. n.eq.nstop &
          .and. lfirst(lefct).ne.lefct) then

         lfirst(lefct)=lefct

         write(t_,235) mnemonic(1:length_char(mnemonic)),lefct
 235     format(a,"_flux_",i1,".nc")

         istatus = NF_CREATE(t_, NF_CLOBBER, ncid) !-YuP: NetCDF-f77
         call check_err(istatus)

         call ncaptc2(ncid,NCGLOBAL,'title',NCCHAR,18, &
          fluxcmpt(lefct),istatus)

!        Put radial flux surface numbers in ll_netcdf(),
!        number of surfaces in n_netcdf, normalized radii in rya_netcdf(),
!        indices of lower and upper trapped-passing boundaries in
!        itl_netcdf() and itu_netcdf().
         if (netcdfvecs.eq."irzplt") then
            n_netcdf=nirzplt
            llcount=0
            do ll=1,lrz
               if (mplot(ll).eq."enabled") then
                  llcount=llcount+1
                  ll_netcdf(llcount)=ll
                  rya_netcdf(llcount)=rya(lrindx(ll))
                  itl_netcdf(llcount)=itl_(lrindx(ll))
                  itu_netcdf(llcount)=itu_(lrindx(ll))
               endif
            enddo
         else
            n_netcdf=lrz
            do ll=1,lrz
               ll_netcdf(ll)=ll
               rya_netcdf(ll)=rya(lrindx(ll))
               itl_netcdf(ll)=itl_(lrindx(ll))
               itu_netcdf(ll)=itu_(lrindx(ll))
           enddo
         endif

         if(igrid.eq.0) then

         istatus= NF_DEF_DIM(ncid, 'jpxydim', jpxy,     jpxydim)  !-YuP: NetCDF-f77
         istatus= NF_DEF_DIM(ncid, 'ipxydim', ipxy,     ipxydim)
         istatus= NF_DEF_DIM(ncid, 'rdim',    n_netcdf, rdim)
         istatus= NF_DEF_DIM(ncid, 'gdim',    ngen,     gdim)
         istatus= NF_DEF_DIM(ncid, 'chardim',     8,    chardim)
         istatus= NF_DEF_DIM(ncid, 'char64dim',  64,    char64dim)

         dims(1)=jpxydim
         dims(2)=ipxydim
         dims(3)=rdim
         dims(4)=gdim

         start(1)=1
         start(2)=1
         start(3)=1
         start(4)=1

         count(1)=jpxy
         count(2)=ipxy
         count(3)=1
         count(1)=1


      vid=ncvdef0(ncid,'mnemonic',NCCHAR,1,char64dim,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,23, &
              'Mnemonic run identifier',istatus)

      vid=ncvdef0(ncid,'grid_type',NCCHAR,1,chardim,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,35, &
              'Indicates xpar-prp or x-theta grid',istatus)

      vid=ncvdef0(ncid,'lrz',NCLONG,0,0,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,34, &
              'Number of FPd radial surface bins',istatus)
         call check_err(istatus)

      vid=ncvdef0(ncid,'n_netcdf',NCLONG,0,0,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,44, &
              'Number of radial surfaces for output (=rdim)',istatus)
         call check_err(istatus)

      vid=ncvdef0(ncid,'ll_netcdf',NCLONG,1,rdim,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,40, &
              'Radial bin numbers where data is stored',istatus)
         call check_err(istatus)

      vid=ncvdef0(ncid,'jpxy',NCLONG,0,0,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,45, &
              'Number of parallel momentum points (=jpxydim)',istatus)
         call check_err(istatus)

      vid=ncvdef0(ncid,'ipxy',NCLONG,0,0,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,41, &
           'Number of perp momentum points (=ipxydim)',istatus)
         call check_err(istatus)

      vid=ncvdef0(ncid,'xll',NCDOUBLE,0,0,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,32, &
              'Lower parallel momentum-per-mass, normalized',istatus)

      vid=ncvdef0(ncid,'xlu',NCDOUBLE,0,0,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,32, &
              'Upper parallel momentum-per-mass, normalized',istatus)

      vid=ncvdef0(ncid,'xpl',NCDOUBLE,0,0,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,40, &
              'Lower perp momentum-per-mass, normalized',istatus)

      vid=ncvdef0(ncid,'xpu',NCDOUBLE,0,0,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,40, &
              'Upper perp momentum-per-mass, normalized',istatus)

      vid=ncvdef0(ncid,'xpar',NCDOUBLE,(1),jpxydim,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,32, &
              'Normalized par momentum-per-mass',istatus)

      vid=ncvdef0(ncid,'xperp',NCDOUBLE,(1),ipxydim,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,33, &
              'Normalized perp momentum-per-mass',istatus)

      vid=ncvdef0(ncid,'rhomax',NCDOUBLE,0,0,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,31, &
              'Generalized plasma minor radius',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,3, &
              'cms',istatus)

      vid=ncvdef0(ncid,'rya_netcdf',NCDOUBLE,(1),rdim,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,37, &
              'Normalized radial mesh at bin centers',istatus)

      if (ngen.eq.1) then !Maintaining backwards compatability
      vid=ncvdef2(ncid,'gamma_par',NCDOUBLE,3,dims,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,28, &
              'Parallel momentum-space flux',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,16, &
              'Normalized units',istatus)
         call ncaptc2(ncid,vid,'comment',NCCHAR,39, &
              'Facility set up only for single species',istatus)

      vid=ncvdef2(ncid,'gamma_perp',NCDOUBLE,3,dims,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,33, &
              'Perpendicular momentum-space flux',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,16, &
              'Normalized units',istatus)
      else  !  ngen.ge.2
      vid=ncvdef2(ncid,'gamma_par',NCDOUBLE,4,dims,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,28, &
              'Parallel momentum-space flux',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,16, &
              'Normalized units',istatus)
         call ncaptc2(ncid,vid,'comment',NCCHAR,39, &
              'Facility set up for one or more species',istatus)

      vid=ncvdef2(ncid,'gamma_perp',NCDOUBLE,4,dims,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,33, &
              'Perpendicular momentum-space flux',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,16, &
              'Normalized units',istatus)
         call ncaptc2(ncid,vid,'comment',NCCHAR,39, &
              'Facility set up for one or more species',istatus)
      endif   !  on ngen

         else  ! igrid.ne.0 case


!     Maximum iy as function of radius is iymax, set in ainsetpa.f:

!BH021028: This storage is set up for constant iy as function of radius.
!BH021028: Needs generalizing. Should we store in reg array to max(iy_)?
!BH021028: For now, simply stop.
         do ll=1,lrz
            if (iy_(ll).ne.iymax) &
                 stop 'netcdfrw2: Cant handle iy.ne.iymax'
         enddo

         istatus= NF_DEF_DIM(ncid, 'xdim',    jx,       xdim)  !-YuP: NetCDF-f77
         istatus= NF_DEF_DIM(ncid, 'ydim',    iy,       ydim)
         istatus= NF_DEF_DIM(ncid, 'rdim',    n_netcdf, rdim)
         istatus= NF_DEF_DIM(ncid, 'gdim',    ngen,     gdim)
         istatus= NF_DEF_DIM(ncid, 'chardim',     8,    chardim)
         istatus= NF_DEF_DIM(ncid, 'char64dim',  64,    char64dim)

         dims(1)=ydim
         dims(2)=xdim
         dims(3)=rdim
         dims(4)=gdim

         start(1)=1
         start(2)=1
         start(3)=1
         start(4)=1

         count(1)=iy
         count(2)=jx
         count(3)=1

         y_dims(1)=ydim
         y_dims(2)=rdim

         y_start(1)=1
         y_start(2)=1

         y_count(1)=iy
         y_count(2)=1


      vid=ncvdef0(ncid,'mnemonic',NCCHAR,1,char64dim,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,23, &
              'Mnemonic run identifier',istatus)

      vid=ncvdef0(ncid,'grid_type',NCCHAR,1,chardim,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,36, &
              'Indicates xpar-prp or x-theta grid',istatus)

      vid=ncvdef0(ncid,'lrz',NCLONG,0,0,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,34, &
              'Number of FPd radial surface bins',istatus)
         call check_err(istatus)

      vid=ncvdef0(ncid,'n_netcdf',NCLONG,0,0,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,44, &
              'Number of radial surfaces for output (=rdim)',istatus)
         call check_err(istatus)

      vid=ncvdef0(ncid,'ll_netcdf',NCLONG,1,rdim,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,40, &
              'Radial bin numbers where data is stored',istatus)
         call check_err(istatus)

      vid=ncvdef0(ncid,'jx',NCLONG,0,0,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,35, &
              'momentum-per-mass dimension (=xdim)',istatus)
         call check_err(istatus)

      vid=ncvdef0(ncid,'x',NCDOUBLE,(1),xdim,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,28, &
              'normalized momentum-per-mass',istatus)

      vid=ncvdef0(ncid,'vnorm',NCDOUBLE,0,0,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,33, &
              'velocity (momentum-per-mass) norm',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,7, &
              'cms/sec',istatus)

      vid=ncvdef0(ncid,'enorm',NCDOUBLE,0,0,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,20, &
              'Energy normalization',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,3, &
              'keV',istatus)

      vid=ncvdef0(ncid,'iy',NCLONG,0,0,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,33, &
              'max Pitch angle dimension (=ydim)',istatus)
         call check_err(istatus)

      vid=ncvdef2(ncid,'y',NCDOUBLE,2,y_dims,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,36, &
              'pitch angle at n_netcdf flux surfaces',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,7, &
              'radians',istatus)

      vid=ncvdef0(ncid,'itl_netcdf',NCLONG,1,rdim,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,32, &
              'lower trapped-passing bndy index',istatus)

      vid=ncvdef0(ncid,'itu_netcdf',NCLONG,1,rdim,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,32, &
              'upper trapped-passing bndy index',istatus)

      vid=ncvdef0(ncid,'rhomax',NCDOUBLE,0,0,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,31, &
              'Generalized plasma minor radius',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,3, &
              'cms',istatus)

      vid=ncvdef0(ncid,'rya_netcdf',NCDOUBLE,(1),rdim,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,37, &
              'Normalized radial mesh at bin centers',istatus)

      if (ngen.eq.1) then  !Maintaining backwards compatability
      vid=ncvdef2(ncid,'gamma_x',NCDOUBLE,3,dims,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,28, &
              'x-compnt momentum-space flux',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,16, &
              'Normalized units',istatus)
         call ncaptc2(ncid,vid,'comment',NCCHAR,39, &
              'Facility set up for one or more species',istatus)

      vid=ncvdef2(ncid,'gamma_theta',NCDOUBLE,3,dims,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,32, &
              'Theta-compnt momentum-space flux',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,16, &
              'Normalized units',istatus)
         call ncaptc2(ncid,vid,'comment',NCCHAR,39, &
              'Facility set up for one or more species',istatus)
      else  !  ngen.ge.2
      vid=ncvdef2(ncid,'gamma_x',NCDOUBLE,4,dims,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,28, &
              'x-compnt momentum-space flux',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,16, &
              'Normalized units',istatus)
         call ncaptc2(ncid,vid,'comment',NCCHAR,39, &
              'Facility set up for one or more species',istatus)

      vid=ncvdef2(ncid,'gamma_theta',NCDOUBLE,4,dims,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,32, &
              'Theta-compnt momentum-space flux',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,16, &
              'Normalized units',istatus)
         call ncaptc2(ncid,vid,'comment',NCCHAR,39, &
              'Facility set up for one or more species',istatus)
      endif  !  on ngen


         endif   !end igrid options



         istatus= NF_ENDDEF(ncid) !-YuP: NetCDF-f77
         call check_err(istatus)
         istatus = NF_CLOSE(ncid) !-YuP: NetCDF-f77

      endif



!     The following calculation of velocity space fluxes is copied
!     from subroutine pltvec.

      do 190 k=1,ngen
         if (tandem.eq."enabled" .and. k.eq.kionn) then
            xll=-xlwr
            xlu=xlwr
            xpl=0.
            xpu=xlwr
            xmaxq=xlwr
            target="ionmesh"
         else
            xll=-xmax
            xlu=xmax
            xpl=0.
            xpu=xmax
            xmaxq=xmax
            target="mainmesh"
         endif

!     If pltlim.ne."disabled", then plot in xpar,xprp-space up
!     to appropriate limit:
         if (pltlim.ne."disabled") then
            if (pltlim.eq.'x') then
               pltlimmm=pltlimm
            elseif (pltlim.eq.'u/!') then
               pltlimmm=pltlimm*cnorm
            elseif (pltlim.eq.'energy') then
               pltlimmm=sqrt((1.+pltlimm/restmkev)**2-1.)*cnorm
            endif
            xll=-pltlimmm
            xlu=pltlimmm
            xpl=0.
            xpu=pltlimmm
            xmaxq=pltlimmm
         endif

!
         call bcast(da,zero,iyjxp1)
         call bcast(db,zero,iyjxp1)
         call bcast(dc,zero,iyjxp1)
         call bcast(dd,zero,iyp1jx)
         call bcast(de,zero,iyp1jx)
         call bcast(df,zero,iyp1jx)
         if (lefct .eq. 1) go to 50
         if (lefct .eq. 2) go to 60
         if (lefct .eq. 3) go to 80
         if (lefct .eq. 4) go to 100
 50      call coeffpad(k)



         go to 120
 60      if (abs(elecfld(lr_)) .lt. 1.e-09) go to 190
         call coefefad(k)



         go to 120
 80      continue
         xrf=0.
         if (n .lt. nonrf(k) .or. n .ge. noffrf(k)) go to 90
         call coefrfad(k,xrf)
 90      continue
         if (xrf.gt.0.) then



         endif
         go to 120
 100     continue
         call coefstup(k)



 120     continue

!...................................................................
!     The coefficients of the equation are currently defined on the
!     same mesh as the distribution function f. The fluxes are best
!     defined (from the point of view of differencing and particle
!     conservation) on mid meshpoints. We undertake here to
!     interpolate the coefficients as needed to either (i,j+1/2)
!     (velocity flux) or to (i+1/2,j) (theta flux).
!     Finally to enforce boundary conditions (zero flux in general
!     except at the pass/trapped boundary) certain coefficients
!     are zeroed out or suitably averaged at specific mesh points.
!     The numbers 1,2,3 appearing in the calls below signify
!     which coefficient is being treated.
!
!     first the velocity flux coefficients for gfi..
!...................................................................

         call coefmidv(da,1)
         call coefmidv(db,2)
         call coefmidv(dc,3)

!...................................................................
!     the theta flux coefficients for hfi..
!...................................................................

         call coefmidt(dd,1)
         call coefmidt(de,2)
         call coefmidt(df,3)
         if (lefct .eq. 3 .and. xrf .eq. 0) go to 190
         call bcast(temp5(0:iy+1,0:jx+1),zero,iyjx2)
         call bcast(temp4(0:iy+1,0:jx+1),zero,iyjx2)
!
!        In the following, tam2 and tam3 are u-space fluxes Gamma_x and
!          sin(theta)*Gamma_theta (in code units)
!          interpolated onto the code theta,x-mesh.
!
!          For igrid.eq.0, temp5 and temp4 are the
!          parallel, perpendicular components of flux Gamma,
!          respectively, on the x,theta-mesh.
!
!          For igrid.eq.1, temp5,temp4 are x,theta components,
!          respectively, of flux Gamma on the theta,x-mesh.
!
         if (implct .eq. "enabled") then
            do 140 i=2,iy-1
               do 141 j=2,jxm1
                  tam2(j)=-(gfi(i,j,k)+gfi(i,j-1,k))*.5/xsq(j)
                  tam3(j)=-(hfi(i,j,k,l_)+hfi(i-1,j,k,l_))*.5*xi(j)
 141           continue
               if (igrid.eq.1) then
                  do j=2,jxm1
                     temp5(i,j)=tam2(j)
                     temp4(i,j)=tam3(j)/sinn(i,l_)
                 enddo
               else
                  do j=2,jxm1
                     temp5(i,j)=tam2(j)*coss(i,l_)-tam3(j)
                     temp4(i,j)=tam2(j)*sinn(i,l_)+tam3(j)/tann(i,l_)
                  enddo
               endif
 140        continue
         else
            call dcopy(iyjx2,  f_(0:iy+1,0:jx+1,k,l_),1, &
                            temp1(0:iy+1,0:jx+1),1)
            call dcopy(iyjx2,fxsp(0:iy+1,0:jx+1,k,l_),1, &
                            temp2(0:iy+1,0:jx+1),1)
            do 240 j=2,jxm1
               do 241 i=2,iy-1
                  temp6(i,j)=-(gfu(i,j,k)+gfu(i,j-1,k))*.5/xsq(j)
 241           continue
 240        continue
            call dcopy(iyjx2,temp2(0:iy+1,0:jx+1),1, &
                             temp1(0:iy+1,0:jx+1),1)
            call dcopy(iyjx2,    f(0:iy+1,0:jx+1,k,l_),1, &
                             temp2(0:iy+1,0:jx+1),1)
            do 242 i=2,iy-1
               do 243 j=2,jxm1
                  tam3(j)=-(hfu(i,j,k,l_)+hfu(i-1,j,k,l_))*.5*xi(j)
 243           continue
               if (igrid.eq.1) then
                  do j=2,jxm1
                     temp5(i,j)=temp6(i,j)
                     temp4(i,j)=tam3(j)/sinn(i,l_)
                  enddo
               else
                  do j=2,jxm1
                     temp5(i,j)=temp6(i,j)*coss(i,l_)-tam3(j)
                     temp4(i,j)=temp6(i,j)*sinn(i,l_)+tam3(j)/tann(i,l_)
                  enddo
               endif
 242        continue
         endif
         do 244 j=1,jx
            temp5(itl,j)=temp5(itl-1,j)
            temp5(itu,j)=temp5(itu+1,j)
            temp4(itl,j)=temp4(itl-1,j)
            temp4(itu,j)=temp4(itu+1,j)
 244     continue

!
!     Above gives:
!     for igrid.eq.0,
!       flux flux_perp(u,theta) ~ temp4, flux_par ~ temp5.
!       These are output to the netcdf file.
!     for igrid.eq.1,
!       flux flux_u(theta,u)~temp5, flux_theta~temp4.
!
!

         if (igrid.eq.0) then

!
!       The following calls to prppr and dcopy put
!       flux_par into xhead, flux_perp into yhead, on an x,y-grid.
!       (xhead(jpxy,ipxy) has temp1 ptr,
!        yhead(jpxy,ipxy) has temp4 ptr),
!        and these are output to the netcdf file.

         call dcopy(iyjx2,temp5(0:iy+1,0:jx+1),1,temp3(0:iy+1,0:jx+1),1)

         call prppr(target,"norm",xll,xlu,xpl,xpu)

         call dcopy(iyjx2,temp2(0:iy+1,0:jx+1),1,temp1(0:iy+1,0:jx+1),1)

         call dcopy(iyjx2,temp4(0:iy+1,0:jx+1),1,temp3(0:iy+1,0:jx+1),1)

         call prppr(target,"norm",xll,xlu,xpl,xpu)

         call dcopy(iyjx2,temp2(0:iy+1,0:jx+1),1,temp4(0:iy+1,0:jx+1),1)

         endif  !end igrid.eq.0




!     gamma_par ~ temp1, gamma_prp ~ temp4
!     For netcdf data:  for each l_, compress into
!     jpxy*ipxy array and output.


         if (netcdfnm.ne."disabled" .and. n.eq.nstop) then


            write(t_,236) mnemonic(1:length_char(mnemonic)),lefct
 236        format(a,"_flux_",i1,".nc")
!     Open the correct netcdf file, to get ncid.
!-YuP            ncid = ncopn(t_,NCWRITE,istatus)
            istatus = NF_OPEN(t_, NF_WRITE, ncid) !-YuP: NetCDF-f77

!     Write some data at head of the file:

            if (igrid.eq.0) then

            istatus= NF_INQ_VARID(ncid,'mnemonic',vid)  !-YuP: NetCDF-f77 get vid
            ll=length_char(mnemonic)
            call ncvptc0(ncid,vid,1,ll,mnemonic,ll,istatus)

            istatus= NF_INQ_VARID(ncid,'grid_type',vid)  !-YuP: NetCDF-f77 get vid
            call ncvptc0(ncid,vid,1,8,'xpar-prp',8,istatus)

            istatus= NF_INQ_VARID(ncid,'lrz',vid)  !-YuP: NetCDF-f77 get vid
            call ncvpt_int0(ncid,vid,1,1,(lrz),istatus)

            istatus= NF_INQ_VARID(ncid,'n_netcdf',vid)  !-YuP: NetCDF-f77 get vid
            call ncvpt_int0(ncid,vid,1,1,(n_netcdf),istatus)

            istatus= NF_INQ_VARID(ncid,'ll_netcdf',vid)  !-YuP: NetCDF-f77 get vid
            call ncvpt_int0(ncid,vid,1,n_netcdf, &
                 ll_netcdf(1:n_netcdf),istatus)

            istatus= NF_INQ_VARID(ncid,'jpxy',vid)  !-YuP: NetCDF-f77 get vid
            call ncvpt_int0(ncid,vid,1,1,(jpxy),istatus)

            istatus= NF_INQ_VARID(ncid,'ipxy',vid)  !-YuP: NetCDF-f77 get vid
            call ncvpt_int0(ncid,vid,1,1,(ipxy),istatus)

            istatus= NF_INQ_VARID(ncid,'xll',vid)  !-YuP: NetCDF-f77 get vid
            call ncvpt_doubl0(ncid,vid,(1),1,(xll),istatus)

            istatus= NF_INQ_VARID(ncid,'xlu',vid)  !-YuP: NetCDF-f77 get vid
            call ncvpt_doubl0(ncid,vid,(1),1,(xlu),istatus)

            istatus= NF_INQ_VARID(ncid,'xpl',vid)  !-YuP: NetCDF-f77 get vid
            call ncvpt_doubl0(ncid,vid,(1),1,(xpl),istatus)

            istatus= NF_INQ_VARID(ncid,'xpu',vid)  !-YuP: NetCDF-f77 get vid
            call ncvpt_doubl0(ncid,vid,(1),1,(xpu),istatus)

            istatus= NF_INQ_VARID(ncid,'xpar',vid)  !-YuP: NetCDF-f77 get vid
            call ncvpt_doubl0(ncid,vid,(1),jpxy,xpar,istatus)

            istatus= NF_INQ_VARID(ncid,'xperp',vid)  !-YuP: NetCDF-f77 get vid
            call ncvpt_doubl0(ncid,vid,(1),ipxy,xperp,istatus)

            istatus= NF_INQ_VARID(ncid,'rhomax',vid)  !-YuP: NetCDF-f77 get vid
            call ncvpt_doubl0(ncid,vid,(1),1,(rhomax),istatus)

            istatus= NF_INQ_VARID(ncid,'rya_netcdf',vid)  !-YuP: NetCDF-f77 get vid
            call ncvpt_doubl0(ncid,vid,(1),n_netcdf, &
                 rya_netcdf(1:n_netcdf), istatus)



!     Write the fluxes:
!     (First, figure out bin number for netcdf file.)
            do ll=1,n_netcdf
               if (l_ .eq. ll_netcdf(ll)) start(3)=ll
            enddo

            call pack21(xhead,1,jpxy,1,ipxy,wkpack,jpxy,ipxy)
            istatus= NF_INQ_VARID(ncid,'gamma_par',vid)  !-YuP: NetCDF-f77 get vid
            call ncvpt_doubl2(ncid,vid,start,count,wkpack,istatus)

            call pack21(yhead,1,jpxy,1,ipxy,wkpack,jpxy,ipxy)
            istatus= NF_INQ_VARID(ncid,'gamma_perp',vid)  !-YuP: NetCDF-f77 get vid
            call ncvpt_doubl2(ncid,vid,start,count,wkpack,istatus)



            else   ! igrid.ne.0 case



            istatus= NF_INQ_VARID(ncid,'mnemonic',vid)  !-YuP: NetCDF-f77 get vid
            ll=length_char(mnemonic)
            call ncvptc0(ncid,vid,1,ll,mnemonic,ll,istatus)

            istatus= NF_INQ_VARID(ncid,'grid_type',vid)  !-YuP: NetCDF-f77 get vid
            call ncvptc0(ncid,vid,1,7,'x-theta',7,istatus)

            istatus= NF_INQ_VARID(ncid,'lrz',vid)  !-YuP: NetCDF-f77 get vid
            call ncvpt_int0(ncid,vid,1,1,(lrz),istatus)

            istatus= NF_INQ_VARID(ncid,'n_netcdf',vid)  !-YuP: NetCDF-f77 get vid
            call ncvpt_int0(ncid,vid,1,1,(n_netcdf),istatus)

            istatus= NF_INQ_VARID(ncid,'ll_netcdf',vid)  !-YuP: NetCDF-f77 get vid
            call ncvpt_int0(ncid,vid,1,n_netcdf, &
                 ll_netcdf(1:n_netcdf),istatus)

            istatus= NF_INQ_VARID(ncid,'jx',vid)  !-YuP: NetCDF-f77 get vid
            call ncvpt_int0(ncid,vid,1,1,(jx),istatus)

            istatus= NF_INQ_VARID(ncid,'x',vid)  !-YuP: NetCDF-f77 get vid
            call ncvpt_doubl0(ncid,vid,(1),jx,x,istatus)

            istatus= NF_INQ_VARID(ncid,'vnorm',vid)  !-YuP: NetCDF-f77 get vid
            call ncvpt_doubl0(ncid,vid,(1),1,(vnorm),istatus)

            istatus= NF_INQ_VARID(ncid,'enorm',vid)  !-YuP: NetCDF-f77 get vid
            call ncvpt_doubl0(ncid,vid,(1),1,(enorm),istatus)

            istatus= NF_INQ_VARID(ncid,'iy',vid)  !-YuP: NetCDF-f77 get vid
            call ncvpt_int0(ncid,vid,1,1,(iy),istatus)

            istatus= NF_INQ_VARID(ncid,'itl_netcdf',vid)  !-YuP: NetCDF-f77 get vid
            call ncvpt_int0(ncid,vid,1,n_netcdf, &
                 itl_netcdf(1:n_netcdf),istatus)

            istatus= NF_INQ_VARID(ncid,'itu_netcdf',vid)  !-YuP: NetCDF-f77 get vid
            call ncvpt_int0(ncid,vid,1,n_netcdf, &
                 itu_netcdf(1:n_netcdf),istatus)

            istatus= NF_INQ_VARID(ncid,'rhomax',vid)  !-YuP: NetCDF-f77 get vid
            call ncvpt_doubl0(ncid,vid,(1),1,(rhomax),istatus)

            istatus= NF_INQ_VARID(ncid,'rya_netcdf',vid)  !-YuP: NetCDF-f77 get vid
            call ncvpt_doubl0(ncid,vid,(1),n_netcdf, &
                 rya_netcdf(1:n_netcdf), istatus)


!     Write the fluxes (and the pitch angle mesh y):
!     (First, figure out bin number for netcdf file.)
            do ll=1,n_netcdf
               if (l_ .eq. ll_netcdf(ll)) then
                  start(3)=ll
                  y_start(2)=ll
               endif
            enddo
            start(4)=k  !For ngen.ge.2 case.

            istatus= NF_INQ_VARID(ncid,'y',vid)  !-YuP: NetCDF-f77 get vid
            call ncvpt_doubl2(ncid,vid,y_start,y_count,y(1:,l_),istatus)

            call pack21(temp5,0,iyp1,0,jxp1,wkpack,iy,jx)
            istatus= NF_INQ_VARID(ncid,'gamma_x',vid)  !-YuP: NetCDF-f77 get vid
            call ncvpt_doubl2(ncid,vid,start,count,wkpack,istatus)

            call pack21(temp4,0,iyp1,0,jxp1,wkpack,iy,jx)
            istatus= NF_INQ_VARID(ncid,'gamma_theta',vid)  !-YuP: NetCDF-f77 get vid
            call ncvpt_doubl2(ncid,vid,start,count,wkpack,istatus)

            endif   !end igrid options

            istatus = NF_CLOSE(ncid) !-YuP: NetCDF-f77

         endif

 190  continue  ! on k=1,ngen

      return
      end subroutine netcdfvec



!=YuP=> ADDED: conversion from Netcdf-2 to Netcdf-3 or higher ==========
! These routines/function convert old routines:

      integer function ncvdef2(NCID,VARNAM,VARTYP,NDIMS,VDIMS,istatus)
      !YuP: For multi-dimensional arrays (2 or more dimensions)
      ! vid=NCVDEF() is renamed to vid=ncvdef2() in *.f files
      implicit none
      INCLUDE 'netcdf.inc'
      CHARACTER*(*) VARNAM
      INTEGER istatus,vid,NCID, VARTYP,XTYPE, NDIMS,VDIMS(*)
      if (VARTYP.eq.NCFLOAT)  XTYPE=NF_FLOAT    ! 32 BITS
      if (VARTYP.eq.NCDOUBLE) XTYPE=NF_DOUBLE   ! 64 BITS
      if (VARTYP.eq.NCCHAR)   XTYPE=NF_CHAR
      if (VARTYP.eq.NCBYTE)   XTYPE=NF_BYTE
      if (VARTYP.eq.NCSHORT)  XTYPE=NF_SHORT    ! 16 BITS
      if (VARTYP.eq.NCLONG)   XTYPE=NF_INT      ! 32 BITS
      istatus = NF_DEF_VAR(NCID,VARNAM,XTYPE,NDIMS,VDIMS,vid)
      ncvdef2 = vid
      end function ncvdef2

      integer function ncvdef0(NCID,VARNAM,VARTYP,NDIMS,VDIMS0,istatus)
      !YuP: for scalars or vectors [note that VDIMS0 is scalar in this case]
      implicit none
      INCLUDE 'netcdf.inc'
      CHARACTER*(*) VARNAM
      INTEGER istatus,vid,NCID, VARTYP,XTYPE, NDIMS,VDIMS0
      if (VARTYP.eq.NCFLOAT)  XTYPE=NF_FLOAT    ! 32 BITS
      if (VARTYP.eq.NCDOUBLE) XTYPE=NF_DOUBLE   ! 64 BITS
      if (VARTYP.eq.NCCHAR)   XTYPE=NF_CHAR
      if (VARTYP.eq.NCBYTE)   XTYPE=NF_BYTE
      if (VARTYP.eq.NCSHORT)  XTYPE=NF_SHORT    ! 16 BITS
      if (VARTYP.eq.NCLONG)   XTYPE=NF_INT      ! 32 BITS
      istatus = NF_DEF_VAR(NCID,VARNAM,XTYPE,NDIMS,VDIMS0,vid)
      ncvdef0 = vid
      end function ncvdef0

      subroutine ncaptc2(NCID, vid, NAME, ATTYPE, LEN, TEXT, istatus)
      implicit none
      INCLUDE 'netcdf.inc'
      CHARACTER*(*) NAME
      CHARACTER*(*) TEXT
      INTEGER istatus,vid,NCID, LEN, ATTYPE
      istatus = NF_PUT_ATT_TEXT(NCID, vid, NAME, LEN, TEXT)
      return
      end subroutine ncaptc2

      subroutine ncvptc0(NCID, vid, START, COUNTS, TEXT, LEN, istatus)
      !YuP: START and COUNT as scalars
      implicit none
      INCLUDE 'netcdf.inc'
      CHARACTER*(*) TEXT
      integer, intent(in)  :: vid,NCID,LEN
      integer, intent(in)  :: START, COUNTS !START(*), COUNTS(*)
      INTEGER, intent(out) :: istatus
      !istatus = NF_PUT_VARA_TEXT(NCID, vid, START, COUNTS, TEXT)
      istatus = NF_PUT_VAR_TEXT(NCID, vid, TEXT)
      return
      end subroutine ncvptc0

      subroutine ncvptc2(NCID, vid, START, COUNTS, TEXT, LEN, istatus)
      !YuP: START and COUNT as vectors
      implicit none
      INCLUDE 'netcdf.inc'
      CHARACTER*(*) TEXT
      integer, intent(in)  :: vid,NCID,LEN
      integer, intent(in)  :: START(*), COUNTS(*)
      INTEGER, intent(out) :: istatus
      !istatus = NF_PUT_VARA_TEXT(NCID, vid, START, COUNTS, TEXT)
      istatus = NF_PUT_VAR_TEXT(NCID, vid, TEXT)
      return
      end subroutine ncvptc2

      subroutine ncvpt_doubl0(NCID, vid, START, COUNTS,  DVALS, istatus)
      !YuP: START, COUNTS are scalars
      implicit none
      INCLUDE 'netcdf.inc'
      INTEGER istatus,vid,NCID, START, COUNTS
      REAL*8 DVALS(*)
      istatus=NF_PUT_VARA_DOUBLE(NCID, vid, START, COUNTS, DVALS)
      return
      end subroutine ncvpt_doubl0

      subroutine ncvpt_doubl2(NCID, vid, START, COUNTS,  DVALS, istatus)
      !YuP: START, COUNTS are vectors
      implicit none
      INCLUDE 'netcdf.inc'
      INTEGER istatus,vid,NCID, START(*), COUNTS(*)
      REAL*8 DVALS(*)
      istatus=NF_PUT_VARA_DOUBLE(NCID, vid, START, COUNTS, DVALS)
      return
      end subroutine ncvpt_doubl2

      subroutine ncvpt_int0(NCID, vid, START, COUNTS,  IVALS, istatus)
      !YuP: START and COUNT as scalars 
      implicit none
      INCLUDE 'netcdf.inc'
      INTEGER istatus,vid,NCID, START, COUNTS
      INTEGER, dimension(COUNTS-START+1) :: IVALS
      !istatus=NF_PUT_VARA_INT(NCID, vid, START, COUNTS, IVALS)
      istatus=NF_PUT_VAR_INT(NCID, vid, IVALS)
      return
      end subroutine ncvpt_int0

      subroutine ncvpt_int2(NCID, vid, START, COUNTS,  IVALS, istatus)
      !YuP: START and COUNT as vectors 
      implicit none
      INCLUDE 'netcdf.inc'
      INTEGER istatus,vid,NCID, START(*), COUNTS(*)
      INTEGER IVALS(*)
      istatus=NF_PUT_VARA_INT(NCID, vid, START, COUNTS, IVALS)
      return
      end subroutine ncvpt_int2


!      subroutine NCPOPT
       !!!-YuP: call NCPOPT() commented in files
       ! No similar routine in NetCDF-3


!     OTHER NF_PUT_*** routines available in NetCDF-3 or higher:

!     VAR entire variable access:
!     INTEGER FUNCTION  NF_PUT_VAR_DOUBLE(NCID, VARID, DVAL)

!     VAR1 single value access:
!     INTEGER FUNCTION  NF_PUT_VAR1_INT(NCID, VARID, INDEX, IVAL)

!     VARA array or array section access:
!     INTEGER FUNCTION  NF_PUT_VARA_INT(NCID, VARID, START, COUNT, IVALS)

!     VARS strided access to a subsample of values:
!     INTEGER FUNCTION  NF_PUT_VARS_INT(NCID, VARID, START, COUNT, STRIDE, IVALS)

!     VARM mapped access to values not contiguous in memory:
!     INTEGER FUNCTION  NF_PUT_VARM_INT(NCID, VARID, START, COUNT, STRIDE, IMAP, IVALS)


!======================================================================
!======================================================================

      subroutine f4dwrite
      use bcast_mod, only : bcast
      use cqlcomm_mod
      use param_mod
      use tdfinterp_mod, only : tdfinterp
      use zcunix_mod, only : terp1
      use zcunix_mod, only : terp2
      implicit none
      integer :: j,k,istat1,istat2,istat3,istat4,istat5,istat6,istat7 !local
      integer :: ir,iz,iv,it ! local
      real*8 dr,dz_,dv,dt,avg ! local scalars
      real*8 rr,zz,ppsi ! local scalars, input for terp2(),terp1()
      real*8 arg1,arg2 ! local scalars, input for atan2()
      real*8 vn,pitch,rhoin,polang ! local scalars, input for tdfinterp()
      save

!MPIINSERT_INCLUDE

!     Write 4D distribution, f(v,pitch,R,Z) to netCDF file.
!     First use of this file is in the FIDA syntheic diagnostic.
!     Work with Deyong Liu, Bill Heidbrink.
!     BH120328

!     nr_f4d,nz_f4d are R,Z grid dimensions set in namelist
!     nv_f4d,nt_f4d are normalized vel and pitch angle grids.
!     Grids, at this time, are taken to be equispaced.

!     The distribution function f4d, and the grids will be allocated
!     according to the namelist input, and will be deallocated at
!     the end of the subroutine.

      real*8,dimension(:),allocatable:: f4dr,f4dz,f4dv,f4dt,f4ddv,f4ddt
      real*8,dimension(:,:,:,:),allocatable:: f4d

!MPIINSERT_IF_RANK_NE_0_RETURN

      WRITE(*,*)
      WRITE(*,*)'Entering subroutine f4dwrite'

!     Only set up for one (the first) general species:
      if (ngen.gt.1) then
         WRITE(*,*)
         WRITE(*,*)'WARNING: Output f4d dist ftn for 1st gen spec only'
         WRITE(*,*)
         k=1 ! First general species only, for now
      endif

!.......................................................................
!     Allocate temporary storage, to be deallocated at end of subroutine
!.......................................................................
      allocate(f4dr(nr_f4d),STAT=istat1)
      allocate(f4dz(nz_f4d),STAT=istat2)
      allocate(f4dv(nv_f4d),STAT=istat3)
      allocate(f4dt(nt_f4d),STAT=istat4)
      allocate(f4ddv(nv_f4d),STAT=istat5)
      allocate(f4ddt(nt_f4d),STAT=istat6)
      allocate(f4d(nr_f4d,nz_f4d,nv_f4d,nt_f4d),STAT=istat7)
      if (istat1.ne.0 .or. istat2.ne.0 .or. istat3.ne.0 .or. &
          istat4.ne.0 .or. istat5.ne.0 .or. &
          istat6.ne.0 .or. istat7.ne.0) then
         WRITE(*,*)'f4dwrite: allocation problem'
         STOP
      endif

!.......................................................................
!     Form grids:  R,Z grids fitting plasma; v in [0,1], pitch in [0,pi]
!     See eqrhopsi of calc of rmaxcon,rmincon,zmaxcon,zmincon, slightly
!     inside the LCFS.
!     Future FOW calcs can move grid out to input the chamber wall.
!.......................................................................

      dr=(rmaxcon-rmincon)/(nr_f4d-1)
      dz_=(zmaxcon-zmincon)/(nz_f4d-1)
      dv=one/(nv_f4d-1)
      dt=pi/(nt_f4d-1)

      do ir=1,nr_f4d
         f4dr(ir)=rmincon+(ir-1)*dr
      enddo
      do iz=1,nz_f4d
         f4dz(iz)=zmincon+(iz-1)*dz_
      enddo
      do iv=1,nv_f4d
         f4dv(iv)=(iv-1)*dv
      enddo
      do it=1,nt_f4d
         f4dt(it)=(it-1)*dt
      enddo

      call bcast(f4d,zero,nr_f4d*nz_f4d*nv_f4d*nt_f4d)

!.......................................................................
!     Form volume elements in velocity space
!.......................................................................

      f4ddv(1)=f4dv(2)**3/24.
      do iv=2,nv_f4d
         f4ddv(iv)=f4dv(iv)**2*dv  !v=0 vol element.....
      enddo
      f4ddt(1)=0.25*pi*f4dt(2)**2
      f4ddt(nt_f4d)=f4ddt(1)
      do it=2,nt_f4d-1
         f4ddt(it)=twopi*sin(f4dt(it))*dt
      enddo


!.......................................................................
!     For each R,Z, find rhoin,polang for use with tdfinterp.
!     Coding is similar to that in tdnpa0.
!     Only set up for eqsource="eqdsk".  Easy to generalize.
!.......................................................................

      if (eqsource.ne."eqdsk") then
         WRITE(*,*)
         WRITE(*,*) 'f4dwrite only set up for eqsource=eqdsk'
         WRITE(*,*) 'Easy to generalize, if useful'
         STOP
      endif


      if(fow.eq.'disabled')then

!        Change sign of eqpsi, to get asceding order necessary
!        for coeff1.  (Remember this when using the spline coeffs.)
         do j=1,nconteqn
            eqpsi(j)=-eqpsi(j)
         enddo
         itab(1)=1 ! can check rho value in tab(1) with debugger.
         itab(2)=0
         itab(3)=0
         do ir=1,nr_f4d
            rr=f4dr(ir)
            do iz=1,nz_f4d
               zz=f4dz(iz)
               !epsi(,) has convention of being max at mag axis.
               !eqpsi is sign reversed from this, giving min a mag axis.
               !Therefore, ppsi value is sign reversed here.
               ppsi=-terp2(rr,zz,nnr,er,nnz,ez,epsi,epsirr,epsizz, &
               epsirz,nnra,0,0)
               !(Check sign of ppsi and eqpsi with debugger)
               call terp1(nconteqn,eqpsi,eqrho,d2eqrho,ppsi,1,tab,itab)
               rhoin=tab(1)/rhomax
!              If rhoin.gt.1, point is outside the LCFS.  Leave f4d=0.
               !  Poloidal angle (rmag is major radius of magnetic axis).
               if (rhoin.lt.one) then
                  arg1=zz
                  arg2=rr-rmag
                  if(arg1.eq.zero .and. arg2.eq.zero)  then
                     polang=0.
                  else
                     polang=atan2(arg1,arg2)
                     if (polang.lt.zero) polang=polang+twopi
                  endif
                !WRITE(*,'(a,3e12.3)')'f4dwrite : rr,zz,rhoin=',rr,zz,rhoin
                  do iv=1,nv_f4d
                     vn=f4dv(iv)*vnorm
                     do it=1,nt_f4d
                        pitch=f4dt(it)
                         call tdfinterp(vn,pitch,rhoin,polang, &
                                        f4d(ir,iz,iv,it))

                     enddo
                  enddo
               endif  ! On rhoin
            enddo ! iz
            !pause
         enddo ! ir

      elseif(fow.eq.'hybrid' .or. fow.eq.'full')then
! Not applicable in this version
!         do ir=1,nr_f4d
!            rloc=f4dr(ir)  ! local R
!         do iz=1,nz_f4d
!            zloc=f4dz(iz)  ! local Z
!            do iv=1,nv_f4d
!               vn=f4dv(iv)*vnorm ! local v
!            do it=1,nt_f4d
!               pitch=f4dt(it)    ! local pitch
!               call fow_tdfinterp(k,vn,pitch,rloc,zloc,f4d(ir,iz,iv,it))
!            enddo ! it
!            enddo ! iv
!         enddo ! iz
!         enddo ! ir
      else
         STOP 'f4dwrite: only setup for fow=disabled'
      endif ! fow




!.......................................................................

!  Some printing:
!$$$      do ir=1,nr_f4d
!$$$         iz=33
!$$$         WRITE(*,*) 'f4d(ir,iz,1,:), ir,iz=', ir,iz
!$$$         WRITE(*,100) (f4d(ir,iz,1,i),i=1,nt_f4d)
!$$$         WRITE(*,*) 'f4d(ir,iz,2,:), ir,iz=', ir,iz
!$$$         WRITE(*,100) (f4d(ir,iz,2,i),i=1,nt_f4d)
!$$$         WRITE(*,*) 'f4d(ir,iz,:,1), ir,iz=', ir,iz
!$$$         WRITE(*,100) (f4d(ir,iz,j,1),j=1,nv_f4d)
!$$$         WRITE(*,*) 'f4d(ir,iz,:,nt_f4d), ir,iz=', ir,iz
!$$$         WRITE(*,100) (f4d(ir,iz,j,nt_f4d),j=1,nv_f4d)
!$$$         iz=40
!$$$         WRITE(*,*) 'f4d(ir,iz,1,:), ir,iz=', ir,iz
!$$$         WRITE(*,100) (f4d(ir,iz,1,i),i=1,nt_f4d)
!$$$         WRITE(*,*) 'f4d(ir,iz,2,:), ir,iz=', ir,iz
!$$$         WRITE(*,100) (f4d(ir,iz,2,i),i=1,nt_f4d)
!$$$         WRITE(*,*) 'f4d(ir,iz,:,1), ir,iz=', ir,iz
!$$$         WRITE(*,100) (f4d(ir,iz,j,1),j=1,nv_f4d)
!$$$         WRITE(*,*) 'f4d(ir,iz,:,nt_f4d), ir,iz=', ir,iz
!$$$         WRITE(*,100) (f4d(ir,iz,j,nt_f4d),j=1,nv_f4d)
!$$$      enddo  ! On ir
!$$$ 100  format(10ES10.2)


!  Fix some odd values:  Need to check tdfinterp
      do ir=1,nr_f4d
         do iz=1,nz_f4d
            iv=1
            avg=zero
            do it=1,nt_f4d
               avg=avg+f4d(ir,iz,iv+1,it)
            enddo  ! On it
            avg=avg/nt_f4d
            do it=1,nt_f4d
               f4d(ir,iz,iv,it)=avg
            enddo  ! On it
            iv=nv_f4d
            it=nt_f4d
            f4d(ir,iz,iv,it)=f4d(ir,iz,iv,it-1)  !Might not be needed
         enddo  ! On iz
      enddo  ! On ir

!$$$      ir=20
!$$$      iz=33
!$$$      it=nt_f4d
!$$$      WRITE(*,*)'f4d(ir,iz,:,it-1):',(f4d(ir,iz,j,it-1),j=1,nv_f4d)
!$$$      WRITE(*,*)'f4d(ir,iz,:,it):',(f4d(ir,iz,j,it),j=1,nv_f4d)

!.......................................................................



!     Change back eqpsi sign.
      do j=1,nconteqn
         eqpsi(j)=-eqpsi(j)
      enddo

      if (eqsym.eq."none") then
         WRITE(*,*)
         WRITE(*,*)'************************************************'
         WRITE(*,*)'WARNING:  NPA CALC ASSUMING UP-DOWN SYMM'
         WRITE(*,*)'          and eqsym=none.  Need to check coding'
         WRITE(*,*)'************************************************'
         WRITE(*,*)
      endif

!.......................................................................
!     Write data into netcdf file f4d.nc
!.......................................................................

      call ncwritef4d(f4dr,f4dz,f4dv,f4dt,f4ddv,f4ddt,f4d)

      deallocate(f4dr,f4dz,f4dv,f4dt,f4ddv,f4ddt,f4d)

      return

      end subroutine f4dwrite

!======================================================================
!======================================================================

      subroutine ncwritef4d(f4dr,f4dz,f4dv,f4dt,f4ddv,f4ddt,f4d)
      use param_mod
      use cqlcomm_mod
      implicit none

! --- include file for netCDF declarations
! --- (obtained from NetCDF distribution)
      include 'netcdf.inc'
!MPIINSERT_INCLUDE

!     nr_f4d,nz_f4d are R,Z grid dimensions set in namelist
!     nv_f4d,nt_f4d are dims of normalized vel and of pitch angle grids.
      real*8, dimension(:) :: f4dr(nr_f4d),f4dz(nz_f4d), &
                              f4dv(nv_f4d),f4dt(nt_f4d), &
                              f4ddv(nv_f4d),f4ddt(nt_f4d)
      real*8, dimension(nr_f4d,nz_f4d,nv_f4d,nt_f4d) :: f4d  !dims in comm.h



      integer ncid,vid,istatus
      integer chardim,char64dim
      integer dim_nr_f4d,dim_nz_f4d,dim_nv_f4d,dim_nt_f4d
      integer dims(4),start(4),count(4)

      character*128 ltitle
      
      integer ll ! local

      data start/1,1,1,1/

!     another way to do this, when the code is a mess.
      character(len=8) :: cqlpmod
      integer :: lrindx(0:lrorsa)
      integer :: lrz
      integer :: lrzmax
      character(len=256) :: mnemonic
      cqlpmod = setup0%cqlpmod
      lrindx = setup0%lrindx
      lrz = setup0%lrz
      lrzmax =  setup0%lrzmax
      mnemonic = setup0%mnemonic

!MPIINSERT_IF_RANK_NE_0_RETURN

      count(1)=nr_f4d
      count(2)=nz_f4d
      count(3)=nv_f4d
      count(4)=nt_f4d

!     Create netCDF file
      write(t_,1000) mnemonic(1:length_char(mnemonic))
 1000 format(a,"_f4d.nc")
      istatus = NF_CREATE(t_, NF_CLOBBER, ncid)
      call check_err(istatus)

!     Define dimensions
      istatus= NF_DEF_DIM(ncid, 'dim_nr_f4d', nr_f4d, dim_nr_f4d)
      istatus= NF_DEF_DIM(ncid, 'dim_nz_f4d', nz_f4d, dim_nz_f4d)
      istatus= NF_DEF_DIM(ncid, 'dim_nv_f4d', nv_f4d, dim_nv_f4d)
      istatus= NF_DEF_DIM(ncid, 'dim_nt_f4d', nt_f4d, dim_nt_f4d)
      istatus= NF_DEF_DIM(ncid, 'chardim',      8,    chardim)
      istatus= NF_DEF_DIM(ncid, 'char64dim',   64,    char64dim)

!     Define vector of dimensions
      dims(1)=dim_nr_f4d
      dims(2)=dim_nz_f4d
      dims(3)=dim_nv_f4d
      dims(4)=dim_nt_f4d

!     Define netCDF variables

      ltitle='CQL3D 4D Distn Function (R,Z,v,pitch): '//version
      if( length_char(ltitle).gt.128 ) stop 'Adjust ltitle in f4dwrite'
      call ncaptc2(ncid,NCGLOBAL,'title',NCCHAR,length_char(ltitle), &
           ltitle,istatus)

      vid=ncvdef0(ncid,'version',NCCHAR,1,char64dim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,20, &
                'CQL3D version number',istatus)

      vid=ncvdef0(ncid,'mnemonic',NCCHAR,1,char64dim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,23, &
                'Mnemonic run identifier',istatus)

      vid=ncvdef0(ncid,'vnorm',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,33, &
                'velocity (momentum-per-mass) norm',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7, &
                           'cms/sec',istatus)

      vid=ncvdef0(ncid,'enorm',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,20, &
                           'Energy normalization',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                           'keV',istatus)

      vid=ncvdef0(ncid,'rmag',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,27, &
                'Magnetic axis major radius',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                 'cms',istatus)

      vid=ncvdef0(ncid,'zmag',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,31, &
                'Magnetic axis vertical position',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                 'cms',istatus)

      vid=ncvdef0(ncid,'eqsym',NCCHAR,1,chardim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,28, &
                'Indicator for symmetrization',istatus)

      vid=ncvdef0(ncid,'zshift',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,39, &
                'Vertical shift of equilibrium per eqsym',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                'cms',istatus)

      vid=ncvdef0(ncid,'rmincon',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,37, &
                'Minimum major radius at outer contour',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                'cms',istatus)

      vid=ncvdef0(ncid,'rmaxcon',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,37, &
                'Maximum major radius at outer contour',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                 'cms',istatus)
      call check_err(istatus)

      vid=ncvdef0(ncid,'zmincon',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,39, &
                'Minimum axial (z-dirn) at outer contour',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                'cms',istatus)
      call check_err(istatus)

      vid=ncvdef0(ncid,'zmaxcon',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,39, &
                'Maximum axial (Z-dirn) at outer contour',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                'cms',istatus)
      call check_err(istatus)

      vid=ncvdef0(ncid,'nr_f4d',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,40, &
                'Major radius grid dimension(=dim_nr_f4d)',istatus)
      call check_err(istatus)

      vid=ncvdef0(ncid,'nz_f4d',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,35, &
                'Axial Z-grid dimension(=dim_nz_f4d)',istatus)
      call check_err(istatus)

      vid=ncvdef0(ncid,'nv_f4d',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,42, &
                'Normalized vel grid dimension(=dim_nv_f4d)',istatus)
      call check_err(istatus)

      vid=ncvdef0(ncid,'nt_f4d',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,39, &
                'Pitch angle grid dimension(=dim_nt_f4d)',istatus)
      call check_err(istatus)

      vid=ncvdef0(ncid,'f4dr',NCDOUBLE,(1),dim_nr_f4d,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,40, &
                'Major Radius mesh (presently equispaced)',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,'cms',istatus)
      call check_err(istatus)

      vid=ncvdef0(ncid,'f4dz',NCDOUBLE,(1),dim_nz_f4d,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,35, &
                'Axial Z-mesh (presently equispaced)',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3,'cms',istatus)
      call check_err(istatus)

      vid=ncvdef0(ncid,'f4dv',NCDOUBLE,(1),dim_nv_f4d,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,41, &
                'Normalized to vnorm vel mesh (equispaced)',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,8,'unitless',istatus)
      call check_err(istatus)

      vid=ncvdef0(ncid,'f4dt',NCDOUBLE,(1),dim_nt_f4d,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,39, &
                'Pitch angle mesh (presently equispaced)',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7,'radians',istatus)
      call check_err(istatus)

      vid=ncvdef0(ncid,'f4ddv',NCDOUBLE,(1),dim_nv_f4d,istatus)
      call ncaptc2(ncid,vid,'long_name1',NCCHAR,7, &
                 'v**2*dv',istatus)
      call ncaptc2(ncid,vid,'long_name2',NCCHAR,29, &
                'Normalized vel volume element',istatus)

      vid=ncvdef0(ncid,'f4ddt',NCDOUBLE,(1),dim_nt_f4d,istatus)
      call ncaptc2(ncid,vid,'long_name1',NCCHAR,22, &
                '2*pi*sin(pitch)*dpitch',istatus)
      call ncaptc2(ncid,vid,'long_name2',NCCHAR,36, &
                'vel space pitch angle volume element',istatus)

      vid=ncvdef2(ncid,'f4d',NCDOUBLE,4,dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,41, &
                'Distribution function on R,Z,v,pitch grid',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,28, &
                'vnorm**3/(cm**3*(cm/sec)**3)',istatus)
      call check_err(istatus)

!     End the define-mode
      istatus= NF_ENDDEF(ncid)
      call check_err(istatus)


!.......................................................................
!     Put data in the file
!.......................................................................

      istatus= NF_INQ_VARID(ncid,'version',vid)
      ll=length_char(version)
      call ncvptc0(ncid,vid,1,ll,version,ll,istatus)

      istatus= NF_INQ_VARID(ncid,'mnemonic',vid)
      ll=length_char(mnemonic)
      call ncvptc0(ncid,vid,1,ll,mnemonic,ll,istatus)

      istatus= NF_INQ_VARID(ncid,'vnorm',vid)
      call ncvpt_doubl0(ncid,vid,(1),1,(vnorm),istatus)

      istatus= NF_INQ_VARID(ncid,'enorm',vid)
      call ncvpt_doubl0(ncid,vid,(1),1,(enorm),istatus)

      istatus= NF_INQ_VARID(ncid,'rmag',vid)
      call ncvpt_doubl0(ncid,vid,(1),1,(rmag),istatus)

      istatus= NF_INQ_VARID(ncid,'zmag',vid)
      call ncvpt_doubl0(ncid,vid,(1),1,(zmag),istatus)

      istatus= NF_INQ_VARID(ncid,'eqsym',vid)
      call ncvptc0(ncid,vid,1,8,eqsym,8,istatus)

      istatus= NF_INQ_VARID(ncid,'zshift',vid)
      call ncvpt_doubl0(ncid,vid,(1),1,(zshift),istatus)

      istatus= NF_INQ_VARID(ncid,'rmincon',vid)
      call ncvpt_doubl0(ncid,vid,(1),1,(rmincon),istatus)

      istatus= NF_INQ_VARID(ncid,'rmaxcon',vid)
      call ncvpt_doubl0(ncid,vid,(1),1,(rmaxcon),istatus)

      istatus= NF_INQ_VARID(ncid,'zmincon',vid)
      call ncvpt_doubl0(ncid,vid,(1),1,(zmincon),istatus)

      istatus= NF_INQ_VARID(ncid,'zmaxcon',vid)
      call ncvpt_doubl0(ncid,vid,(1),1,(zmaxcon),istatus)

      istatus= NF_INQ_VARID(ncid,'nr_f4d',vid)
      call ncvpt_int0(ncid,vid,1,1,(nr_f4d),istatus)

      istatus= NF_INQ_VARID(ncid,'nz_f4d',vid)
      call ncvpt_int0(ncid,vid,1,1,(nz_f4d),istatus)

      istatus= NF_INQ_VARID(ncid,'nv_f4d',vid)
      call ncvpt_int0(ncid,vid,1,1,(nv_f4d),istatus)

      istatus= NF_INQ_VARID(ncid,'nt_f4d',vid)
      call ncvpt_int0(ncid,vid,1,1,(nt_f4d),istatus)

      istatus= NF_INQ_VARID(ncid,'f4dr',vid)
      call ncvpt_doubl0(ncid,vid,(1),nr_f4d,f4dr,istatus)

      istatus= NF_INQ_VARID(ncid,'f4dz',vid)
      call ncvpt_doubl0(ncid,vid,(1),nz_f4d,f4dz,istatus)

      istatus= NF_INQ_VARID(ncid,'f4dv',vid)
      call ncvpt_doubl0(ncid,vid,(1),nv_f4d,f4dv,istatus)

      istatus= NF_INQ_VARID(ncid,'f4dt',vid)
      call ncvpt_doubl0(ncid,vid,(1),nt_f4d,f4dt,istatus)

      istatus= NF_INQ_VARID(ncid,'f4ddv',vid)
      call ncvpt_doubl0(ncid,vid,(1),nv_f4d,f4ddv,istatus)

      istatus= NF_INQ_VARID(ncid,'f4ddt',vid)
      call ncvpt_doubl0(ncid,vid,(1),nt_f4d,f4ddt,istatus)

      istatus= NF_INQ_VARID(ncid,'f4d',vid)
      call ncvpt_doubl2(ncid,vid,start,count,f4d,istatus)

      return

      end subroutine ncwritef4d


!-----------------------------------------------------------------------
end module netcdfrw2_mod