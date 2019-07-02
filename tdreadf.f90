module tdreadf_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use bcast_mod, only : ibcast
  use cqlcomm_mod
  use cqlconf_mod, only : read_all_conf_nml
  !XXXXuse pack21_mod, only : unpack21
  use param_mod
  use tdnflxs_mod, only : tdnflxs
  use zcunix_mod, only : allocate_error
  use zcunix_mod, only : coeff1
  use zcunix_mod, only : terp1

  !---END USE


  external unpack21  !XXXXX

!
!
!      module tdreadf_mod
!      use iso_c_binding.only :: c_double

!      !common/tdreadf_com/
!      real(c_double), private :: x_rstrt,v_rstrt,v_rstrt2,y_rstrt, &
!          f_rstrt_ln,v,v2,aa,bb,d2fparr,workk,tam2r,cint2r, &
!          f_rstrt,jf_rstrt
!      save
!      contains


contains

  subroutine tdreadf(kopt)
    use cqlconf_mod, only : setup0
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!.......................................................................
!     Read namelist, distribution function and spatial source from
!       a previous run.
!     kopt = 1: Read namelists from text file distrfunc (should be first
!               call, to place file pointer at beginning of f).
!               Only called if setup0%nlrestrt="enabled", otherwise distrfunc
!               not needed.
!          = 2: Read f
!.......................................................................


      character*1 blank
      character*8 ilrestrt
      include 'frname_decl.h'
      include 'frname.h'
!
! --- include file for netCDF declarations
! --- (obtained from NetCDF distribution)
      include 'netcdf.inc'

!MPIINSERT_INCLUDE

      integer ncid,istatus
      integer xdim,ydim,rdim,gen_species_dim,vid, dimid, dimlen
!BH180517      integer count(3),start(3)
      integer count(4),start(4)
      integer tdim,r00dim,tdim_rstrt,r00dim_rstrt
      integer r00_count(2)
      character*128 name

!     Pointers for dynamic memory allocation, local usage:
      real(c_double), dimension(:), pointer :: x_rstrt,v_rstrt,v_rstrt2,y_rstrt
      real(c_double), dimension(:), pointer :: f_rstrt_ln,v,v2,aa,bb,d2fparr
      real(c_double), dimension(:), pointer :: workk,tam2r,cint2r
      real(c_double), dimension(:,:,:), pointer :: f_rstrt
      integer, dimension(:), pointer :: jf_rstrt
      real(c_double), allocatable :: wkpack(:) !local working array for pack21 YuP[2019]
      real(c_double), allocatable :: temp_rstrt(:,:) !(0:iy_rstrt+1,0:jx_rstrt+1)
!      save foverf,renorm_f,k,iy_rstrt,jx_rstrt,lrz_rstrt

!.......................................................................

      iunwrif=19

      if (kopt .eq. 2) go to 200  !line 109

      open(unit=iunwrif,file='distrfunc')
      rewind(iunwrif)

!.......................................................................
!l    1. Read number of lines to skip and skip them
!.......................................................................

      read(iunwrif,*) inbline
!MPIINSERT_IF_RANK_EQ_0
      WRITE(*,*)'tdreadf:  inbline= ',inbline
!MPIINSERT_ENDIF_RANK
      do 102 l=1,inbline
        read(iunwrif,'(a)') blank
 102  continue

!.......................................................................
!l    1.2 Read namelists   [Originally: No read frsetup?]
!         It is possible that namelist entries are to be changed
!         from cqlinput namelist file, for example, modified
!         plasma profiles.
!.......................................................................

!BH080204:  Will work towards putting all necessary restart data
!BH080204:  into the netcdf file, for netcdf restart case.
      if (setup0%nlrestrt.ne."ncdfdist" .and. setup0%nlrestrt.ne."ncregrid") then !BH080204

         ilrestrt=setup0%nlrestrt
         call read_all_conf_nml(iunwrif)
         !NMLXXX, frsetup nml not converted yet
         read(iunwrif,frsetup)
         setup0%nlrestrt=ilrestrt

         if ((lrors.ne.setup0%lrz .and. setup0%cqlpmod.ne."enabled") .or. &
              (lrors.ne.setup0%ls  .and. setup0%cqlpmod.eq."enabled")) then
            !MPIINSERT_IF_RANK_EQ_0
            PRINT *,' current lrors=',lrors,' does not correspond to old ', &
                 'setup0%lrz=',setup0%lrz,' or setup0%ls=',setup0%ls
            !MPIINSERT_ENDIF_RANK
            stop 'tdreadf'
         endif

      endif ! on setup0%nlrestrt

      return

!.......................................................................
!l    2. Read f(i,j,k,l) from text distfunc file.
!     l_.eq.lrors clause detects first of l_=lrors,1,-1 calls
!.......................................................................

 200  continue ! kopt .eq. 2   handle
      ! Note: tdreadf(2) is called for each l_ , which is in comm.h

      if (setup0%nlrestrt.eq."enabled" .and. l_.eq.lrors) then  !elseif, ln 145

!MPIINSERT_IF_RANK_EQ_0
         WRITE(*,*)'tdreadf:  Reading distfunc (text f)'
!MPIINSERT_ENDIF_RANK
         do 210 k=1,ngen
            do 211 il=1,lrors
               do 212 j=1,jx
                  read(iunwrif,9210) (f(i,j,k,il),i=1,iy_(il))
 212           continue
 211        continue
 210     continue  !On k

!.......................................................................
!l    2.2 Read spasou(i,j,k,l);  only works for setup0%nlrestrt.eq.'enabled'
!.......................................................................

         if (setup0%cqlpmod.eq."enabled") then
            do 220 k=1,ngen
               do 221 il=1,lrors
                  do 222 j=1,jx
                     read(iunwrif,9210) (spasou(i,j,k,il),i=1,iy_(il))
 222              continue
 221           continue
 220        continue
         endif

         close(unit=iunwrif)




!.......................................................................
!     Read distribution function from NetCDF file, v-mesh UNchanged
!.......................................................................

      elseif (setup0%nlrestrt.eq."ncdfdist" .and. l_.eq.lrors) then
                                                       !elseif, line 290

         !YuP[2019-02-07] added wkpack and temp_rstrt
         if (.NOT.ALLOCATED(wkpack)) then !allocate working array for pack21
           nwkpack=(iy+2)*(jx+2)
           allocate(wkpack(nwkpack),STAT=istat)
           if(istat.ne.0) &
              call allocate_error("wkpack, sub tdreadf",0,istat)
           wkpack=zero
         endif
         if (.NOT.ALLOCATED(temp_rstrt)) then !allocate working array for pack21
           allocate(temp_rstrt(0:iy+1,0:jx+1),STAT=istat)
           if(istat.ne.0) &
              call allocate_error("temp_rstrt, sub tdreadf",0,istat)
           temp_rstrt=zero
         endif
!$$$BH180517:
!$$$         if (ngen.gt.1) then
!$$$MPIINSERT_IF_RANK_EQ_0
!$$$            WRITE(*,*)'tdreadf: need ngen.gt.1 reading of distfunc.nc'
!$$$MPIINSERT_ENDIF_RANK
!$$$            stop 'in tdreadf'
!$$$         endif

!     Open existing netCDF file [Typically, the distrfunc.nc is linked
!     to the setup0%mnemonic.nc file from an earlier cql3d run.]
!MPIINSERT_IF_RANK_EQ_0
         WRITE(*,*) &
         'tdreadf[setup0%nlrestrt="ncdfdist"]: Reading distrfunc.nc (netcdf)'
!MPIINSERT_ENDIF_RANK

         istatus = NF_OPEN('distrfunc.nc', 0, ncid)
         if (istatus .NE. NF_NOERR) then
            t_=trim(setup0%mnemonic(1:length_char(setup0%mnemonic)))//'.nc' !YuP added trim
!MPIINSERT_IF_RANK_EQ_0
            WRITE(*,*)'tdreadf: distrfunc.nc missing/problem'
            WRITE(*,*)'tdreadf: CHECKING for setup0%mnemonic file : ',t_
            WRITE(*,*)'tdreadf: BUT BE CAREFUL, the file will be '
            WRITE(*,*)'tdreadf:     overwritten at end of present run. '
!MPIINSERT_ENDIF_RANK
            ! YuP[03-01-2016] Added: if distrfunc.nc is missing,
            ! try reading the setup0%mnemonic.nc file.
            ! But be careful: the file will be overwritten by the end
            ! of present run !
            istatus = NF_OPEN(t_, 0, ncid)
            if (istatus .NE. NF_NOERR) then
!MPIINSERT_IF_RANK_EQ_0
               WRITE(*,*)
               WRITE(*,*)'tdreadf: Cannot open/missing setup0%mnemonic.nc'// &
                    ' for which opening has been attempted in lieu'// &
                    ' of distrfunc.nc'
               WRITE(*,*)
!MPIINSERT_ENDIF_RANK
            else ! found the file
!MPIINSERT_IF_RANK_EQ_0
               WRITE(*,*)'tdreadf: found setup0%mnemonic file: ',t_
!MPIINSERT_ENDIF_RANK
            endif
         endif


!     read in dimension IDs and sizes
         istatus= NF_INQ_DIMID(ncid,'xdim',xdim) !-YuP: NetCDF-f77
         istatus= NF_INQ_DIMID(ncid,'ydim',ydim) !-YuP: NetCDF-f77
         istatus= NF_INQ_DIMID(ncid,'rdim',rdim) !-YuP: NetCDF-f77
         istatus= NF_INQ_DIMID(ncid,'gen_species_dim',gen_species_dim)

!     --- inquire about dimension sizes ---
         istatus= NF_INQ_DIM(ncid,ydim,name,iy_rstrt)  !-YuP: NetCDF-f77
         istatus= NF_INQ_DIM(ncid,xdim,name,jx_rstrt)  !-YuP: NetCDF-f77
         istatus= NF_INQ_DIM(ncid,rdim,name,lrz_rstrt) !-YuP: NetCDF-f77
         istatus= NF_INQ_DIM(ncid,gen_species_dim,name,ngen_rstrt)

         if (iy_rstrt.ne.iy &
              .or. jx_rstrt.ne.jx &
              .or. lrz_rstrt.ne.setup0%lrz .or. ngen_rstrt.ne.ngen) then
!MPIINSERT_IF_RANK_EQ_0
            WRITE(*,*)'Problem with distrfunc.nc file'
            WRITE(*,*)'  iy,jx,setup0%lrz,ngen=',iy,jx,setup0%lrz,ngen, &
                 ' iy_rstrt,jx_rstrt,lrz_rstrt,ngen_rstrt=', &
                 iy_rstrt,jx_rstrt,lrz_rstrt,ngen_rstrt
!MPIINSERT_ENDIF_RANK
            STOP ' in tdreadf'
         endif

!-----pitch angle variable y
!-YuP:         vid = ncvid(ncid,'iy_',istatus)
         istatus= NF_INQ_VARID(ncid,'iy_',vid)  !-YuP: NetCDF-f77 get vid
!-YuP:         call ncvgt(ncid,vid,1,setup0%lrz,iy_,istatus)
         istatus= NF_GET_VARA_INT(ncid,vid,(1),(setup0%lrz),iy_) !-YuP: NetCDF-f77
         do ll=1,setup0%lrz
            if (iy_(ll).ne.iy) then
!MPIINSERT_IF_RANK_EQ_0
               WRITE(*,*)'tdreadf:Variable iy_(ll), but only cnst setup'
!MPIINSERT_ENDIF_RANK
               stop 'in tdtreadf'
            endif
         enddo

!-----elecfld, restore if ampfmod=enabled
         if (ampfmod.eq."enabled") then
            istatus= NF_INQ_DIMID(ncid,'tdim',tdim)
            istatus= NF_INQ_DIMID(ncid,'r00dim',r00dim)
            istatus= NF_INQ_DIM(ncid,tdim,name,tdim_rstrt)
            istatus= NF_INQ_DIM(ncid,r00dim,name,r00dim_rstrt)
            istatus= NF_INQ_VARID(ncid,'elecfld',vid)
            start(1)=1
            start(2)=tdim_rstrt
            r00_count(1)=r00_rstrt
            r00_count(2)=1
            istatus= NF_GET_VARA_DOUBLE( &
                     ncid,vid,start,r00_count,elecfld)
            elecfldc=elecfld(0)
            istatus=NF_INQ_VARID(ncid,'elecfldb',vid)
            istatus= NF_GET_VAR1_DOUBLE( &
                     ncid,vid,(1),elecfldb)
         endif

!-----distribution function f: restore after checking that
!-----there is only one time-set.  [BH181025: Could adjust code
!-----for netcdfshort multiple time f, but not much need.]

         istatus=NF_INQ_VARID(ncid,'netcdfshort',vid)
         istatus= NF_GET_VAR_TEXT(ncid,vid,netcdfshort)
         if (netcdfshort.eq."enabled" .or. &
             netcdfshort.eq."longer_f" .or. &
             netcdfshort.eq."lngshrtf") then
            WRITE(*,*)
            WRITE(*,*)'STOP: distrfunc.nc not setup for single time f'
            WRITE(*,*)
            stop
         endif


         !The structure of 'f' in setup0%mnemonic.nc file can be different
         !depending on ngen and on netcdfshort settings.
         !YuP[2019-02-07] See netcdfrw2: 'f' can be 3D, or 4D, or 5D.
         !
         !---1--- Case of saving f at each or several selected time steps,
         !      (netcdfshort.eq.'longer_f').or.(netcdfshort.eq.'lngshrtf')
         !if (ngen.eq.1) then
         ! 'f' is 4D, NCDOUBLE,4, dimsf(1:4)={ydim,xdim,rdim, tdim}
         !                        start1(3)=ll (in do loop)
         !                        start1(4)=numrec1 (or numrecsave)
         !                        count1(1:4)={iy,jx,1,1} (in ll loop)
         !else  !ngen.ge.2
         ! 'f' is 5D, NCDOUBLE,5, dimsf(1:5)={ydim,xdim,rdim, gdim,tdim}
         !                        startg(3)=ll (in do loop)
         !                        startg(4)=k  (in do loop)
         !                        startg(5)=numrec1 (or numrecsave)
         !                        countg(1:5)={iy,jx,1,1,1} (in k,ll loops)
         !endif  !on ngen
         !The above case is difficult: we need to know the value of numrec1
         !(it can be less than nstop)
         !
         !---2--- Case of saving f only at the end of run
         !        (other values of netcdfshort, e.g. 'disabled',
         !         but not 'enabled' when f is not saved at all)
         !if (ngen.eq.1) then
         ! 'f' is 3D, NCDOUBLE,3, dimsf(1:3)={ydim,xdim,rdim}
         !                        start1(3)=ll (in do loop)
         !                        count1(1:3)={iy,jx,1} (in ll loop)
         !else  !ngen.ge.2
         ! 'f' is 4D, NCDOUBLE,4, dimsf(1:4)={ydim,xdim,rdim, gdim}
         !                        startg(3)=ll (in do loop)
         !                        startg(4)=k  (in do loop)
         !                        countg(1:4)={iy,jx,1,1} (in k,ll loops)
         !endif  !on ngen
         !Note that ngen=1 and ngen>1 cases
         !could be combined during reading.
         !Just use
         !                        start(3)=ll (in do ll loop)
         !                        start(4)=k  (in do k loop)
         !                        count(1:4)={iy,jx,1,1} (in k,ll loops)
         !                        In case of ngen=1 the values of
         !                        start(4) and count(4) will not be used.
         !Note: gdim is same as gen_species_dim (=ngen)
         ! Search "temp1(i,j)=f(" in netcdfrw2.


         ! The reading of f below is only valid in case when f was recorded
         ! at the last time step only (netcdfshort='disabled').
         start(1)=1
         start(2)=1
         count(1)=iy
         count(2)=jx
         istatus= NF_INQ_VARID(ncid,'f',vid)  !-YuP: NetCDF-f77 get vid
         istatus= NF_INQ_DIMID(ncid,'f', dimid)
         istatus= NF_INQ_DIMLEN(ncid, dimid, dimlen)
!MPIINSERT_IF_RANK_EQ_0
         WRITE(*,*)'tdreaf[nlrestrt="ncdfdist"]: dimid=',dimid
         WRITE(*,*)'tdreaf[nlrestrt="ncdfdist"]: dimlen=',dimlen
!MPIINSERT_ENDIF_RANK
         call bcast(f,zero,(iy+2)*(jx+2)*ngen*lrors)  !This f is f_code

!         do k=1,ngen
!         do ll=1,setup0%lrz
!            start(3)=ll
!            start(4)=k
!            count(3)=1
!            count(4)=1
!            istatus= NF_GET_VARA_DOUBLE(ncid,vid,start,count,tem1)
!            ij=0
!            do j=1,jx
!               do i=1,iy
!                  ij=ij+1
!                  f(i,j,k,ll)=tem1(ij)
!               enddo
!            enddo
!         enddo  !on ll
!         enddo  !on k=1,ngen

         !YuP[2019-02-07] Replaced by this version
         !(the above does not work).
         !In netcdfrw2, f() was packed as
         !call pack21(temp1,0,iyp1,0,jxp1,wkpack,iy,jx)
         !for each given ll, each k,
         !and each time step (if more than one).
         !So now we unpack it in exactly same way.
         do k=1,ngen
         do ll=1,setup0%lrz
           ! These values (3-4) can be different, dep. on cases
           start(3)=ll
           start(4)=k !in case of ngen=1, start(4)=1,
           !actually start(4) and count(4) are not needed in case ngen=1
           count(3)=1
           count(4)=1
           istatus= NF_GET_VARA_DOUBLE(ncid,vid,start,count,wkpack)
           !YuP[2019-02-07] unpack wkpack(:) into temp_rstrt(i,j)
           call unpack21(temp_rstrt,0,iy+1,0,jx+1,wkpack,iy,jx)
           do j=1,jx
           do i=1,iy
              f(i,j,k,ll)=temp_rstrt(i,j)
           enddo
           enddo
         enddo  !on ll
         enddo  !on k=1,ngen


         istatus=NF_CLOSE(ncid)
!MPIINSERT_IF_RANK_EQ_0
         WRITE(*,*)'tdreaf[nlrestrt="ncdfdist"]:Done reading f'
         WRITE(*,*)'tdreaf_360: For checkup SUM(f)=', SUM(f)
         !pause
!MPIINSERT_ENDIF_RANK






!.......................................................................
!     Read netcdf distribution function, v-mesh changed
!.......................................................................

!BH110125. See cqlinput_help
!_rstrt refers to restart distn,

!In this section, read only one flux surface at a time, since
!need flux surface geometry quantities zmaxpsi(lr_) and tau(lr_)
!which are also calculated one flux surface at a time.

      elseif (setup0%nlrestrt.eq."ncregrid") then !re-grid mod option
                                           !endif, line 645

         if (l_.eq.lrors) then  !That is, first l_ call of tdreadf

!$$$BH180518
!$$$         if (ngen.gt.1) then
!$$$MPIINSERT_IF_RANK_EQ_0
!$$$            WRITE(*,*)
!$$$            WRITE(*,*)'================================================'
!$$$            WRITE(*,*)'tdreadf: Re-grid of rstrt distribution:'
!$$$            WRITE(*,*)'requires ngen.le.1'
!$$$            WRITE(*,*)'================================================'
!$$$MPIINSERT_ENDIF_RANK
!$$$            stop
!$$$         else
!$$$            k=1
!$$$         endif


!     Open existing NetCDF file
!MPIINSERT_IF_RANK_EQ_0
         WRITE(*,*) &
         'tdreadf[nlrestrt.eq."ncregrid"]: Reading distfunc.nc (netcdf)'
!MPIINSERT_ENDIF_RANK

         istatus = NF_OPEN('distrfunc.nc', 0, ncid)
         if (istatus .NE. NF_NOERR) then
!MPIINSERT_IF_RANK_EQ_0
            WRITE(*,*)'tdreadf:  distrfunc.nc missing/problem'
!MPIINSERT_ENDIF_RANK
         endif


!     read in dimension IDs and sizes
         istatus= NF_INQ_DIMID(ncid,'xdim',xdim)
         istatus= NF_INQ_DIMID(ncid,'ydim',ydim)
         istatus= NF_INQ_DIMID(ncid,'rdim',rdim)
         istatus= NF_INQ_DIMID(ncid,'gen_species_dim',gen_species_dim)

!     --- inquire about dimension sizes ---
         istatus= NF_INQ_DIM(ncid,xdim,name,jx_rstrt)
         istatus= NF_INQ_DIM(ncid,ydim,name,iy_rstrt)
         istatus= NF_INQ_DIM(ncid,rdim,name,lrz_rstrt)
         istatus= NF_INQ_DIM(ncid,gen_species_dim,name,ngen_rstrt)

!        for nlrestrt.eq.'ncregrid', jx_rstrt can be different from jx
!        Also enforce that ngen is .le. gen_species_dim.
         if (iy_rstrt.ne.iy .or. lrz_rstrt.ne.setup0%lrz .or. ngen_rstrt.ne.ngen) then
!     1        .or. jx_rstrt.ne.jx

!MPIINSERT_IF_RANK_EQ_0
            WRITE(*,*)'Problem with distrfunc.nc file'
            WRITE(*,*)'  iy,setup0%lrz,ngen=',iy,setup0%lrz,ngen, &
                 ' iy_rstrt,lrz_rstrt,ngen_rstrt=', &
                 iy_rstrt,lrz_rstrt,ngen_rstrt
!MPIINSERT_ENDIF_RANK
            STOP ' in tdreadf'
         endif

!-----pitch angle variable y
         istatus= NF_INQ_VARID(ncid,'iy_',vid)
         istatus= NF_GET_VARA_INT(ncid,vid,(1),(setup0%lrz),iy_)
         do ll=1,setup0%lrz
            if (iy_(ll).ne.iy) then
!MPIINSERT_IF_RANK_EQ_0
               WRITE(*,*)'tdreadf: Variable iy_(ll), only cnst setup'
!MPIINSERT_ENDIF_RANK
               stop 'in tdtreadf'
            endif
         enddo

!-----distribution function f(i,j,ll) [vnorm**3/(cm**3*(cm/sec)**3)]
!-----to be obtained by extrapolation and interpolation onto new
!-----x-grid, BH110121.
!-----general species: for coding simplicity, all gen species read
!                      in and restored.

!     Allocate temporary space
!         allocate (f_rstrt(iy,jx_rstrt,setup0%lrz,ngen),STAT = istat)   !Can't do this: f_rstrt has 3 dims
         allocate (f_rstrt(iy,jx_rstrt,setup0%lrz),STAT = istat)
         if(istat .ne. 0) &
              call allocate_error("f_rstrt, sub tdreadf",0,istat)
         call bcast(f_rstrt,zero,SIZE(f_rstrt))

         allocate (x_rstrt(jx_rstrt),STAT = istat)
         if(istat .ne. 0) &
              call allocate_error("x_rstrt, sub tdreadf",0,istat)
         call bcast(x_rstrt,zero,SIZE(x_rstrt))

         allocate (v_rstrt(jx_rstrt),STAT = istat)
         if(istat .ne. 0) &
              call allocate_error("v_rstrt, sub tdreadf",0,istat)
         call bcast(v_rstrt,zero,SIZE(v_rstrt))

         allocate (v_rstrt2(jx_rstrt),STAT = istat)
         if(istat .ne. 0) &
              call allocate_error("v_rstrt2, sub tdreadf",0,istat)
         call bcast(v_rstrt2,zero,SIZE(v_rstrt2))

         allocate (y_rstrt(iy),STAT = istat)  !Could just use y()
         if(istat .ne. 0) &
              call allocate_error("y_rstrt, sub tdreadf",0,istat)
         call bcast(y_rstrt,zero,SIZE(y_rstrt))

         allocate (jf_rstrt(iy),STAT = istat)
         if(istat .ne. 0) &
              call allocate_error("jf_rstrt, sub tdreadf",0,istat)
         call ibcast(jf_rstrt,0,SIZE(jf_rstrt))

         allocate (aa(iy),STAT = istat)
         if(istat .ne. 0) &
              call allocate_error("aa, sub tdreadf",0,istat)
         call bcast(aa,zero,SIZE(aa))

         allocate (bb(iy),STAT = istat)
         if(istat .ne. 0) &
              call allocate_error("bb, sub tdreadf",0,istat)
         call bcast(bb,zero,SIZE(bb))

         allocate (f_rstrt_ln(jx_rstrt),STAT = istat)
         if(istat .ne. 0) &
              call allocate_error("f_rstrt_ln, sub tdreadf",0,istat)
         call bcast(f_rstrt_ln,zero,SIZE(f_rstrt_ln))

         allocate (d2fparr(jx_rstrt),STAT = istat)
         if(istat .ne. 0) &
              call allocate_error("d2fparr, sub tdreadf",0,istat)
         call bcast(d2fparr,zero,SIZE(d2fparr))

         allocate (workk(3*jx_rstrt+1),STAT = istat)
         if(istat .ne. 0) &
              call allocate_error("workk, sub tdreadf",0,istat)
         call bcast(workk,zero,SIZE(workk))

         allocate (v(jx),STAT = istat)
         if(istat .ne. 0) &
              call allocate_error("v, sub tdreadf",0,istat)
         call bcast(v,zero,SIZE(v))

         allocate (v2(jx),STAT = istat)
         if(istat .ne. 0) &
              call allocate_error("v2, sub tdreadf",0,istat)
         call bcast(v2,zero,SIZE(v2))

         !YuP[2019-02-07] added wkpack and temp_rstrt
         if (.NOT.ALLOCATED(wkpack)) then !allocate working array for pack21
           nwkpack=(iy_rstrt+2)*(jx_rstrt+2) !and iy_rstrt=iy, see above
           allocate(wkpack(nwkpack),STAT=istat)
           if(istat.ne.0) &
              call allocate_error("wkpack, sub tdreadf",0,istat)
           wkpack=zero
         endif
         if (.NOT.ALLOCATED(temp_rstrt)) then !allocate working array for pack21
           allocate(temp_rstrt(0:iy_rstrt+1,0:jx_rstrt+1),STAT=istat)
           if(istat.ne.0) &
              call allocate_error("temp_rstrt, sub tdreadf",0,istat)
           temp_rstrt=zero
         endif

         istatus= NF_INQ_VARID(ncid,'enorm',vid)
         istatus= NF_GET_VAR1_DOUBLE(ncid,vid,(1),enorm_rstrt)

         istatus= NF_INQ_VARID(ncid,'vnorm',vid)
         istatus= NF_GET_VAR1_DOUBLE(ncid,vid,(1),vnorm_rstrt)

         istatus= NF_INQ_VARID(ncid,'x',vid)
         istatus= NF_GET_VARA_DOUBLE(ncid,vid,(1),(jx_rstrt),x_rstrt)

!MPIINSERT_IF_RANK_EQ_0
         WRITE(*,*)'tdreadf: re-grid, vnorm, vnorm_rstrt =', &
              vnorm,vnorm_rstrt
!MPIINSERT_ENDIF_RANK

         vnorm_rstrt2=vnorm_rstrt*vnorm_rstrt !Restart grid, momntm/mass
         do j=1,jx_rstrt
            v_rstrt(j)=vnorm_rstrt*x_rstrt(j)
            v_rstrt2(j)=vnorm_rstrt2*x_rstrt(j)*x_rstrt(j)
         enddo

         do j=1,jx
            v(j)=x(j)*vnorm  !Code grid momentum-per-mass
            v2(j)=v(j)*v(j)
         enddo

         call bcast(f,zero,(iy+2)*(jx+2)*ngen*lrors)

!        Allocate space and read data for density calc below.
         allocate (tam2r(jx_rstrt),STAT = istat)
         if(istat .ne. 0) &
              call allocate_error("tam2r, sub tdreadf",0,istat)

         allocate (cint2r(jx_rstrt),STAT = istat)
         if(istat .ne. 0) &
              call allocate_error("cint2r, sub tdreadf",0,istat)
         call bcast(cint2r,zero,jx_rstrt)
         istatus= NF_INQ_VARID(ncid,'cint2',vid)
         istatus= NF_GET_VARA_DOUBLE(ncid,vid,(1),(jx_rstrt),cint2r)

         endif  !  On l_.eq.lrors

!.......................................................................
!        Read in and process f_rstrt one gen species at a time
!        f_rstrt dimensioned (1:iy,1:jx,1:setup0%lrz)
!BH180602:  Read does not exactly mirror the write in netcdfrw2.
!BH180602:  Check that it is OK!
!.......................................................................
         do k=1,ngen   !down to line 692

         ! First two indexes in 'f' array:
         start(1)=1
         start(2)=1
         count(1)=iy_rstrt   !=iy
         count(2)=jx_rstrt   !possibly diff from jx

         !YuP[2019-02-07] See netcdfrw2: 'f' can be 3D, or 4D, or 5D :
         !---1--- Case of saving f at each or several selected time steps
         !if (ngen.eq.1) then
         ! 'f' is 4D, NCDOUBLE,4, dimsf(1:4)={ydim,xdim,rdim, tdim}
         !                        start1(3)=ll (in do loop)
         !                        start1(4)=numrec1 (or numrecsave)
         !                        count1(1:4)={iy,jx,1,1} (in ll loop)
         !else  !ngen.ge.2
         ! 'f' is 5D, NCDOUBLE,5, dimsf(1:5)={ydim,xdim,rdim, gdim,tdim}
         !                        startg(3)=ll (in do loop)
         !                        startg(4)=k  (in do loop)
         !                        startg(5)=numrec1 (or numrecsave)
         !                        countg(1:5)={iy,jx,1,1,1} (in k,ll loops)
         !endif  !on ngen
         !
         !---2--- Case of saving f only at the end of run
         ! (and it is favg that was saved into setup0%mnemonic.nc file)
         !if (ngen.eq.1) then
         ! 'f' is 3D, NCDOUBLE,3, dimsf(1:3)={ydim,xdim,rdim}
         !                        start1(3)=ll (in do loop)
         !                        count1(1:3)={iy,jx,1} (in ll loop)
         !else  !ngen.ge.2
         ! 'f' is 4D, NCDOUBLE,4, dimsf(1:4)={ydim,xdim,rdim, gdim}
         !                        startg(3)=ll (in do loop)
         !                        startg(4)=k  (in do loop)
         !                        countg(1:4)={iy,jx,1,1} (in k,ll loops)
         !endif  !on ngen
         !Note: gdim is same as gen_species_dim (=ngen)
         ! Search "temp1(i,j)=f(" in netcdfrw2.

         ! These values can be different, dep. on cases above
         start(3)=l_ ! l_ is from each call of tdreadf(2)
         start(4)=k
         count(3)=1
         count(4)=1
         istatus= NF_INQ_VARID(ncid,'f',vid)
         !YuP: had errors with using f_rstrt; changed to wkpack
         istatus= NF_GET_VARA_DOUBLE(ncid,vid,start,count,wkpack)
                                     !This netCDF 'f' is a 4D object,
                                     !gen_species_dim, rdim, xdim, ydim.
                                     !Read it in to 3D f_rstrt
                                     !and process to code f(i,j,k,ll),
                                     !one gen species at at time.

         !YuP[2019-02-07] unpack wkpack(:) into temp_rstrt(i,j)
         call unpack21(temp_rstrt,0,iy_rstrt+1,0,jx_rstrt+1, &
                       wkpack,iy_rstrt,jx_rstrt)
         do j=1,jx_rstrt
            do i=1,iy_rstrt
               f_rstrt(i,j,l_)=temp_rstrt(i,j)
            enddo
         enddo

!        Code distribution function is f_cgs_units*vnorm**3.
!        The restart distribution is similarly normalized:
!        f_cgs_units_rstrt*vnorm_rstrt**3.
!        Thus, f_code=vnorm**3 * fcode_rstrt/vnorm_rstrt**3.
         renorm_f=(vnorm/vnorm_rstrt)**3

!        If vnorm_rstrt.lt.(vnorm+1d-10), will extrapolate beyond
!        edge of vnorm_rstrt grid, or beyond point where f(v)/f(v=0)
!        = 1.d-16.  Extrapolate linearly in ln(f_rstrt)
!        versus v_rstrt, using j_rstrt values at 75% and 95% of
!        jx_rstrt, separately for each pitch angle index i.
!        For values of code grid x(j) on the x_rstrt()-grid, linearly
!        interpolate in ln(f) vs v**2, to maintain a Maxwellian at
!        low velocity. Do one flux surface at a time.

!BH180606         endif  !On l_.eq.lrors)

!BH180606         if ( l_.le.lrors ) then !Calc f, one flux surface per call

!......................................................................
!        Cutoff ratio for extended distribution function, measured
!        from f(v=0).
!         foverf=1.e-16         ! Explore other values?
         foverf=1.e-24
!......................................................................

         ll=l_

            fj0=f_rstrt(1,1,ll)
            fof_fj0=foverf*fj0
            do i=1,iy
!              Find greatest jf_rstrt such that
!              f_rstrt(j)/f_rstrt(j=0).ge.foverf.
               jf_rstrt(i)=1
!               if (ll.eq.31) then
!                  write(*,*)'f_rstrt(1,j=1:10,31)=',
!     1              (f_rstrt(i,j,31),j=1,10)
!               endif
               do j=1,jx_rstrt-1
                  if (f_rstrt(i,j,ll).lt.fof_fj0) go to 10
                  jf_rstrt(i)=jf_rstrt(i)+1
               enddo
 10            continue
               if (jf_rstrt(i).lt.10) then
!MPIINSERT_IF_RANK_EQ_0
                  WRITE(*,*) 'tdreadf: i,jf_rstrt(i)=',i,jf_rstrt(i)
                  WRITE(*,*) 'tdreadf: Need jf_rstrt(i).ge.10'
!MPIINSERT_ENDIF_RANK
                  stop
               endif

               j1=0.75*jf_rstrt(i)
               j2=0.95*jf_rstrt(i)
               call bcast(f_rstrt_ln,zero,jx_rstrt)
               do j=1,j2
                  f_rstrt_ln(j)=log(f_rstrt(i,j,ll))
               enddo

!              Obtain extrapolation formula: ln(f)=aa*x+bb
               f1=f_rstrt_ln(j1)
               f2=f_rstrt_ln(j2)
!BH110525:  Could try v_rstrt2() in aa/bb, and mod extrapolation
!BH110525:  below.   There was an inconsistency.
               aa(i)=(f1-f2)/(v_rstrt(j1)-v_rstrt(j2))
               bb(i)=0.5*((f1+f2)-aa(i)*(v_rstrt(j1)+v_rstrt(j2)))
!
!              Use extrapolated values beyond 0.95 of jf_rstrt(i)
!              For velocities on the current code mesh less than
!              v(j2), will interpolate using cubic splines of
!              ln(f_rstrt(i,j=1:j2)) versus v**2.
!
               i1p(1)=4
               i1p(2)=4
               call coeff1(j2,v_rstrt2,f_rstrt_ln,d2fparr,i1p,1,workk)
               itab(1)=1  !tab(1) will contain interpolated value
               itab(2)=0
               itab(3)=0

               do j=1,jx
                  if (v2(j).lt.v_rstrt2(j2)) then
                     call terp1(j2,v_rstrt2,f_rstrt_ln,d2fparr,v2(j), &
                          1,tab,itab)  !interpolate
                     f(i,j,k,ll)=renorm_f*exp(tab(1))
                  else  !extrapolate
!BH110525                     f(i,j,1,ll)=renorm_f*exp(aa(i)*v2(j)+bb(i))
                     f(i,j,k,ll)=renorm_f*exp(aa(i)*v(j)+bb(i))
                  endif
                  !Set floor on extrapolated/interp'd f:
                  if (f(i,j,k,ll).lt.em100) f(i,j,k,ll)=em100
               enddo

            enddo  !On i

!.......................................................................
!     Calculate FSA density/energy in rstrt distribution and in code
!     distribution after interpolation and extrapolation.
!     Renormalize code distribution to rstrt FSA density.
!     Follow coding in diaggnde.f.
!     The y(i,l_) mesh is same for rstrt and code grids; thus
!     can use code cynt2, coss, tau, for density calc from f_rstrt.
!.......................................................................


            call tdnflxs(ll)
            call bcast(tam2r,zero,jx_rstrt)
            do i=1,iy
               do j=1,jx_rstrt
                  tam2r(j)=tam2r(j)+f_rstrt(i,j,l_)*cynt2(i,l_)* &
                       abs(coss(i,lmdpln_))*tau(i,lr_)
               enddo
            enddo
            hnr=zero
            snr=zero
            do j=1,jx_rstrt
               hnr=hnr+tam2r(j)*cint2r(j)
               snr=snr+tam2r(j)*cint2r(j)*x_rstrt(j)**2  !NonRel Energy
            enddo
            redenr=hnr/zmaxpsi(lr_)  !FSA density from f_rstrt
            senergyr=snr/hnr*fions(k)*vnorm_rstrt**2/vnorm2
                                       !FSA energy density from f_rstrt

            call bcast(tam2,zero,jx)
            do i=1,iy
               do j=1,jx
                  tam2(j)=tam2(j)+f(i,j,k,l_)*cynt2(i,l_)* &
                       abs(coss(i,lmdpln_))*tau(i,lr_)
               enddo
            enddo
            hn=zero
            sn=zero
            do j=1,jx
               hn=hn+tam2(j)*cint2(j)
               sn=sn+tam2(j)*cint2(j)*x(j)**2  !NonRel
            enddo
            reden_code=hn/zmaxpsi(lr_) !FSA density from code f, after
                                       !interpolation and extrapolation.
            reden_eps=(reden_code-redenr)/redenr
            senergy_code=sn/hn*fions(k)
                                           !FSA en density from code f.
            senergy_eps=(senergy_code-senergyr)/senergyr
!MPIINSERT_IF_RANK_EQ_0
            WRITE(*,*)'tdreadf:k=',k,', l_=',l_,' reden_code=', &
                 reden_code, &
                 ' Fractional density chng, reden_eps_=',reden_eps
            WRITE(*,*)'tdreadf: senergy_code=', &
                 senergy_code, &
                ' Fractional en density chng, senergy_eps_=',senergy_eps
!MPIINSERT_ENDIF_RANK
!           Renormalize code f to FSA density of restart f:
            reden_rat=redenr/reden_code
            do i=1,iy
               do j=1,jx
                  f(i,j,k,l_)=reden_rat*f(i,j,k,l_)
               enddo
            enddo

!BH180606         endif  !On l_.le.lrors


      enddo  !On k=1,ngen


      if (l_.eq.1) then

         istatus=NF_CLOSE(ncid)
         deallocate(f_rstrt,f_rstrt_ln,x_rstrt,v_rstrt,v_rstrt2, &
              y_rstrt,jf_rstrt,aa,bb,d2fparr,workk,v,v2, &
              tam2r,cint2r,STAT=ist1)
         if (ist1.ne.0) then
!     MPIINSERT_IF_RANK_EQ_0
            WRITE(*,*)'rdc_multi: dallocation error'
!     MPIINSERT_ENDIF_RANK
            STOP
         endif

      endif  !On l_.eq.1

      endif  !On nlrestrt

      return

!.......................................................................
!H080201 9210 format(1p10e13.6)
 9210 format(1p10e14.6)
!.......................................................................

      end subroutine tdreadf

!      end module tdreadf_mod
end module tdreadf_mod
