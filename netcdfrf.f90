module netcdfrf_mod

!
!
      module netcdfrf_mod
      use param_mod, only : nmodsa

contains

      integer, private ::  numrec(nmodsa) !-YuP-> added: as a function of krf

      save

      contains

      subroutine netcdfrf(kopt,krf)
      use param_mod
      use comm_mod
      implicit integer (i-n), real*8 (a-h,o-z)
      save


!     This subroutine uses netCDF-2 subroutines.
!     http://www.unidata.ucar.edu/packages/netcdf/index.html
!     (Had a problem with linking to netCDF-3 subroutines
!      on account of embedded underscores in their names).
!
!     This subroutine either reads an rf data file (see below),
!       or creates/writes a netCDF file with name "mnemonic"_rf.nc,
!       for RF data from cql3d.
!      (Additional standard netCDF output is given
!       in file "mnemonic".nc).
!
!     Action of the subroutine is controlled by kopt,krf,
!       time step counter n:
!       kopt=0, n=0, initialize "mnemonic"_rf.nc and write data
!       kopt=1, n.gt.0, write data (only time-dep data, unless n.eq.nstop)
!       kopt=2, read basic ray data from file rffile(krf)
!       kopt=3 is a limited part of kopt=2,
!            to only read the values of nray(krf),nharm(krf),freqcy(krf)
!       kopt=4 is a limited part of kopt=2
      ! kopt=2 part is called by urfread with krf=1:mrf to read all data.
      ! BUT if kopt=3 or 4, this part is called by urfsetup
      ! to only read the values of nray(krf),nharm(krf),freqcy(krf),
      ! and also nrayelt(:,krf) (when kopt=4).
!
!     krf= wave-type number, presently for reading input files
!            for up to 3 modes, with names specified by rffile(1:3);
!            by default, rffile()=rayech.nc, and/or rayrfw.nc,
!                                 and/or raylh.nc, depending on
!                                 namelist values of ech,lh,fw
!                                 See cqlinput for alternate namings.
!     BH060314:  Generalized to multi-wave-types/multi-harmonic, up to
!                maximum specified by parameter nmodsa.
!
!     YuP-100404: writing output ray data is now set up for all krf=1:mrf.
!     The data is saved into files "mnemonic"_krf###.nc, for each krf.

!     Data written consists of  ray data and ql diffusion
!     coefficient urfb (if netcdfshort.ne.'enabled' .and. .ne.'longer_b'
!     .and. .ne. 'long_urf'),
!     in which case urfb will be stored on the last step.
!     If netcdfshort.eq.'longer_b', the urfb stored at every time step.
!     if netcdfshort='long_urf', then all urf coefs (urfb,urfc,
!     urfe,urff) are stored at last time step.
!     Only data for one mode (or harmonic) is stored,
!     for the time being.
!     Also, some mesh data is included,  for possible
!     convenience in using this data file without
!     the general "mnemonic".nc file, e.g., with idl or ncbrowse.
!
!
!     NOTE:  Some of code comments are vestigial to previous code
!            fragments used in this code construction.
!            Take them with a grain of salt.
!
!
! --- include file for netCDF declarations
! --- (obtained from NetCDF distribution)
      include 'netcdf.inc'
!MPIINSERT_INCLUDE

! --- some stuff for netCDF file ---
      character*128 name
      character*3  krf_index               !-YuP-> to form file name
      integer ncid,vid,istatus
      integer xdim,ydim,rdim,r0dim,twodim,tdim
      integer nraydim,neltdim
      integer cdim

      integer dims(4),count(4),start(4)
      integer count1(4),start1(4)
      integer ray_dims(3),ray_count(3)
      integer r0_dims(2),r0_count(2)
      integer y_dims(2),y_count(2)
      integer tau_dims(2),tau_count(2)

      complex*16 cei

      data start/1,1,1,1/, start1/1,1,1,1/

!MPIINSERT_IF_RANK_NE_0_RETURN
! This subroutine is only called from MPI rank=0.

      cei=(0.,1.)



      WRITE(*,*)'Begin of netcdfrf, kopt,krf,n:',kopt,krf,n

      if(n.eq.0) call ibcast(numrec,1,SIZE(numrec)) ! counter for data recording



! --- begin if ---
      if ((kopt.eq.0 .or. kopt.eq.1) .and. &
          (irffile(krf).ne."combine2"))  then         !To line 1076

      !---- Form file name for each krf:
      WRITE(krf_index(1:3),'(3I1)') 0,0,0
      IF(                 krf.LE.9  ) WRITE(krf_index(3:3),'(I1)') krf
      IF(krf.GE.10  .AND. krf.LE.99 ) WRITE(krf_index(2:3),'(I2)') krf
      IF(krf.GE.100 .AND. krf.LE.999) WRITE(krf_index(1:3),'(I3)') krf
      if(krf.GE.1000) stop "Not setup for krf>999"
      t_= mnemonic(1:length_char(mnemonic))//'_krf'//krf_index//'.nc'
      !---- This file is created-then-closed at n=0,
      !---- opened-then-closed for writing data at n>0,
      !---- and opened-then-closed at n=nstop

!-----------------------------------------------------------------------
!     Only set up for cqlpmod.ne."enabled",ngen=1, for the time being.
!-----------------------------------------------------------------------
      if (cqlpmod.eq."enabled" ) then  !-YuP: removed  .or. ngen.ne.1   ???
         WRITE(*,*) 'WARNING: netcdfrf subroutine not implemented'
         return
      endif

!     Maximum iy as function of radius:
!$$$      iy=0
!$$$      do l=1,lrz
!$$$         iy=max(iy,iy_(l))
!$$$      enddo

!     Maximum number of ray elements per ray, for this specific krf (wave type):
      neltmax=0
      do iray=1,nray(krf)  !-YuP: nray(1) -> nray(krf)
         neltmax=max(neltmax,nrayelt(iray,krf)) ! nrayelt(*,1)->nrayelt(*,krf)
      enddo

!     Following counting vectors set up to facilitate writing
!     of the various arrays. Set up here to ensure initialization
!     in each call to the subroutine.

!$$$      count(1)=iy
      count(1)=iymax
      count(2)=jx
      count(3)=lrz
      count(4)=1

!$$$      count1(1)=iy
      count1(1)=iymax
      count1(2)=jx
      count1(3)=1
      count1(4)=1

      ray_count(1)=neltmax ! depends on krf
      ray_count(2)=nray(krf) !-YuP: nray(1) -> nray(krf)
      ray_count(3)=2

      r0_count(1)=lrzmax
      r0_count(2)=1

!$$$      y_count(1)=iy
      y_count(1)=iymax
      y_count(2)=lrz

      tau_count(1)=iy
      tau_count(2)=lrzmax


!.......................................................................
!l    1. Initialize part, creating new netcdf file
!

! --- begin if ---
      !WRITE(*,*)'netcdfrf: kopt,krf,n= ',kopt,krf,n
      if ((kopt.eq.0) .and. (n.eq.0)) then      !Down to line 634

!-----------------------------------------------------------------------
!
!l     1.1 create netCDF file and define dimensions,variables
!          and attributes
!

!.......................................................................
!l    1.1.1 create netCDF filename (Entering define mode.)
!     integer function nccre(filename,overwrite?,error_code)
!     Ref to page 46 of NetCDF-2 manual.
!     CLOBber old file, if it exists.
!     istatus is 0, if no errors.

!-YuP:      ncid=nccre(t_,NCCLOB,istatus)
      istatus = NF_CREATE(t_, NF_CLOBBER, ncid) !-YuP: NetCDF-f77
      call check_err(istatus)

!.......................................................................
!l    1.1.2 define dimensions
!     p. 67 of netcdf-2 manual
!     integer function ncddef(ncid,dim_name,dim_siz,error_code)
!       returns dimension id.

!-YuP:      xdim=ncddef(ncid,'x',jx,istatus)
!-YuP:      ydim=ncddef(ncid,'y',iymax,istatus)
!-YuP:      rdim=ncddef(ncid,'r',lrz,istatus)
!-YuP:      r0dim=ncddef(ncid,'r0',lrzmax,istatus)
!-YuP:      twodim=ncddef(ncid,'two',2,istatus)
!-YuP:      neltdim=ncddef(ncid,'neltmax',neltmax,istatus)
!-YuP:      nraydim=ncddef(ncid,'nrays',nray(1),istatus)
!-YuP:      cdim=ncddef(ncid,'cdim',8,istatus)

      istatus= NF_DEF_DIM(ncid, 'x',      jx ,      xdim)  !-YuP: NetCDF-f77
      istatus= NF_DEF_DIM(ncid, 'y',      iymax,    ydim)
      istatus= NF_DEF_DIM(ncid, 'r',      lrz,      rdim)
      istatus= NF_DEF_DIM(ncid, 'r0',     lrzmax,   r0dim)
      istatus= NF_DEF_DIM(ncid, 'two',    2,        twodim)
      istatus= NF_DEF_DIM(ncid, 'neltmax',neltmax,  neltdim) ! depends on krf
      istatus= NF_DEF_DIM(ncid, 'nrays',  nray(krf),nraydim) ! (1) -> (krf)
      istatus= NF_DEF_DIM(ncid, 'cdim',   8,        cdim)  !-YuP: NetCDF-f77

!     unlimited dimension for time, dimension name= 'time'
!     ncddef(ncid,'time',NF_UNLIMITED,error_code)
!-YuP:      tdim=ncddef(ncid,'time',NCUNLIM,istatus)
      istatus= NF_DEF_DIM(ncid, 'time',NF_UNLIMITED,tdim) !-YuP: NetCDF-f77

!     define vector of dimensions [for urfb(1:iy,1:jx,1:lrz,krf)], unlimited last
      dims(1)=ydim
      dims(2)=xdim
      dims(3)=rdim
      dims(4)=tdim

      ray_dims(1)=neltdim ! depends on krf
      ray_dims(2)=nraydim ! depends on krf
      ray_dims(3)=twodim

      r0_dims(1)=r0dim
      r0_dims(2)=tdim

      y_dims(1)=ydim
      y_dims(2)=rdim

      tau_dims(1)=ydim
      tau_dims(2)=r0dim

!.......................................................................
!l    1.1.3 define variables

!     Note, the variable IDs (denoted below as vid) are
!     not saved here in this subroutine; rather, the IDs
!     are retrieved from the netCDF data file, as needed,
!     by calling the netCDF routine ncvid.

!     netCDF variable_define:
!     integer function ncvdef2(ncid,variable_name,variable_type,
!                number_of_dimensions,
!                vector_for_length_per_dimension,error_code)
!       returns varid.
!     Note: Unlimited dimension must be last
!     Refer to p. 77, netcdf-2 manual

!     NCDOUBLE for REAL*8.  This is independent of the
!       internal representation of data on a particular
!       machine, e.g., 32- or 64-bit based calculations.


!--------------------------
!     Ray data
!--------------------------

      vid=ncvdef2(ncid,'nray',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,14, &
                  'Number of rays',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'freqcy',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,14, &
                  'Wave frequency',istatus)
      call check_err(istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,2, &
                 'Hz',istatus)

      vid=ncvdef2(ncid,'lh',NCCHAR,1,cdim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,37, &
                'Indicates raylh LH ray data file used',istatus)

      vid=ncvdef2(ncid,'fw',NCCHAR,1,cdim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,37, &
                'Indicates rayfw FW ray data file used',istatus)

      vid=ncvdef2(ncid,'ech',NCCHAR,1,cdim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,39, &
                'Indicates rayech ECH ray data file used',istatus)

      vid=ncvdef2(ncid,'nharm1',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,21, &
                  'First harmonic number',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'nharms',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name0',NCCHAR,39, &
                  'Number of harmonics, starting at nharm1',istatus)
      call ncaptc2(ncid,vid,'long_name1',NCCHAR,27, &
                  ' Only .gt.1 for one rf mode',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'nrayelt',NCLONG,1,nraydim,istatus) ! dep. on krf
      call ncaptc2(ncid,vid,'long_name',NCCHAR,35, &
                  'Number of ray elements for each ray',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'ws',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,29, &
                 'poloidal distance along a ray',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                 'cms',istatus)

      vid=ncvdef2(ncid,'spsi',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,20, &
                 'radial-like variable',istatus)

      vid=ncvdef2(ncid,'wr',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,12, &
                 'major radius',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                 'cms',istatus)

      vid=ncvdef2(ncid,'wphi',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,14, &
                 'toroidal angle',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                 'radians',istatus)

      vid=ncvdef2(ncid,'wz',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,15, &
                 'vertical height',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                 'cms',istatus)

      vid=ncvdef2(ncid,'wnpar',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,25, &
                 'parallel refractive index',istatus)

      vid=ncvdef2(ncid,'wnper',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30, &
                 'perpendicular refractive index',istatus)

      vid=ncvdef2(ncid,'delpwr',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,20, &
                 'power in ray channel',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,29, &
                 'erg/sec [BH180207:MKS, check]',istatus)

      vid=ncvdef2(ncid,'sdpwr',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,23, &
                 'power to ions (if incl)',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,8, &
                 'ergs/sec',istatus)

      vid=ncvdef2(ncid,'wdnpar',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,23, &
                 'effective width in npar',istatus)

!     Added 3rd dimension equal to 2 accomodates complex data.
      vid=ncvdef2(ncid,'cwexde',NCDOUBLE,3,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,25, &
                 'Complex Ex/E Polarization',istatus)

      vid=ncvdef2(ncid,'cweyde',NCDOUBLE,3,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,25, &
                 'Complex Ey/E Polarization',istatus)

      vid=ncvdef2(ncid,'cwezde',NCDOUBLE,3,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,25, &
                 'Complex Ez/E Polarization',istatus)

      vid=ncvdef2(ncid,'fluxn',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,23, &
                 'fluxn, Stix norm, |E|=1',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,13, &
                 'ergs/sec/cm^2',istatus)

      vid=ncvdef2(ncid,'sbtot',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,23, &
                 'Magnetic field strength',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5, &
                 'guass',istatus)

      vid=ncvdef2(ncid,'sene',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,17, &
                 'Density along ray',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,14, &
                 'particles/cm^3',istatus)

      vid=ncvdef2(ncid,'salphac',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30, &
                 'Collisional damping wavenumber',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5, &
                 '1/cms',istatus)

      vid=ncvdef2(ncid,'salphal',NCDOUBLE,2,ray_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,25, &
                 'Linear damping wavenumber',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,5, &
                 '1/cms',istatus)


!--------------------------
!     Additional Time-independent data
!--------------------------


      vid=ncvdef2(ncid,'lrzmax',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,25, &
                  'Number of radial surfaces',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'rya',NCDOUBLE,1,r0dim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,22, &
                 'Normalized radial mesh',istatus)

      vid=ncvdef2(ncid,'rhomax',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                 'cms',istatus)

      vid=ncvdef2(ncid,'lrz',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,29, &
                  'Number of FPd radial surfaces',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'lrindx',NCLONG,1,rdim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30, &
                 'Radial indices of FPd surfaces',istatus)

      vid=ncvdef2(ncid,'jx',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,27, &
                  'momentum-per-mass dimension',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'mrf',NCLONG,0,0,istatus) !-YuP: added
      call ncaptc2(ncid,vid,'long_name',NCCHAR,23, &
                  'number of rf wave types',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'mrfn',NCLONG,0,0,istatus) !-YuP[2017]added
      call ncaptc2(ncid,vid,'long_name',NCCHAR,59, &
       'number of rf modes (sum over all wave types and all nharms)', &
       istatus)
      call check_err(istatus)
      vid=ncvdef2(ncid,'krf',NCLONG,0,0,istatus) !-YuP: added
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30, &
                  'This file: rf wave type number',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'x',NCDOUBLE,1,xdim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,28, &
                 'normalized momentum-per-mass',istatus)

      vid=ncvdef2(ncid,'vnorm',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,33, &
                 'velocity (momentum-per-mass) norm',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7, &
                           'cms/sec',istatus)

      vid=ncvdef2(ncid,'enorm',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,20, &
                           'Energy normalization',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                           'keV',istatus)

      vid=ncvdef2(ncid,'iy',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,21, &
                  'Pitch angle dimension',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'y',NCDOUBLE,2,y_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,11, &
                 'pitch angle',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7, &
                 'radians',istatus)

      vid=ncvdef2(ncid,'iy_',NCLONG,1,rdim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,36, &
                  'Pitch angle dimension at each radius',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'itl',NCLONG,1,rdim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,26, &
                 'lower trapped-passing bndy',istatus)

      vid=ncvdef2(ncid,'itu',NCLONG,1,rdim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,26, &
                 'upper trapped-passing bndy',istatus)

      vid=ncvdef2(ncid,'tau',NCDOUBLE,2,tau_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,23, &
                 'tau_bounce * abs(speed)',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                           'cms',istatus)


!--------------------------
!     Time-dependent data
!--------------------------

      vid=ncvdef2(ncid,'time',NCDOUBLE,1,tdim,istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7, &
                           'seconds',istatus)


      if (netcdfshort.eq.'enabled') then
!         Do nothing: no storage defined.
      elseif (netcdfshort.eq.'longer_b') then
         vid=ncvdef2(ncid,'urfb',NCDOUBLE,4,dims,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,41, &
              'u^2 * BA QL Diff Coeff * v_par*tau_bounce',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,12, &
              'cgs/vnorm**4',istatus)
      elseif (netcdfshort.eq.'long_urf') then
         vid=ncvdef2(ncid,'urfb',NCDOUBLE,3,dims,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,42, &
              'u^2 * BA<QL Duu Coeff> * v_par0*tau_bounce',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,12, &
              'cgs/vnorm**4',istatus)
         vid=ncvdef2(ncid,'urfc',NCDOUBLE,3,dims,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,51, &
           'u * BA<c/(psi**0.5*c0) QL Dtu> * v_par0*tau_bounce',istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,39, &
              'Above: c=cos(theta),c0=c at B0, t=theta',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,12, &
              'cgs/vnorm**3',istatus)
!         vid=ncvdef2(ncid,'urfe',NCDOUBLE,3,dims,istatus)
!         call ncaptc2(ncid,vid,'long_name',NCCHAR,42,
!     +        'u * BA<c*s/(psi*c0) QL Dut> * v_par0*tau_bounce',istatus)
!         call ncaptc2(ncid,vid,'units',NCCHAR,12,
!     +        'cgs/vnorm**3',istatus)
!         vid=ncvdef2(ncid,'urff',NCDOUBLE,3,dims,istatus)
!         call ncaptc2(ncid,vid,'long_name',NCCHAR,55,
!     +        'BA<c**2*s/(psi**(3/2)*c0**2) QL Dtt> * v_par*tau_bounce',
!     +        istatus)
!         call ncaptc2(ncid,vid,'units',NCCHAR,12,
!     +        'cgs/vnorm**2',istatus)
!YuP[03/18/2015] urfe,urff are expressed through urfb,urfc
! No need to save them.
      else ! netcdfshort='disabled'; Standard o/p: f at last time step
         vid=ncvdef2(ncid,'urfb',NCDOUBLE,3,dims,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,42, &
              'u^2 * BA<QL Duu Coeff> * v_par0*tau_bounce',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,25, &
              'code units:  cgs/vnorm**4',istatus)
      endif


!.......................................................................
!l    1.1.4 end the define-mode and start the data-mode
!     p. 51-2 of netcdf manual

      istatus= NF_ENDDEF(ncid) !-YuP: NetCDF-f77
      call check_err(istatus)

!.......................................................................

!l    1.2 Write data
!
      start(4)=numrec(krf) ! counter for data recording: initially 1

! --- initialize data file ---
!     First get variable_id:
!        integer function ncvid(ncid,variable_name,error_code)
!        returns varid
!     Then write data with nc variable_put:
!        ncvpt(ncid,varid,start_indices,counts,values,error_code)
!
!     p. 79, 89 of netcdf-2 manual
!

!--------------------------
!     Time-independent data
!--------------------------

!-YuP:      vid=ncvid(ncid,'lrzmax',istatus)
      istatus= NF_INQ_VARID(ncid,'lrzmax',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_int2(ncid,vid,1,1,lrzmax,istatus)

!-YuP:      vid=ncvid(ncid,'rya',istatus)
      istatus= NF_INQ_VARID(ncid,'rya',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,1,lrzmax,rya(1),istatus)

!-YuP:      vid=ncvid(ncid,'rhomax',istatus)
      istatus= NF_INQ_VARID(ncid,'rhomax',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,1,1,rhomax,istatus)

!-YuP:      vid=ncvid(ncid,'lrz',istatus)
      istatus= NF_INQ_VARID(ncid,'lrz',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_int2(ncid,vid,1,1,lrz,istatus)

!-YuP:      vid=ncvid(ncid,'lrindx',istatus)
      istatus= NF_INQ_VARID(ncid,'lrindx',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_int2(ncid,vid,1,lrz,lrindx(1),istatus)

!-YuP:      vid=ncvid(ncid,'jx',istatus)
      istatus= NF_INQ_VARID(ncid,'jx',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_int2(ncid,vid,1,1,jx,istatus)

      istatus= NF_INQ_VARID(ncid,'mrf',vid)  !-YuP: added
      call ncvpt_int2(ncid,vid,1,1,mrf,istatus)

      istatus= NF_INQ_VARID(ncid,'mrfn',vid)  !-YuP[2017]added
      call ncvpt_int2(ncid,vid,1,1,mrfn,istatus)
      !number of rf modes (sum over all wave types and all nharms)
      istatus= NF_INQ_VARID(ncid,'krf',vid)  !-YuP: added
      call ncvpt_int2(ncid,vid,1,1,krf,istatus)

!-YuP:      vid=ncvid(ncid,'x',istatus)
      istatus= NF_INQ_VARID(ncid,'x',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,1,jx,x,istatus)

!-YuP:      vid=ncvid(ncid,'vnorm',istatus)
      istatus= NF_INQ_VARID(ncid,'vnorm',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,1,1,vnorm,istatus)

!-YuP:      vid=ncvid(ncid,'enorm',istatus)
      istatus= NF_INQ_VARID(ncid,'enorm',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,1,1,enorm,istatus)

!-YuP:      vid=ncvid(ncid,'iy',istatus)
      istatus= NF_INQ_VARID(ncid,'iy',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_int2(ncid,vid,1,1,iymax,istatus)

      call pack21(y,1,iy,1,lrors,tem1,iymax,lrors)
!-YuP:      vid=ncvid(ncid,'y',istatus)
      istatus= NF_INQ_VARID(ncid,'y',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start,y_count,tem1,istatus)

!-YuP:      vid=ncvid(ncid,'iy_',istatus)
      istatus= NF_INQ_VARID(ncid,'iy_',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_int2(ncid,vid,1,count(3),iy_,istatus)

!-YuP:      vid=ncvid(ncid,'itl',istatus)
      istatus= NF_INQ_VARID(ncid,'itl',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_int2(ncid,vid,1,count(3),itl_,istatus)

!-YuP:      vid=ncvid(ncid,'itu',istatus)
      istatus= NF_INQ_VARID(ncid,'itu',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_int2(ncid,vid,1,count(3),itu_,istatus)

      call pack21(tau,1,iy,1,lrzmax,tem1,iy,lrzmax)
!-YuP:      vid=ncvid(ncid,'tau',istatus)
      istatus= NF_INQ_VARID(ncid,'tau',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start,tau_count,tem1,istatus)

      endif  !kopt.eq.0 .and. n.eq.0

!-----------------------------------------------------------
!     Time-Dependent data:   (for  kopt=0,1 and all n.gte.0)
!-----------------------------------------------------------

!-YuP:      vid=ncvid(ncid,'time',istatus)
      istatus= NF_INQ_VARID(ncid,'time',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start(4),1,timet,istatus)

      if ( netcdfshort.eq.'longer_b' ) then
!-YuP:          vid=ncvid(ncid,'urfb',istatus)
         istatus= NF_INQ_VARID(ncid,'urfb',vid)  !-YuP: NetCDF-f77 get vid
         do ll=1,lrz
            start1(3)=ll
            call ncvpt_doubl2(ncid,vid,start1,count1, &
                 urfb(1,1,lrindx(ll),krf), istatus) !-YuP: (1)->(krf)
         enddo
      endif

!UNNECESSARY      call ncclos(ncid,istatus)
      istatus=NF_CLOSE(ncid)


!$$$C-----------------------------------------------------------------------
!$$$c
!$$$cl    2. Restart from previous run
!$$$c        The following data is NOT read in:
!$$$c          vnorm,x,y,z,dvol,flux,press,hflux,pdens,flux_a,flux_b,
!$$$c          hflux_a,hflux_b,pdens0
!$$$c
!$$$
!$$$c --- begin if ---
!$$$      if ((kopt.ne.0) .and. (n.eq.0)) then
!$$$
!$$$c.......................................................................
!$$$cl    2.1 Open previous netCDF file
!$$$
!$$$      ncid = ncopn(filename,NCWRITE,istatus)
!$$$
!$$$c.......................................................................
!$$$cl    2.2 read in dimension IDs and sizes
!$$$
!$$$      xdim = ncdid(ncid,'x',istatus)
!$$$      ydim = ncdid(ncid,'y',istatus)
!$$$      rdim = ncdid(ncid,'z',istatus)
!$$$      tdim = ncdid(ncid,'tau',istatus)
!$$$
!$$$c --- inquire about dimension sizes ---
!$$$c     ncdinq(netCDF_id, dimension_id_from_ncdid, returned_dim_name,
!$$$c     returned_dim_size)
!$$$c     Note: for unlimited dimension, returned_dim_size=current maximum
!$$$c     which is the same as the maximum record number
!$$$
!$$$      call ncdinq(ncid,ydim,name,iyp,istatus)
!$$$      call ncdinq(ncid,xdim,name,jxp,istatus)
!$$$      call ncdinq(ncid,rdim,name,kzp,istatus)
!$$$
!$$$c --- stop if dimension sizes don't match ---
!$$$      if ((iyp.ne.iy) .or. (jxp.ne.jx) .or. (kzp.ne.kz)) then
!$$$         write(6,*) "set the following parameter values in RUN_FPET"
!$$$         write(6,*) "  iy = ",iyp,"  jx = ",jxp, "  kz = ",kzp
!$$$         stop "non-matching dimensions in NCDFWRITE"
!$$$      endif
!$$$
!$$$c --- set the time-step counter ==> numrec(krf)
!$$$      call ncdinq(ncid,tdim,name,numrec(krf),istatus)
!$$$      start(4)=numrec(krf)
!$$$
!$$$c.......................................................................
!$$$cl    2.3 Read data
!$$$c
!$$$c     Here we read in only what's needed to re-start
!$$$c     the run from the previous time-step.
!$$$c
!$$$c     ncvgt(netCDF_id, variable_id_from_ncvid, vector_of_starting_index,
!$$$c     vector_of_lengths, returned_variable, integer_info)
!$$$c
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
!$$$c --- endif ((kopt.ne.0) .and. (n.eq.0)) ---
!$$$      endif

!-----------------------------------------------------------------------
!
!l    3. Periodic save at each time-steps
!
! --- begin if ---
      if (n.gt.0 .or. nstop.eq.0) then    !Down to line 759
      istatus = NF_OPEN(t_, NF_WRITE, ncid) !-YuP: NetCDF-f77
!.......................................................................
!l    3.1 increment the counter for data recording:
      numrec(krf)=numrec(krf)+1
      start(4)=numrec(krf)
      start1(4)=numrec(krf)
!.......................................................................
!l    3.2 Variables saved at each time-step

!-YuP:      vid=ncvid(ncid,'time',istatus)
      istatus= NF_INQ_VARID(ncid,'time',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start(4),1,timet,istatus)

      if ( netcdfshort.eq.'longer_b' ) then
!-YuP:          vid=ncvid(ncid,'urfb',istatus)
         istatus= NF_INQ_VARID(ncid,'urfb',vid)  !-YuP: NetCDF-f77 get vid
         do ll=1,lrz
            start1(3)=ll
            call ncvpt_doubl2(ncid,vid,start1,count1, &
                 urfb(1,1,lrindx(ll),krf), istatus) !-YuP: (1)->(krf)
         enddo
      endif
!UNNECESSARY      call ncclos(ncid,istatus)
      istatus=NF_CLOSE(ncid)
! --- endif (n.gt.0 .or. nstop.eq.0) ---
      endif

!.......................................................................
!l    3.3 Variables saved only at last time-step

! --- begin if ---
      if (n.eq.nstop) then        !Down to line 1071
      istatus = NF_OPEN(t_, NF_WRITE, ncid) !-YuP: NetCDF-f77

!.....................................................................
!     Adjust data to standard external format specified in urfread_.f.
!     (See urfread_.f. Perform inverse transformation)
!.....................................................................
!
        do iray=1,nray(krf)

        do  is=1,nrayelt(iray,krf)
           wdnpar(is,iray,krf)=abs(wdnpar(is,iray,krf))/wdscale(krf)
        enddo

!     Change sign of sbtot and wnpar, if bsign.lt.0.
!     (see urfwrite_.f and cqlinput_help)
        if (eqsource.eq."eqdsk" .and. bsign.lt.0.) then
           sbsign=sign(one,bsign)
           do is=1,nrayelt(iray,krf)
              sbtot(is,iray,krf)=sbtot(is,iray,krf)/bsign1(krf)
              wnpar(is,iray,krf)=sbsign*wnpar(is,iray,krf)
              cwexde(is,iray,krf)=sbsign*cwexde(is,iray,krf)
              cweyde(is,iray,krf)=sbsign*cweyde(is,iray,krf)
           enddo
        endif

!     fluxn is renormalized:
        do is=1,nrayelt(iray,krf)
           fluxn(is,iray,krf)=fluxn(is,iray,krf)*(8.*pi)/clight
        enddo  !  is

!     shift z-position of ray data, if eqdsk has been
!     vertically shifted in subroutine equilib:
        if (zshift.ne.0) then
           do is=1,nrayelt(iray,krf)
              wz(is,iray,krf)=wz(is,iray,krf)+zshift
           enddo
        endif

        enddo     !  iray

!..................................................................

!--------------------------
!     Put data
!--------------------------
!-YuP:      vid=ncvid(ncid,'nray',istatus)
      istatus= NF_INQ_VARID(ncid,'nray',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_int2(ncid,vid,1,1,nray(krf),istatus)

!-YuP:      vid=ncvid(ncid,'freqcy',istatus)
      istatus= NF_INQ_VARID(ncid,'freqcy',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,1,1,freqcy(krf),istatus)

!-YuP:      vid=ncvid(ncid,'lh',istatus)
      istatus= NF_INQ_VARID(ncid,'lh',vid)  !-YuP: NetCDF-f77 get vid
      call ncvptc2(ncid,vid,1,8,lh,8,istatus)

!-YuP:      vid=ncvid(ncid,'fw',istatus)
      istatus= NF_INQ_VARID(ncid,'fw',vid)  !-YuP: NetCDF-f77 get vid
      call ncvptc2(ncid,vid,1,8,fw,8,istatus)

!-YuP:      vid=ncvid(ncid,'ech',istatus)
      istatus= NF_INQ_VARID(ncid,'ech',vid)  !-YuP: NetCDF-f77 get vid
      call ncvptc2(ncid,vid,1,8,ech,8,istatus)

!-YuP:      vid=ncvid(ncid,'nharm1',istatus)
      istatus= NF_INQ_VARID(ncid,'nharm1',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_int2(ncid,vid,1,1,nharm1(krf),istatus)

!-YuP:      vid=ncvid(ncid,'nharms',istatus)
      istatus= NF_INQ_VARID(ncid,'nharms',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_int2(ncid,vid,1,1,nharms(krf),istatus)

!-YuP:      vid=ncvid(ncid,'nrayelt',istatus)
      istatus= NF_INQ_VARID(ncid,'nrayelt',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_int2(ncid,vid,1,nray(krf),nrayelt,istatus)

!BH170920: Added (1,1,krf) in below section of coding down to l 982, to get
!BH170920: krf dependent data going in to the respective krf data file.
      call pack21 &
        (ws(1,1,krf),1,nrayelts,1,nrayn,urftmp,neltmax,nray(krf)) ! depends on krf
!-YuP:      vid=ncvid(ncid,'ws',istatus)
      istatus= NF_INQ_VARID(ncid,'ws',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start,ray_count,urftmp,istatus)

      call pack21 &
        (spsi(1,1,krf),1,nrayelts,1,nrayn,urftmp,neltmax,nray(krf)) ! dep. on krf
!-YuP:      vid=ncvid(ncid,'spsi',istatus)
      istatus= NF_INQ_VARID(ncid,'spsi',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start,ray_count,urftmp,istatus)

      call pack21 &
        (wr(1,1,krf),1,nrayelts,1,nrayn,urftmp,neltmax,nray(krf)) ! depends on krf
!-YuP:      vid=ncvid(ncid,'wr',istatus)
      istatus= NF_INQ_VARID(ncid,'wr',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start,ray_count,urftmp,istatus)

      call pack21 &
        (wphi(1,1,krf),1,nrayelts,1,nrayn,urftmp,neltmax,nray(krf)) ! depends on krf
!-YuP:      vid=ncvid(ncid,'wphi',istatus)
      istatus= NF_INQ_VARID(ncid,'wphi',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start,ray_count,urftmp,istatus)

      call pack21 &
        (wz(1,1,krf),1,nrayelts,1,nrayn,urftmp,neltmax,nray(krf)) ! depends on krf
!-YuP:      vid=ncvid(ncid,'wz',istatus)
      istatus= NF_INQ_VARID(ncid,'wz',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start,ray_count,urftmp,istatus)

      call pack21 &
        (wnpar(1,1,krf),1,nrayelts,1,nrayn,urftmp,neltmax,nray(krf)) ! depends on krf
!-YuP:      vid=ncvid(ncid,'wnpar',istatus)
      istatus= NF_INQ_VARID(ncid,'wnpar',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start,ray_count,urftmp,istatus)

      call pack21 &
        (wnper(1,1,krf),1,nrayelts,1,nrayn,urftmp,neltmax,nray(krf)) ! depends on krf
!-YuP:      vid=ncvid(ncid,'wnper',istatus)
      istatus= NF_INQ_VARID(ncid,'wnper',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start,ray_count,urftmp,istatus)

      call pack21 &
        (delpwr(1,1,krf),1,nrayelts,1,nrayn,urftmp,neltmax,nray(krf)) ! depends on krf
!-YuP:      vid=ncvid(ncid,'delpwr',istatus)
      istatus= NF_INQ_VARID(ncid,'delpwr',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start,ray_count,urftmp,istatus)

      call pack21 &
        (sdpwr(1,1,krf),1,nrayelts,1,nrayn,urftmp,neltmax,nray(krf)) ! depends on krf
!-YuP:      vid=ncvid(ncid,'sdpwr',istatus)
      istatus= NF_INQ_VARID(ncid,'sdpwr',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start,ray_count,urftmp,istatus)

      call pack21 &
        (wdnpar(1,1,krf),1,nrayelts,1,nrayn,urftmp,neltmax,nray(krf)) ! depends on krf
!-YuP:      vid=ncvid(ncid,'wdnpar',istatus)
      istatus= NF_INQ_VARID(ncid,'wdnpar',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start,ray_count,urftmp,istatus)

      ii=0
      do j=1,nray(krf)
         do i=1,neltmax ! depends on krf
            ii=ii+1
            urftmp(ii)=0.5*(cwexde(i,j,krf)+conjg(cwexde(i,j,krf)))
         enddo
      enddo
      do j=1,nray(krf)
         do i=1,neltmax ! depends on krf
            ii=ii+1
            urftmp(ii)=-cei*0.5*(cwexde(i,j,krf)-conjg(cwexde(i,j,krf)))
         enddo
      enddo
!-YuP:      vid=ncvid(ncid,'cwexde',istatus)
      istatus= NF_INQ_VARID(ncid,'cwexde',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start,ray_count,urftmp,istatus)

      ii=0
      do j=1,nray(krf)
         do i=1,neltmax ! depends on krf
            ii=ii+1
            urftmp(ii)=0.5*(cweyde(i,j,krf)+conjg(cweyde(i,j,krf)))
         enddo
      enddo
      do j=1,nray(krf)
         do i=1,neltmax ! depends on krf
            ii=ii+1
            urftmp(ii)=-cei*0.5*(cweyde(i,j,krf)-conjg(cweyde(i,j,krf)))
         enddo
      enddo
!-YuP:      vid=ncvid(ncid,'cweyde',istatus)
      istatus= NF_INQ_VARID(ncid,'cweyde',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start,ray_count,urftmp,istatus)

      ii=0
      do j=1,nray(krf)
         do i=1,neltmax ! depends on krf
            ii=ii+1
            urftmp(ii)=0.5*(cwezde(i,j,krf)+conjg(cwezde(i,j,krf)))
         enddo
      enddo
      do j=1,nray(krf)
         do i=1,neltmax ! depends on krf
            ii=ii+1
            urftmp(ii)=-cei*0.5*(cwezde(i,j,krf)-conjg(cwezde(i,j,krf)))
         enddo
      enddo
!-YuP:      vid=ncvid(ncid,'cwezde',istatus)
      istatus= NF_INQ_VARID(ncid,'cwezde',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start,ray_count,urftmp,istatus)


      call pack21 &
        (fluxn(1,1,krf),1,nrayelts,1,nrayn,urftmp,neltmax,nray(krf)) ! depends on krf
!-YuP:      vid=ncvid(ncid,'fluxn',istatus)
      istatus= NF_INQ_VARID(ncid,'fluxn',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start,ray_count,urftmp,istatus)

      call pack21 &
        (sbtot(1,1,krf),1,nrayelts,1,nrayn,urftmp,neltmax,nray(krf)) ! depends on krf
!-YuP:      vid=ncvid(ncid,'sbtot',istatus)
      istatus= NF_INQ_VARID(ncid,'sbtot',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start,ray_count,urftmp,istatus)

      call pack21 &
        (sene(1,1,krf),1,nrayelts,1,nrayn,urftmp,neltmax,nray(krf)) ! depends on krf
!-YuP:      vid=ncvid(ncid,'sene',istatus)
      istatus= NF_INQ_VARID(ncid,'sene',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start,ray_count,urftmp,istatus)

      call pack21 &
        (salphac(1,1,krf),1,nrayelts,1,nrayn,urftmp,neltmax,nray(krf)) ! depends on krf
!-YuP:      vid=ncvid(ncid,'salphac',istatus)
      istatus= NF_INQ_VARID(ncid,'salphac',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start,ray_count,urftmp,istatus)

      call pack21 &
        (salphal(1,1,krf),1,nrayelts,1,nrayn,urftmp,neltmax,nray(krf)) ! depends on krf
!-YuP:      vid=ncvid(ncid,'salphal',istatus)
      istatus= NF_INQ_VARID(ncid,'salphal',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start,ray_count,urftmp,istatus)


! --- full update of data file (last time-step only) ---

      if (netcdfshort.ne.'enabled' .and. netcdfshort.ne.'longer_b') then
         istatus= NF_INQ_VARID(ncid,'urfb',vid)  !-YuP: NetCDF-f77 get vid
         do ll=1,lrz
            start1(3)=ll
            call ncvpt_doubl2(ncid,vid,start1,count1, &
                 urfb(1,1,lrindx(ll),krf),istatus) !-YuP: (1)->(krf)
         enddo

         if (netcdfshort.eq.'long_urf') then
            istatus= NF_INQ_VARID(ncid,'urfc',vid) !-YuP: NetCDF-f77 get vid
            do ll=1,lrz
               start1(3)=ll
               call ncvpt_doubl2(ncid,vid,start1,count1, &
                    urfc(1,1,lrindx(ll),krf),istatus) !-YuP: (1)->(krf)
            enddo
!            istatus= NF_INQ_VARID(ncid,'urfe',vid) !-YuP: NetCDF-f77 get vid
!            do ll=1,lrz
!               start1(3)=ll
!               call ncvpt_doubl2(ncid,vid,start1,count1,
!     +              urfe(1,1,lrindx(ll),krf),istatus) !-YuP: (1)->(krf)
!            enddo
!            istatus= NF_INQ_VARID(ncid,'urff',vid) !-YuP: NetCDF-f77 get vid
!            do ll=1,lrz
!               start1(3)=ll
!               call ncvpt_doubl2(ncid,vid,start1,count1,
!     +              urff(1,1,lrindx(ll),krf),istatus) !-YuP: (1)->(krf)
!            enddo
!YuP[03/18/2015] urfe,urff are expressed through urfb,urfc
! No need to save them.
         endif  !On netcdfshort
      endif  !On netcdfshort

!..................................................................
!     Re-Adjust data to internal format,
!     in case this is not last step.
!     (See urfread_.f)
!..................................................................
!
      do iray=1,nray(krf)

        do  is=1,nrayelt(iray,krf)
           wdnpar(is,iray,krf)=wdscale(krf)*abs(wdnpar(is,iray,krf))
        enddo

!     Change sign of sbtot and wnpar, if bsign.lt.0.
!     (see cqlinput_help)
        if (eqsource.eq."eqdsk" .and. bsign.lt.0.) then
           sbsign=sign(one,bsign)
           do is=1,nrayelt(iray,krf)
              sbtot(is,iray,krf)=bsign1(krf)*sbtot(is,iray,krf)
              if (sbtot(is,iray,krf).lt.0.) &
                 WRITE(*,*)'urfread: Sign Problem with sbtot:is,iray=', &
                 is,iray
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
        if (zshift.ne.0) then
           do is=1,nrayelt(iray,krf)
              wz(is,iray,krf)=wz(is,iray,krf)-zshift
           enddo
        endif

      enddo  !  iray


!..................................................................
      istatus = NF_CLOSE(ncid) !-YuP: NetCDF-f77
      call check_err(istatus)


! --- endif (n.eq.nstop) ---
      endif


!-----------------------------------------------------------------------
! --- endif ((kopt.eq.0) .or. (kopt.eq.1))
      endif




!-----------------------------------------------------------------------
!
!     READ NETCDF RAY DATA FILE.
!     (Parallels urfread_.f)
!
!-----------------------------------------------------------------------
! --- begin if ---
      if (kopt.eq.2  .or.  kopt.eq.3  .or.  kopt.eq.4) then !Down to 1443
      ! read data;
      ! This part is called by urfread with krf=1:mrf to read all data.
      ! BUT if kopt=3 or 4, this part is called by urfsetup
      ! to only read the values of nray(krf),nharm(krf),freqcy(krf),
      ! and also nrayelt(:,krf) (when kopt=4).

!     Open existing netCDF file

!-YuP      ncid=ncopn(rffile(krf),NCNOWRITE,istatus)
!-YuP      write(*,*)'after ncopn ncid=',ncid,'istatus',istatus
!-YuP      if (istatus.ne.0) then
!-YuP         write(*,*)'   ***   Problem opening rf .nc data file   ***'
!-YuP         Stop
!-YuP      endif
      istatus = NF_OPEN(rffile(krf), 0, ncid) !-YuP: NetCDF-f77
      if (istatus .NE. NF_NOERR) then         !-YuP: NetCDF-f77
         WRITE(*,*)'   ***   Problem opening rf .nc data file   ***'
         Stop
      endif                                   !-YuP: NetCDF-f77

!     Get dimension ID from dimension name
!-YuP:      neltdim=ncdid(ncid,'neltmax',istatus)
!-YuP:      nraydim=ncdid(ncid,'nrays',istatus)
!-YuP:      twodim=ncdid(ncid,'two',istatus)
      istatus= NF_INQ_DIMID(ncid,'neltmax',neltdim) ! dep. on krf
      istatus= NF_INQ_DIMID(ncid,'nrays',  nraydim) ! dep. on krf
      istatus= NF_INQ_DIMID(ncid,'two',    twodim)  !-YuP: NetCDF-f77 get twodim

!     Query netcdf file for dimensions:
!-YuP:      call ncdinq(ncid,neltdim,name,neltmax,istatus)
!-YuP:      call ncdinq(ncid,nraydim,name,nrays,istatus)
!-YuP:      call ncdinq(ncid,twodim,name,ntwo,istatus)
      istatus= NF_INQ_DIM(ncid,neltdim,name,neltmax)  ! dep. on krf
      istatus= NF_INQ_DIM(ncid,nraydim,name,nrays)    ! dep. on krf
      istatus= NF_INQ_DIM(ncid,twodim,name,ntwo)     !-YuP: NetCDF-f77 get ntwo

      !-YuP: Why needed ???      COMMENTING it out:
!BH110329  Now using urftmp temporary storage rather than temp1, etc.
!BH110329      if (neltmax*nrays*2 .gt. (iy+2)*(jx+2)) then
!BH110329        WRITE(*,*)'STOP in netcdfrf, neltmax*nrays*2.gt.(iy+2)*(jx+2)'
!BH110329        WRITE(*,*)'netcdfrf: neltmax,neltmax*nrays*2,(iy+2)*(jx+2)=',
!BH110329     +                       neltmax,neltmax*nrays*2,(iy+2)*(jx+2)
!BH110329c        STOP 'Increase iy or jx'
!BH110329      endif

!     Read data:
!     The following presupposes that the rank (number of dimensions)
!     is known for the input data. That is, the data conforms
!     to the above output data format.
!     (Checks could be placed in the code.)

!-YuP:      vid=ncvid(ncid,'nray',istatus)
      istatus= NF_INQ_VARID(ncid,'nray',vid)  !-YuP: NetCDF-f77 get vid
!-YuP:      call ncvgt(ncid,vid,1,1,nray(krf),istatus)
      istatus= NF_GET_VAR1_INT(ncid,vid,(1),nray(krf)) !-YuP: NetCDF-f77

!     Variable nharm in toray.nc and 3d.nc, is called nharm1 in
!     the cql3d mnemonic_rf.nc file.  Allow for this.
      !!!-YuP       call NCPOPT(NCVERBOSL)
!-YuP:      vid=ncvid(ncid,'nharm',istatus)
      istatus= NF_INQ_VARID(ncid,'nharm',vid)  !-YuP: NetCDF-f77 get vid
      !!!-YuP       call NCPOPT(NCVERBOS+NCFATAL)
      if (istatus.ne.0) then
!-YuP:          vid=ncvid(ncid,'nharm1',istatus)
      istatus= NF_INQ_VARID(ncid,'nharm1',vid)  !-YuP: NetCDF-f77 get vid
      endif
!-YuP:      call ncvgt(ncid,vid,1,1,nharm(krf),istatus)
      istatus= NF_GET_VAR1_INT(ncid,vid,(1),nharm(krf)) !-YuP: NetCDF-f77

!-YuP:      vid=ncvid(ncid,'freqcy',istatus)
      istatus= NF_INQ_VARID(ncid,'freqcy',vid)  !-YuP: NetCDF-f77 get vid
!-YuP:      call ncvgt(ncid,vid,1,1,freqcy(krf),istatus)
      istatus= NF_GET_VAR1_DOUBLE(ncid,vid,(1),freqcy(krf)) !-YuP: NetCDF-f77

      omega(krf)=2.*pi*freqcy(krf)


      if(kopt.eq.3) goto 300 !-> Finish reading (for urfsetup)
                             !   Close file. ----------------------

!..................................................................
!     Check code dimensions are sufficient for given ray data
!     Send message if nray.gt.nrayn
!..................................................................
!     Several files, corresponding to ray data from different wave
!     modes, may be read.
      if (nray(krf).gt.nrayn) then
        nr=nrayn
        WRITE(*,200) nray(krf),nr
        nray(krf)=nrayn
      endif
 200  format("Ray tracing code provides more rays than CQL parameter",/, &
        "nrayn is equipped to handle. Recompile CQL with nrayn=",i5,/, &
        "CQL will proceed with the calculation with nray reset to",i5)
!..................................................................

      ray_count(1)=neltmax ! depends on krf
      ray_count(2)=nrays   ! depends on krf
      ray_count(3)=2

!-YuP:      vid=ncvid(ncid,'nrayelt',istatus)
      istatus= NF_INQ_VARID(ncid,'nrayelt',vid)  !-YuP: NetCDF-f77 get vid
!-YuP:      call ncvgt(ncid,vid,1,ray_count(2),nrayelt(1,krf),istatus)
      istatus= NF_GET_VARA_INT(ncid,vid,(1),ray_count(2),nrayelt(:,krf))


      if(kopt.eq.4) goto 300 !-> Finish reading (for urfsetup)
                             !   Close file. ----------------------


      if (nrays.gt.nrayn .or. neltmax.gt.nrayelts) then
         WRITE(*,*)'netcdfrf:  nrays,neltmax= ', nrays,neltmax
         WRITE(*,*)'netcdfrf:  Need to be .le. nrayn,nrayelts'
         WRITE(*,*)'netcdfrf:  nrayn,nrayelts= ',nrayn,nrayelts
         STOP
      endif

!-YuP:      vid=ncvid(ncid,'ws',istatus)
      istatus= NF_INQ_VARID(ncid,'ws',vid)  !-YuP: NetCDF-f77 get vid
!-YuP:      call ncvgt(ncid,vid,start,ray_count,urftmp,istatus)
      istatus= NF_GET_VARA_DOUBLE(ncid,vid,start,ray_count,urftmp) !-YuP: NetCDF-f77
      call unpack21(ws(1,1,krf),1,nrayelts,1,nrayn, &
           urftmp,neltmax,nray(krf)) ! depends on krf

!-YuP:      vid=ncvid(ncid,'spsi',istatus)
      istatus= NF_INQ_VARID(ncid,'spsi',vid)  !-YuP: NetCDF-f77 get vid
!-YuP:      call ncvgt(ncid,vid,start,ray_count,urftmp,istatus)
      istatus= NF_GET_VARA_DOUBLE(ncid,vid,start,ray_count,urftmp) !-YuP: NetCDF-f77
      call unpack21(spsi(1,1,krf),1,nrayelts,1,nrayn, &
           urftmp,neltmax,nray(krf)) ! depends on krf

!-YuP:      vid=ncvid(ncid,'wr',istatus)
      istatus= NF_INQ_VARID(ncid,'wr',vid)  !-YuP: NetCDF-f77 get vid
!-YuP:      call ncvgt(ncid,vid,start,ray_count,urftmp,istatus)
      istatus= NF_GET_VARA_DOUBLE(ncid,vid,start,ray_count,urftmp) !-YuP: NetCDF-f77
      call unpack21(wr(1,1,krf),1,nrayelts,1,nrayn, &
           urftmp,neltmax,nray(krf)) ! depends on krf

!-YuP:      vid=ncvid(ncid,'wphi',istatus)
      istatus= NF_INQ_VARID(ncid,'wphi',vid)  !-YuP: NetCDF-f77 get vid
!-YuP:      call ncvgt(ncid,vid,start,ray_count,urftmp,istatus)
      istatus= NF_GET_VARA_DOUBLE(ncid,vid,start,ray_count,urftmp) !-YuP: NetCDF-f77
      call unpack21(wphi(1,1,krf),1,nrayelts,1,nrayn, &
           urftmp,neltmax,nray(krf)) ! depends on krf

!-YuP:      vid=ncvid(ncid,'wz',istatus)
      istatus= NF_INQ_VARID(ncid,'wz',vid)  !-YuP: NetCDF-f77 get vid
!-YuP:      call ncvgt(ncid,vid,start,ray_count,urftmp,istatus)
      istatus= NF_GET_VARA_DOUBLE(ncid,vid,start,ray_count,urftmp) !-YuP: NetCDF-f77
      call unpack21(wz(1,1,krf),1,nrayelts,1,nrayn, &
           urftmp,neltmax,nray(krf)) ! depends on krf

!-YuP:      vid=ncvid(ncid,'wnpar',istatus)
      istatus= NF_INQ_VARID(ncid,'wnpar',vid)  !-YuP: NetCDF-f77 get vid
!-YuP:      call ncvgt(ncid,vid,start,ray_count,urftmp,istatus)
      istatus= NF_GET_VARA_DOUBLE(ncid,vid,start,ray_count,urftmp) !-YuP: NetCDF-f77
      call unpack21(wnpar(1,1,krf),1,nrayelts,1,nrayn, &
           urftmp,neltmax,nray(krf)) ! depends on krf

!-YuP:      vid=ncvid(ncid,'wnper',istatus)
      istatus= NF_INQ_VARID(ncid,'wnper',vid)  !-YuP: NetCDF-f77 get vid
!-YuP:      call ncvgt(ncid,vid,start,ray_count,urftmp,istatus)
      istatus= NF_GET_VARA_DOUBLE(ncid,vid,start,ray_count,urftmp) !-YuP: NetCDF-f77
      call unpack21(wnper(1,1,krf),1,nrayelts,1,nrayn, &
           urftmp,neltmax,nray(krf)) ! depends on krf

!-YuP:      vid=ncvid(ncid,'delpwr',istatus)
      istatus= NF_INQ_VARID(ncid,'delpwr',vid)  !-YuP: NetCDF-f77 get vid
!-YuP:      call ncvgt(ncid,vid,start,ray_count,urftmp,istatus)
      istatus= NF_GET_VARA_DOUBLE(ncid,vid,start,ray_count,urftmp) !-YuP: NetCDF-f77
      call unpack21(delpwr(1,1,krf),1,nrayelts,1,nrayn, &
           urftmp,neltmax,nray(krf)) ! depends on krf

!     The toray .nc file puts out variable sdpwri rather
!     than sdpwr.  But toray always has sdpwri=0.
!     Trap around the difficulty, if there is no sdpwr variable.
      !!!-YuP       call NCPOPT(NCVERBOS)
!-YuP:      vid=ncvid(ncid,'sdpwr',istatus)
      istatus= NF_INQ_VARID(ncid,'sdpwr',vid)  !-YuP: NetCDF-f77 get vid
      if (istatus.ne.0) then
         WRITE(*,*)
         WRITE(*,*)'***************************************************'
         WRITE(*,*)'netcdfrf: sdpwr istatus =',istatus
         WRITE(*,*)'sdpwr variable not found, but not needed from toray'
         WRITE(*,*)'***************************************************'
         WRITE(*,*)
      else
!-YuP:         call ncvgt(ncid,vid,start,ray_count,urftmp,istatus)
      istatus= NF_GET_VARA_DOUBLE(ncid,vid,start,ray_count,urftmp) !-YuP: NetCDF-f77
         call unpack21(sdpwr(1,1,krf),1,nrayelts,1,nrayn, &
              urftmp,neltmax,nray(krf)) ! depends on krf
      endif
      !!!-YuP       call NCPOPT(NCVERBOS+NCFATAL)

!     The toray .nc file does not contain variable wdnpar.
!     If there is trouble reading wdnpar (istatus.ne.0),
!     then wdnpar is set equal to the standard value from
!     toray, corresponding to the DIII-D half angle at half-power EC
!     power pattern equal to 1.7 degrees.
      !!!-YuP       call NCPOPT(NCVERBOS)
!-YuP:      vid=ncvid(ncid,'wdnpar',istatus)
      istatus= NF_INQ_VARID(ncid,'wdnpar',vid)  !-YuP: NetCDF-f77 get vid
      if (istatus.ne.0) then
         WRITE(*,*)
         WRITE(*,*)'***************************************************'
         WRITE(*,*)'netcdfrf: wdnpar istatus =',istatus
         WRITE(*,*)'wdnpar variable not found, set equal to standard'
         WRITE(*,*)'DIII-D value 0.006297 for half-width at half-power'
         WRITE(*,*)'equal to 1.7 degrees.'
         WRITE(*,*)'***************************************************'
         WRITE(*,*)
         do j=1,nray(krf)
            do i=1,neltmax ! depends on krf
               wdnpar(i,j,krf)=0.006297
            enddo
         enddo
      else
!-YuP:         call ncvgt(ncid,vid,start,ray_count,urftmp,istatus)
         istatus= NF_GET_VARA_DOUBLE(ncid,vid,start,ray_count,urftmp) !-YuP: NetCDF-f77
         call unpack21(wdnpar(1,1,krf),1,nrayelts,1,nrayn, &
              urftmp,neltmax,nray(krf)) ! depends on krf
      endif
      !!!-YuP       call NCPOPT(NCVERBOS+NCFATAL)

!-YuP:      vid=ncvid(ncid,'cwexde',istatus)
      istatus= NF_INQ_VARID(ncid,'cwexde',vid)  !-YuP: NetCDF-f77 get vid
!-YuP:      call ncvgt(ncid,vid,start,ray_count,urftmp,istatus)
      istatus= NF_GET_VARA_DOUBLE(ncid,vid,start,ray_count,urftmp) !-YuP: NetCDF-f77
      ii=0
      len_ii=nray(krf)*neltmax ! depends on krf
      do j=1,nray(krf)
         do i=1,neltmax ! depends on krf
            ii=ii+1
            cwexde(i,j,krf)=cmplx(urftmp(ii),urftmp(ii+len_ii))
         enddo
      enddo

!-YuP:      vid=ncvid(ncid,'cweyde',istatus)
      istatus= NF_INQ_VARID(ncid,'cweyde',vid)  !-YuP: NetCDF-f77 get vid
!-YuP:      call ncvgt(ncid,vid,start,ray_count,urftmp,istatus)
      istatus= NF_GET_VARA_DOUBLE(ncid,vid,start,ray_count,urftmp) !-YuP: NetCDF-f77
      ii=0
      len_ii=nray(krf)*neltmax ! depends on krf
      do j=1,nray(krf)
         do i=1,neltmax ! depends on krf
            ii=ii+1
            cweyde(i,j,krf)=cmplx(urftmp(ii),urftmp(ii+len_ii))
         enddo
      enddo

!-YuP:      vid=ncvid(ncid,'cwezde',istatus)
      istatus= NF_INQ_VARID(ncid,'cwezde',vid)  !-YuP: NetCDF-f77 get vid
!-YuP:      call ncvgt(ncid,vid,start,ray_count,urftmp,istatus)
      istatus= NF_GET_VARA_DOUBLE(ncid,vid,start,ray_count,urftmp) !-YuP: NetCDF-f77
      ii=0
      len_ii=nray(krf)*neltmax ! depends on krf
      do j=1,nray(krf)
         do i=1,neltmax ! depends on krf
            ii=ii+1
            cwezde(i,j,krf)=cmplx(urftmp(ii),urftmp(ii+len_ii))
         enddo
      enddo

!-YuP:      vid=ncvid(ncid,'fluxn',istatus)
      istatus= NF_INQ_VARID(ncid,'fluxn',vid)  !-YuP: NetCDF-f77 get vid
!-YuP:      call ncvgt(ncid,vid,start,ray_count,urftmp,istatus)
      istatus= NF_GET_VARA_DOUBLE(ncid,vid,start,ray_count,urftmp) !-YuP: NetCDF-f77
      call unpack21(fluxn(1,1,krf),1,nrayelts,1,nrayn, &
           urftmp,neltmax,nray(krf)) ! depends on krf

!-YuP:      vid=ncvid(ncid,'sbtot',istatus)
      istatus= NF_INQ_VARID(ncid,'sbtot',vid)  !-YuP: NetCDF-f77 get vid
!-YuP:      call ncvgt(ncid,vid,start,ray_count,urftmp,istatus)
      istatus= NF_GET_VARA_DOUBLE(ncid,vid,start,ray_count,urftmp) !-YuP: NetCDF-f77
      call unpack21(sbtot(1,1,krf),1,nrayelts,1,nrayn, &
           urftmp,neltmax,nray(krf)) ! depends on krf

!-YuP:      vid=ncvid(ncid,'sene',istatus)
      istatus= NF_INQ_VARID(ncid,'sene',vid)  !-YuP: NetCDF-f77 get vid
!-YuP:      call ncvgt(ncid,vid,start,ray_count,urftmp,istatus)
      istatus= NF_GET_VARA_DOUBLE(ncid,vid,start,ray_count,urftmp) !-YuP: NetCDF-f77
      call unpack21(sene(1,1,krf),1,nrayelts,1,nrayn, &
           urftmp,neltmax,nray(krf)) ! depends on krf

!-YuP:      vid=ncvid(ncid,'salphac',istatus)
      istatus= NF_INQ_VARID(ncid,'salphac',vid)  !-YuP: NetCDF-f77 get vid
!-YuP:      call ncvgt(ncid,vid,start,ray_count,urftmp,istatus)
      istatus= NF_GET_VARA_DOUBLE(ncid,vid,start,ray_count,urftmp) !-YuP: NetCDF-f77
      call unpack21(salphac(1,1,krf),1,nrayelts,1,nrayn, &
           urftmp,neltmax,nray(krf)) ! depends on krf

!-YuP:      vid=ncvid(ncid,'salphal',istatus)
      istatus= NF_INQ_VARID(ncid,'salphal',vid)  !-YuP: NetCDF-f77 get vid
!BH050318      call unpack21(salphal(1,1,krf),1,nrayelts,1,nrayn,
!BH050318     1     urftmp,neltmax,nray(krf))
!-YuP:      call ncvgt(ncid,vid,start,ray_count,urftmp,istatus)
      istatus= NF_GET_VARA_DOUBLE(ncid,vid,start,ray_count,urftmp) !-YuP: NetCDF-f77
      call unpack21(salphal(1,1,krf),1,nrayelts,1,nrayn, &
           urftmp,neltmax,nray(krf)) ! depends on krf

!..................................................................
!     Adjust data to internal cql3d format.
!     (See urfread_.f)
!..................................................................
!
      do iray=1,nray(krf)

        do  is=1,nrayelt(iray,krf)
           wdnpar(is,iray,krf)=wdscale(krf)*abs(wdnpar(is,iray,krf))
        enddo

!     Change sign of sbtot and wnpar, if bsign.lt.0.
!     (see cqlinput_help and further explanation in urfread_.f)
        if (eqsource.eq."eqdsk" .and. bsign.lt.0.) then
           sbsign=sign(one,bsign)
           if (iray.eq.1) then
              bsign1(krf)=sign(one,sbtot(1,iray,krf))
              WRITE(*,*)'netcdf,kopt=2: eqsource,bsign,bsign1 = ', &
                   eqsource,bsign,bsign1
           endif
           do is=1,nrayelt(iray,krf)
              sbtot(is,iray,krf)=bsign1(krf)*sbtot(is,iray,krf)
              if (sbtot(is,iray,krf).lt.0.) &
                WRITE(*,*)'urfread: Sign Problem with sbtot, is,iray=', &
                is,iray
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
        if (zshift.ne.0) then
           do is=1,nrayelt(iray,krf)
              wz(is,iray,krf)=wz(is,iray,krf)-zshift
           enddo
        endif

      enddo   ! iray

300   istatus = NF_CLOSE(ncid) !-YuP: close file
      call check_err(istatus)

! --- endif --- (kopt = 2  or  3  or  4)
      endif


      return
      end subroutine





!==========================================================================



      subroutine netcdf_rdcb(krf)
      use param_mod
      use comm_mod
      implicit integer (i-n), real*8 (a-h,o-z)
      save

!MPIINSERT_INCLUDE

!     This subroutine uses netCDF-2 subroutines.
!     http://www.unidata.ucar.edu/packages/netcdf/index.html
!     (Had a problem with linking to netCDF-3 subroutines
!      on account of embedded underscores in their names).
!
!     This subroutine writes an DC file of urfb0 diffusion coeffs.
!       into a netCDF file with name "mnemonic"_rdc."krf".nc,
!       where "krf" is numeric value of krf.
!
! --- include file for netCDF declarations
! --- (obtained from NetCDF distribution)
      include 'netcdf.inc'

! --- some stuff for netCDF file ---
      character*128 name
      integer ncid,vid,istatus
      integer xdim,ydim,rdim,r0dim,tdim
      integer dims(4),start(4),count(4)
      integer start1(3),count1(3)
      integer y_dims(2),y_count(2)
      integer tau_dims(2),tau_count(2)

      data start/1,1,1,1/, start1/1,1,1/

!MPIINSERT_IF_RANK_NE_0_RETURN

      WRITE(*,*)'Begin of netcdf_rdc, krf=',krf

!$$$C-----------------------------------------------------------------------
!$$$C     Only set up for cqlpmod.ne."enabled",ngen=1, for the time being.
!$$$C-----------------------------------------------------------------------
!$$$
!$$$      if (ngen.ne.1) then
!$$$         WRITE(*,*) 'WARNING: netcdf_rdcb subroutine not implemented'
!$$$         WRITE(*,*) '         for ngen.gt.1'
!$$$         return
!$$$      endif

!BH180515:  Valid for multiple rdc files and ngen.ge.1

!     Following counting vectors set up to facilitate writing
!     of the various arrays. Set up here to ensure initialization
!     in each call to the subroutine.

      count(1)=iymax
      count(2)=jx
      count(3)=lrz
      count(4)=1

      count1(1)=iymax
      count1(2)=jx
      count1(3)=1

      y_count(1)=iymax
      y_count(2)=lrz

      tau_count(1)=iy
      tau_count(2)=lrzmax


!.......................................................................
!l    1.1.1 create netCDF filename (Entering define mode.)
!     integer function nccre(filename,overwrite?,error_code)
!     Ref to page 46 of NetCDF-2 manual.
!     CLOBber old file, if it exists.
!     istatus is 0, if no errors.

      write(t_,100) mnemonic(1:length_char(mnemonic)),krf
 100  format(a,"_rdc.",i1,".nc")
      WRITE(*,*)'netcdf_rdcb:t_ = ',t_(1:length_char(t_))
!-YuP:      ncid=nccre(t_,NCCLOB,istatus)
      istatus = NF_CREATE(t_, NF_CLOBBER, ncid) !-YuP: NetCDF-f77
      call check_err(istatus)

!.......................................................................
!l    1.1.2 define dimensions
!     p. 67 of netcdf-2 manual
!     integer function ncddef(ncid,dim_name,dim_siz,error_code)
!       returns dimension id.

!-YuP:      xdim=ncddef(ncid,'x',jx,istatus)
!-YuP:      ydim=ncddef(ncid,'y',iymax,istatus)
!-YuP:      rdim=ncddef(ncid,'r',lrz,istatus)
!-YuP:      r0dim=ncddef(ncid,'r0',lrzmax,istatus)
      istatus= NF_DEF_DIM(ncid, 'x',  jx ,    xdim)  !-YuP: NetCDF-f77
      istatus= NF_DEF_DIM(ncid, 'y',  iymax,  ydim)
      istatus= NF_DEF_DIM(ncid, 'r',  lrz,    rdim)
      istatus= NF_DEF_DIM(ncid, 'r0', lrzmax, r0dim) !-YuP: NetCDF-f77


!     unlimited dimension for time, dimension name= 'time'
!     ncddef(ncid,'time',NF_UNLIMITED,error_code)
!-YuP:      tdim=ncddef(ncid,'time',NCUNLIM,istatus)
      istatus= NF_DEF_DIM(ncid, 'time',NF_UNLIMITED,tdim) !-YuP: NetCDF-f77

!     define vector of dimensions [for urfb(iy,jx,lrz,mrfn)], unlimited last
      dims(1)=ydim
      dims(2)=xdim
      dims(3)=rdim
      dims(4)=tdim

      y_dims(1)=ydim
      y_dims(2)=rdim

      tau_dims(1)=ydim
      tau_dims(2)=r0dim

!--------------------------
!     Additional Time-independent data
!--------------------------


      vid=ncvdef2(ncid,'lrzmax',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,25, &
                  'Number of radial surfaces',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'rya',NCDOUBLE,1,r0dim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,22, &
                 'Normalized radial mesh',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'rpconz',NCDOUBLE,1,r0dim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,39, &
                 'Major radius at outside of flux surface',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'rhomax',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                 'cms',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'lrz',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,29, &
                  'Number of FPd radial surfaces',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'lrindx',NCLONG,1,rdim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30, &
                 'Radial indices of FPd surfaces',istatus)

      vid=ncvdef2(ncid,'jx',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,27, &
                  'momentum-per-mass dimension',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'x',NCDOUBLE,1,xdim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,28, &
                 'normalized momentum-per-mass',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'vnorm',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,33, &
                 'velocity (momentum-per-mass) norm',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7, &
                           'cms/sec',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'enorm',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,20, &
                           'Energy normalization',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                           'keV',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'iy',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,21, &
                  'Pitch angle dimension',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'y',NCDOUBLE,2,y_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,11, &
                 'pitch angle',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7, &
                 'radians',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'iy_',NCLONG,1,rdim,istatus)

      vid=ncvdef2(ncid,'itl',NCLONG,1,rdim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,26, &
                 'lower trapped-passing bndy',istatus)

      vid=ncvdef2(ncid,'itu',NCLONG,1,rdim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,26, &
                 'upper trapped-passing bndy',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'tau',NCDOUBLE,2,tau_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,23, &
                 'tau_bounce * abs(speed)',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                           'cms',istatus)


      vid=ncvdef2(ncid,'rdcb',NCDOUBLE,3,dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,40, &
                 'u^2 * BA QL Diff Coeff * v_par*tau_bounce',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,12, &
                 'cgs/vnorm**4',istatus)
      call check_err(istatus)

      if (netcdfshort.eq.'long_urf') then

         vid=ncvdef2(ncid,'rdcc',NCDOUBLE,3,dims,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,51, &
           'u * BA<c/(psi**0.5*c0) QL Dtu> * v_par0*tau_bounce',istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,39, &
              'Above: c=cos(theta),c0=c at B0, t=theta',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,12, &
              'cgs/vnorm**3',istatus)
         vid=ncvdef2(ncid,'rdce',NCDOUBLE,3,dims,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,42, &
              'u * BA<c*s/(psi*c0) QL Dut> * v_par0*tau_bounce',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,12, &
              'cgs/vnorm**3',istatus)
         vid=ncvdef2(ncid,'rdcf',NCDOUBLE,3,dims,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,55, &
              'BA<c**2*s/(psi**(3/2)*c0**2) QL Dtt> * v_par*tau_bounce', &
              istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,12, &
              'cgs/vnorm**2',istatus)
         call check_err(istatus)

      endif  ! On netcdfshort

!.......................................................................
!l    1.1.4 end the define-mode and start the data-mode
!     p. 51-2 of manual

!-YuP:      call ncendf(ncid,istatus)
      istatus= NF_ENDDEF(ncid) !-YuP: NetCDF-f77
      call check_err(istatus)


!-----------------------------------------------------------------------
!l    1.2 Write data
!

! --- set the time-step counter ==> numrec_rdcb
      numrec_rdcb=1
      start(4)=numrec_rdcb

! --- initialize data file ---
!     First get variable_id:
!        integer function ncvid(ncid,variable_name,error_code)
!        returns varid
!     Then write data with nc variable_put:
!        ncvpt(ncid,varid,start_indices,counts,values,error_code)
!
!     p. 79, 89 of netcdf-2 manual


!--------------------------
!     Time-independent data
!--------------------------

!-YuP:      vid=ncvid(ncid,'lrzmax',istatus)
      istatus= NF_INQ_VARID(ncid,'lrzmax',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_int2(ncid,vid,1,1,lrzmax,istatus)

!-YuP:      vid=ncvid(ncid,'rya',istatus)
      istatus= NF_INQ_VARID(ncid,'rya',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,1,tau_count(2),rya(1),istatus)

!-YuP:      vid=ncvid(ncid,'rpconz',istatus)
      istatus= NF_INQ_VARID(ncid,'rpconz',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,1,tau_count(2),rpconz(1),istatus)

!-YuP:      vid=ncvid(ncid,'rhomax',istatus)
      istatus= NF_INQ_VARID(ncid,'rhomax',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,1,1,rhomax,istatus)

!-YuP:      vid=ncvid(ncid,'lrz',istatus)
      istatus= NF_INQ_VARID(ncid,'lrz',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_int2(ncid,vid,1,1,lrz,istatus)

!-YuP:      vid=ncvid(ncid,'lrindx',istatus)
      istatus= NF_INQ_VARID(ncid,'lrindx',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_int2(ncid,vid,1,lrz,lrindx(1),istatus)

!-YuP:      vid=ncvid(ncid,'jx',istatus)
      istatus= NF_INQ_VARID(ncid,'jx',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_int2(ncid,vid,1,1,jx,istatus)

!-YuP:      vid=ncvid(ncid,'x',istatus)
      istatus= NF_INQ_VARID(ncid,'x',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,1,jx,x,istatus)

!-YuP:      vid=ncvid(ncid,'vnorm',istatus)
      istatus= NF_INQ_VARID(ncid,'vnorm',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,1,1,vnorm,istatus)

!-YuP:      vid=ncvid(ncid,'enorm',istatus)
      istatus= NF_INQ_VARID(ncid,'enorm',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,1,1,enorm,istatus)

!-YuP:      vid=ncvid(ncid,'iy',istatus)
      istatus= NF_INQ_VARID(ncid,'iy',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_int2(ncid,vid,1,1,iymax,istatus)

      if (iy*lrors.gt.iyjx2) stop 'netcdfrf:  Need to set jx>lrza'
      call pack21(y,1,iy,1,lrors,tem1,iymax,lrors)
!-YuP:      vid=ncvid(ncid,'y',istatus)
      istatus= NF_INQ_VARID(ncid,'y',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start,y_count,tem1,istatus)

!-YuP:      vid=ncvid(ncid,'iy_',istatus)
      istatus= NF_INQ_VARID(ncid,'iy_',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_int2(ncid,vid,1,lrz,iy_,istatus)

!-YuP:      vid=ncvid(ncid,'itl',istatus)
      istatus= NF_INQ_VARID(ncid,'itl',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_int2(ncid,vid,1,lrz,itl_,istatus)

!-YuP:      vid=ncvid(ncid,'itu',istatus)
      istatus= NF_INQ_VARID(ncid,'itu',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_int2(ncid,vid,1,lrz,itu_,istatus)

      call pack21(tau,1,iy,1,lrzmax,tem1,iy,lrzmax)
!-YuP:      vid=ncvid(ncid,'tau',istatus)
      istatus= NF_INQ_VARID(ncid,'tau',vid)  !-YuP: NetCDF-f77 get vid
      call ncvpt_doubl2(ncid,vid,start,tau_count,tem1,istatus)

!     rdcb coefficient write

!-YuP:      vid=ncvid(ncid,'rdcb',istatus)
      istatus= NF_INQ_VARID(ncid,'rdcb',vid)  !-YuP: NetCDF-f77 get vid
      do ll=1,lrz
         start1(3)=ll
         call ncvpt_doubl2(ncid,vid,start1,count1,rdcb(1,1,lrindx(ll), &
                    krf),istatus)
         if (istatus.ne.NF_NOERR) WRITE(*,*) 'netcdf_rdc, ll= :',ll
         call check_err(istatus)
      enddo

!     Additional rfc.... coefficients write
      if (netcdfshort.eq.'long_urf') then

      istatus= NF_INQ_VARID(ncid,'rdcc',vid)  !-YuP: NetCDF-f77 get vid
      do ll=1,lrz
         start1(3)=ll
         call ncvpt_doubl2(ncid,vid,start1,count1,rdcc(1,1,lrindx(ll), &
                    krf),istatus)
         if (istatus.ne.NF_NOERR) WRITE(*,*) 'netcdf_rdc, ll= :',ll
         call check_err(istatus)
      enddo

      istatus= NF_INQ_VARID(ncid,'rdce',vid)  !-YuP: NetCDF-f77 get vid
      do ll=1,lrz
         start1(3)=ll
         call ncvpt_doubl2(ncid,vid,start1,count1,rdce(1,1,lrindx(ll), &
                    krf),istatus)
         if (istatus.ne.NF_NOERR) WRITE(*,*) 'netcdf_rdc, ll= :',ll
         call check_err(istatus)
      enddo

      istatus= NF_INQ_VARID(ncid,'rdcf',vid)  !-YuP: NetCDF-f77 get vid
      do ll=1,lrz
         start1(3)=ll
         call ncvpt_doubl2(ncid,vid,start1,count1,rdcf(1,1,lrindx(ll), &
                    krf),istatus)
         if (istatus.ne.NF_NOERR) WRITE(*,*) 'netcdf_rdc, ll= :',ll
         call check_err(istatus)
      enddo


      endif  ! On netcdfshort

!     Close netcdf file
!-YuP:      call ncclos(ncid,istatus)
      istatus = NF_CLOSE(ncid) !-YuP: NetCDF-f77

      call check_err(istatus)

      return
      end subroutine

      end module netcdfrf_mod
end module netcdfrf_mod
