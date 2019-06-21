module tdtrdfus_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double
  use netcdfrw2_mod, only : ncvdef0,ncvdef2,ncaptc2,ncvptc0, &
                            length_char,check_err

  use bcast_mod, only : bcast
  use cqlcomm_mod
  !XXXXX use pack21_mod, only : pack21
  !XXXXX use pack21_mod, only : unpack21
  use r8subs_mod, only : luf
  use tdnflxs_mod, only : tdnflxs

  !---END USE


  external pack21  !XXXXX
  external unpack21

!
!

contains

      subroutine tdtrdfus
        use param_mod
        use cqlconf_mod, only : setup0
      use cqlcomm_mod, only : d_rr
      use cqlcomm_mod, only : transp, difus_io, difus_type
      use cqlcomm_mod, only : enerkev
      use cqlcomm_mod, only : ryain
      use cqlcomm_mod, only : difin
      use cqlcomm_mod, only : y

      implicit integer (i-n), real(c_double) (a-h,o-z)

!..............................................................
!     This routine computes the radial diffusion coefficients as
!     a function of mu or y and of rho (radius) - at a
!     given v and k (species index).
!     Present possibilities:  see cqlinput_help
!     There are three general possibilities:
!        1) Simplified ion neoclassical model taking pitch
!           angle variation uniform (should actually be strongly peaked
!           at tr-pass bndry) and proportional to 1/vel above vth_ion.
!        2) User specified Drr giving phenomenological radius and
!           velocity dependences as specified in cqlinput_help.
!        3) Read in a .nc file (named drr_in) which contains the diffusion
!           coefficient (d_rr) on a compatible v,theta,rho grid.
!    Alternatively, a .nc file is created (drr_out_<mnenomic>) containing
!    a d_rr computed in this subroutine and the grids.
!    (Also, d_r advection is subsequently added.)
!..............................................................

!%OS
      dimension drshape(0:setup0%lrz), rhotr(0:setup0%lrz)

      if (transp.eq."disabled") return

      d_rr=zero !YuP[2019-06-08]was call bcast(d_rr,zero,iyjx2*ngen*(lrz+1))
      drshape=zero !YuP[2019-06-08]was call bcast(drshape(0),zero,lrz+1)

!      write(*,*)'tdtrdfus: difus_type(),',(difus_type(k),k=1,ngen)

      if (difus_io(1).ne."drrin") then

      do k=1,ngen

      if ( difus_type(k) .eq. "neo_smpl" &
              .or. difus_type(k) .eq. "neo_plus" &
              .or. difus_type(k) .eq. "neo_trap" &
              .or. difus_type(k) .eq. "neo_trpp" ) then

      do l=0,setup0%lrz-1
         ilr=setup0%lrindx(l)+1        !Need to check why index is 0,setup0%lrz-1 below.
                                !BH051101: bndry condition is zero flux at
                                !rho=0, and Maxwl distn at rho=a.
                                !Zero flux at rho=0 is implemented with
                                !Drr=0 at l=0.  See McCoy notes, 10/12/90,
                                !pp. 8-10.
         vth_cnst=3.*sqrt(2.)*vth(kionn,ilr) !Matching v-dependence in tdrmshst
         rooteps=eps(ilr)**0.5  !This could cause too large diff as eps==>0.
         call tdnflxs(ilr)
         do j=1,jx

            if (x(j)*vnorm .lt. vth_cnst) then
               do i=1,iymax
                  d_rr(i,j,k,ilr)=drr_gs(ilr)
               enddo
            else
               do i=1,iymax
                  d_rr(i,j,k,ilr)=drr_gs(ilr)*vth_cnst/(x(j)*vnorm)
               enddo
            endif

            if (difus_type(k) .eq. "neo_trap") then
               do i=1,itl-1
                  d_rr(i,j,k,ilr)=0.
               enddo
               do i=itl,itu
                  d_rr(i,j,k,ilr)=d_rr(i,j,k,ilr)/rooteps
               enddo
               do i=itu+1,iy
                  d_rr(i,j,k,ilr)=0.
               enddo
            endif

         enddo
!     c        Checking calc of d_rr vel dependence against seperate
!     c        tdtrdfus calc of d_rr_b  ==> OK
         eighty=80.             ! getting consistent precision
         jkev=luf(eighty,enerkev,jx)
         jkev=jkev-1
      if ( l.eq.0) write(*,*)'tdtrdfus: enerkev', &
         (enerkev(j,k),j=1,jx) !YuP[2018-01-08] added 2nd index (k)
      if ( l.eq.0) write(*,*)'tdtrdfus: d_rr',(d_rr(iyh,j,k,ilr),j=1,jx)
      write(*,*)'tdtrdfus: jkev,d_rr(iyh,1,k,ilr),d_rr(iyh,jkev,k,ilr) ' &
              ,jkev,d_rr(iyh,1,k,ilr),d_rr(iyh,jkev,k,ilr)
      enddo  !  on l

      endif  ! on difus_type

      if (difus_type(k) .eq. "specify" &
              .or. difus_type(k) .eq. "neo_plus" &
              .or. difus_type(k) .eq. "neo_trpp") then

!     Test for radial diffusion coeff cnst in velocity.
      icnst=0
      do ii=1,4
         if (difus_vshape(ii).ne.zero) icnst=1
      enddo

!   check if difin non-zero to test if should be used
      rsum = 0.
      do ii = 1,njene
        rsum = rsum+abs(difin(ii))
      enddo
!   mesh to compute drshape on from difin(ryain)
      do ii=0,setup0%lrz
        rhotr(ii) = rrz(ii)/radmin
      enddo

      rdefr = -1.
      if (rsum.gt.0.) then
        rdefr = 1.
        call ryaintorz(njene,ryain,difin,setup0%lrz,rhotr(0),drshape(0))
      endif

      do l=0,setup0%lrz-1
         ilr=setup0%lrindx(l)
         if (l.eq.0) then   !setup0%lrindx(0)=0
            difusr1=0.
            rshape=0.
            rshape1=0.
         else
            difusr1=difusr
            if (rdefr.lt.0.5) then
              rshape=tdtrrshape(ilr)
            else
              rshape=drshape(ilr)
            endif
            if (difus_rshape(8).ne.zero.and. &
                         difus_type(k) .eq. "neo_plus") then
               rshape1=tdtrrshape1(ilr)
            else
               rshape1=0.
            endif
         endif
!         write(*,*)'tdtrdfus, rshape,rshape1:',rshape,rshape1
         if (l.ne.0) call tdtrvshape(k,l) !Returns shape in temp1
         do j=1,jx
            if (l.eq.0 .or. icnst.eq.0) then
               do i=1,iymax
!BH050921 Substantial bug fix, adding  *rshape     d_rr(i,j,k,l)=difusr1
                  d_rr(i,j,k,l)=d_rr(i,j,k,l) &
                              + difusr1*rshape &
                              + difusr1*rshape1
               enddo
            else
               do i=1,iytr(l)
                  d_rr(idx(i,l),j,k,l)=d_rr(idx(i,l),j,k,l) &
                       + difusr1*rshape*temp1(idx(i,l),j) &
                       + difusr1*rshape1
               enddo
            endif
         enddo  !On j

!        write(*,*)'tdtrdfus:l,(temp1(1,j),j=1,jx)',l,(temp1(1,j),j=1,jx)
!        write(*,*)'tdtrdfus: l, d_rr',l,(d_rr(1,j,k,l),j=1,jx)

       enddo  ! on l

      endif   ! on difus_type
      enddo   ! on k

!$$$  This will be done at end of code, after final calc of d_r.
!$$$      elseif (difus_io(1).eq."drrout") then
!$$$
!$$$         call diffus_io(1)  !Writes d_rr for k=1 to .nc file
!$$$
!$$$      elseif (difus_io(1).eq."drrdrout") then
!$$$
!$$$         call diffus_io(2)  !Writes d_rr for k=1 to .nc file
!$$$
!$$$  This will be done a beginning of the code, as alternative to
!$$$  call tdtrdfus.
!$$$      elseif (difus_io(1).eq."drrin") then
!$$$
!$$$         call diffus_io(3)  !Read d_rr for k=1 to .nc file
!$$$
!$$$      elseif (difus_io(1).eq."drrdrin") then
!$$$
!$$$         call diffus_io(4) !Reads d_rr and d_r for k=1 from .nc file

      endif  !On difus_io(1).ne."drrin"

      return
      end subroutine tdtrdfus
!
!
      real(c_double) function tdtrrshape(lr)
      use param_mod
      use cqlcomm_mod, only : zeff, temp, reden, rrz
      use cqlcomm_mod, only : difus_rshape
      implicit integer (i-n), real(c_double) (a-h,o-z)

!.......................................................................
!     This routine computes additional radial shape function for
!     the radial diffusion coefficient.
!     Coefficients in the following expression are input
!     throught the namelist variable difus_rshape(1:7).
!.......................................................................

!
      if (lr.eq.0) then

         tdtrrshape=0.

      elseif (lr.le.setup0%lrzmax-1) then

         if (kelecm.ne.0) then  !kelecm=0 when colmodl=1
            kk=kelecm
         else
            kk=kelecg
         endif

         tdtrrshape=(difus_rshape(1)+difus_rshape(2)* &
                    (rrz(lr)/radmin)**difus_rshape(3))**difus_rshape(4)
         tdtrrshape=tdtrrshape*(0.5*(reden(kk,lr)+reden(kk,lr+1))/ &
                    reden(kk,0))**difus_rshape(5)
         tdtrrshape=tdtrrshape*(0.5*(temp(kk,lr)+temp(kk,lr+1))/ &
                    temp(kk,0))**difus_rshape(6)
         tdtrrshape=tdtrrshape*(0.5*(zeff(lr)+zeff(lr+1))/ &
                    zeff(1))**difus_rshape(7)

      else

         tdtrrshape=(difus_rshape(1)+difus_rshape(2)* &
                    (rrz(lr)/radmin)**difus_rshape(3))**difus_rshape(4)
!        Simply use last (setup0%lrzmax) value of density....
         tdtrrshape=tdtrrshape* &
                    (reden(kk,lr)/reden(kk,0))**difus_rshape(5)
         tdtrrshape=tdtrrshape* &
                    (temp(kk,lr)/temp(kk,0))**difus_rshape(6)
         tdtrrshape=tdtrrshape* &
                    (zeff(lr)/zeff(1))**difus_rshape(7)

      endif

      return
      end function tdtrrshape
!
!
      real(c_double) function tdtrrshape1(lr)
      use param_mod
      use cqlcomm_mod, only : rya, difus_rshape
      implicit integer (i-n), real(c_double) (a-h,o-z)

!.......................................................................
!     This routine computes a radial shape function for
!     the radial diffusion coefficient.
!     Coefficients in the following expression are input
!     throught the namelist variable difus_rshape(1:7).
!.......................................................................

!
      if (lr.eq.0) then

         tdtrrshape1=0.

      elseif (lr.le.setup0%lrzmax-1) then

         ravg=0.5*(rya(lr)+rya(lr+1))

         if (ravg.le.difus_rshape(8)) then
            tdtrrshape1=1.
         elseif (ravg.le.(1.1*difus_rshape(8))) then
            tdtrrshape1=0.5*(1.+cos(pi*(ravg-difus_rshape(8))/ &
                                        (0.1*difus_rshape(8))))
         else
            tdtrrshape1=0.
         endif

      else

         ravg=rya(lr)

         if (ravg.le.difus_rshape(8)) then
            tdtrrshape1=1.
         elseif (ravg.le.1.1*difus_rshape(8)) then
           tdtrrshape1=0.5*(1.+cos(pi*(ravg-difus_rshape(8))/ &
                                       (0.1*difus_rshape(8))))
         else
            tdtrrshape1=0.
         endif

      endif

      return
      end function tdtrrshape1
!
!
      subroutine tdtrvshape(k,l)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!.......................................................................
!     This routine computes a velocity-space shape function
!     for the radial diffusion coefficient.
!     Coefficients in the following expression are input
!       throught the namelist variable difus_vshape(1:4).
!     The shape function is stored in temp1.
!     BH180921:
!     Note that apart from the coll_cutoff and gamm factor,
!     that the velocity dependence is in terms of velocity
!     (not momentum-per-mass) normalized to central thermal
!     velocity vth(k,1), with value 1.0 at the thermal vel.
!.......................................................................


      real(c_double) l_autocorr,lambda_mfp

      lr=setup0%lrindx(l)
      l_autocorr=pi*qsafety(lr)*radmaj

      call bcast(temp1(0:iy+1,0:jx+1),zero,iyjx2)
!      write(*,*)'tdtrdfus:difus_vshape',difus_vshape
      do  j=1,jx
         vel=x(j)*vnorm/gamma(j)
         coll=coll_freq(vel,k,lr)
         gamm=gamma(j)**difus_vshape(4)
         do  i=1,iytr(l)
            vpar=abs(vel*coss(idx(i,l),l))
            vprp=vel*sinn(idx(i,l),l)
            lambda_mfp=max(em100,vpar/coll)
            coll_cutoff=(1.+l_autocorr/lambda_mfp)**difus_vshape(2)
            vpar=vpar/vth(k,1)
            vprp=vprp/vth(k,1)
            vpar=vpar**difus_vshape(1)
            vprp=vprp**difus_vshape(3)
!      write(*,*)'tdtrvshape, l_autocorr,lambda_mfp,coll,vpar,vprp',
!     +                       l_autocorr,lambda_mfp,coll,vpar,vprp
            temp1(idx(i,l),j)=vpar/coll_cutoff*vprp/gamm
!      write(*,*)'tdtrvshape,idx(i,l),i,j,l,temp1(idx(i,l),j),coll_cut',
!     +                     idx(i,l),i,j,l,temp1(idx(i,l),j),coll_cutoff

         enddo
      enddo

      return
      end subroutine tdtrvshape




!
!
      real(c_double) function coll_freq(vel,k,lr)
      use param_mod
      use cqlcomm_mod, only : bnumb, zeff, temp, reden, fmass
      implicit integer (i-n), real(c_double) (a-h,o-z)

!.......................................................................
!     This routine computes
!     Approx. electron or ion collision freq from NRL
!     perp. deflection time.
!     Relativistic effects are neglected, since this
!     term is only significant in tdtrvshape for lambda_mfp
!     less than the magnetic fluctuation autocorrelation length,
!     generally for velocites much less that vth.
!.......................................................................


      real(c_double) lnlambda

!     Approx. electron or ion collision freq from NRL
!     perp. deflection time

      if (vel.eq.zero) then
         coll_freq=ep100
         return
      endif


      lnlambda=24.-log(reden(kelec,lr)**0.5/temp(kelec,lr))

      if (k.eq.kelecg) then

         energyev=0.5*fmass(kelec)*vel**2/ergtkev*1000.
         tempev=temp(kelec,lr)*1000.
         fnu_perp1=5.8e-6/(tempev**0.5*energyev)
         fnu_perp2=7.7e-6/energyev**1.5
         coll_freq1=min(fnu_perp1, fnu_perp2)*reden(k,lr)*lnlambda

         tempev=temp(kionn,lr)*1000.
         fmu=fmass(kionn)/proton
         fnu_perp1=2.5e-4*fmu**0.5/(tempev**0.5*energyev)
         fnu_perp2=7.7e-6/energyev**1.5
         coll_freq2=min(fnu_perp1, fnu_perp2)* &
                   reden(kelec,lr)*zeff(lr)*lnlambda

         coll_freq=coll_freq1+coll_freq2

      else                                           !ion case

         energyev=0.5*fmass(kionn)*vel**2/ergtkev*1000.
         tempev=temp(kionn,lr)*1000.
         fmu=fmass(kionn)/proton
         fnu_perp1=1.4e-7/(fmu**0.5*tempev**0.5*energyev)
         fnu_perp2=1.8e-7/(fmu**0.5*energyev**1.5)
         coll_freq=min(fnu_perp1, fnu_perp2)* &
                   (reden(kelec,lr)*zeff(lr)*bnumb(kionn)**2*lnlambda)

      endif

      return
      end function coll_freq


      subroutine ryaintorz(npts_in,oldx,oldf,npts,ynewx,ynewf)
      use param_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!.......................................................................
!     This routine interpolates difin vector in oldx, given in ryain
!     in cqlinput file to the transport mesh ynewx.
!     Boundary conditions are:
!                 1st derivative = 0 at plasma centre
!                 2nd derivative = 0 at plasma edge
!
!     At this stage uses linear interpolation
!.......................................................................

      dimension oldx(1),oldf(1),ynewx(1),ynewf(1)

      do jj = 1,npts_in-1
        do ll = 0,npts
          if (ynewx(ll).ge.oldx(jj) .and. ynewx(ll).le.oldx(jj+1)) then
            ynewf(ll) = &
              oldf(jj)+(oldf(jj+1)-oldf(jj))/(oldx(jj+1)-oldx(jj))* &
              (ynewx(ll)-oldx(jj))
          endif
        enddo
      enddo

      return
      end subroutine ryaintorz


      subroutine diffus_io(kopt)
      use param_mod
      use cqlcomm_mod, only : difus_io, d_rr, d_r, rya, rpconz
      use cqlcomm_mod, only : temp1, difus_io_file, t_, iy, iy_, y
      use cqlcomm_mod, only : tem1
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!.......................................................................
!     Write or read radial diffusion coeff and advective/pinch term,
!     and relate grids.
!     kopt=0:  Define output file & variables (called before kopt=1 or 2)
!     kopt=1:  write d_rr
!     kopt=2:  write d_rr,d_r
!     kopt=3:  read d_rr
!     kopt=4:  read d_rr,d_r
!.......................................................................

      include 'netcdf.inc'

! --- some stuff for netCDF file ---
      integer ncid,vid,istatus
      integer xdim,ydim,rdim,r0dim,gdim
      integer char64dim
      character*128 name,ltitle

      integer dims(3),count1(3),start1(3)  !for ngen=1 case
      integer dimsg(4),countg(4),startg(4)  !for ngen.gt.1 case
      integer y_dims(2),y_count(2)

      real(c_double), allocatable :: wkpack(:) ! local working array for pack21

      data start1/1,1,1/
      data startg/1,1,1,1/

      WRITE(*,*) 'In diffus_io, kopt=',kopt

!     Maximum iy as function of radius:
      iyy=0
      do l=1,setup0%lrz
         iyy=max(iyy,iy_(l))
      enddo
      if(iyy.gt.iy) stop 'difus_io: iy_(l) should not exceed iy'

      if (.NOT.ALLOCATED(wkpack)) then ! allocate working array for pack21
         nwkpack=iyjx2 +10 ! +10 just in case
         allocate(wkpack(nwkpack),STAT=istat)
         call bcast(wkpack,zero,SIZE(wkpack))
      endif

      if (kopt.eq.0) then  !Down to line 682

! ngen=1 case
      count1(1)=iymax
      count1(2)=jx
      count1(3)=1  !radial index

! ngen.gt.1 case
      countg(1)=iymax
      countg(2)=jx
      countg(3)=1  !radial index, 1 at at time
      countg(4)=1  !species index, 1 at at time

      y_count(1)=iymax
      y_count(2)=setup0%lrz  ! setup0%lrz radii at a time

!.......................................................................
!l    1.1.1 create netCDF filename.
!     CLOBBER old file, if it exists.
!     istatus is 0, if no errors.
!     ncid is created file id.

      if (difus_io_file.eq."setup0%mnemonic") then
         write(t_,1000) setup0%mnemonic(1:length_char(setup0%mnemonic))
 1000    format(a,"_difus_io.nc")
      else
!         write(t_,1001)
! 1001    format(a,difus_io_file)
         t_=difus_io_file
      endif

      istatus = NF_CREATE(t_, NF_CLOBBER, ncid)
      call check_err(istatus)

!     Define dimensions

      istatus= NF_DEF_DIM(ncid, 'xdim',     jx,       xdim)
      istatus= NF_DEF_DIM(ncid, 'ydim',     iymax,    ydim)
      istatus= NF_DEF_DIM(ncid, 'rdim',     setup0%lrz,      rdim)
      istatus= NF_DEF_DIM(ncid, 'r0dim', setup0%lrzmax, r0dim)
!     Number of general species treated by the _difus_io.nc file
      n_d_rr=0
      do k=1,ngen
         if (difus_io(k).ne."disabled") n_d_rr=n_d_rr+1
      enddo
      istatus= NF_DEF_DIM(ncid, 'drr_gen_species_dim', n_d_rr,  gdim)

      istatus= NF_DEF_DIM(ncid, 'char64dim',   64,    char64dim)

!     Define vectors of dimensions
      dims(1)=ydim
      dims(2)=xdim
      dims(3)=rdim

      dimsg(1)=ydim
      dimsg(2)=xdim
      dimsg(3)=rdim
      dimsg(4)=gdim

      y_dims(1)=ydim
      y_dims(2)=rdim


!     Define variables
!     Note, the variable IDs (denoted below as vid) are
!     not saved here in this subroutine; rather, the IDs
!     are retrieved from the netCDF data file, as needed,
!     by calling the netCDF routine ncvid.

      ltitle='NetCDF output/input of diffusion/pinch coeffs'
      if( length_char(ltitle).gt.128 ) stop 'Adjust ltitle in difus_io'
      call ncaptc2(ncid,NCGLOBAL,'title',NCCHAR,length_char(ltitle), &
           ltitle,istatus)

      vid=ncvdef0(ncid,'version',NCCHAR,1,char64dim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,20, &
                'CQL3D version number',istatus)

      vid=ncvdef0(ncid,'setup0%mnemonic',NCCHAR,1,char64dim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,23, &
                'Mnemonic run identifier',istatus)

!  Mesh related quantities

      vid=ncvdef0(ncid,'setup0%lrzmax',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,25, &
                  'Number of radial surfaces',istatus)
      call check_err(istatus)

      vid=ncvdef0(ncid,'rya',NCDOUBLE,1,r0dim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,22, &
                 'Normalized radial mesh',istatus)
      call check_err(istatus)

      vid=ncvdef0(ncid,'rpconz',NCDOUBLE,1,r0dim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,39, &
                 'Major radius at outside of flux surface',istatus)
      call check_err(istatus)

      vid=ncvdef0(ncid,'rhomax',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                 'cms',istatus)
      call check_err(istatus)

      vid=ncvdef0(ncid,'setup0%lrz',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,29, &
                  'Number of FPd radial surfaces',istatus)
      call check_err(istatus)

      vid=ncvdef0(ncid,'setup0%lrindx',NCLONG,1,rdim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,30, &
                 'Radial indices of FPd surfaces',istatus)

      vid=ncvdef0(ncid,'jx',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,27, &
                  'momentum-per-mass dimension',istatus)
      call check_err(istatus)

      vid=ncvdef0(ncid,'x',NCDOUBLE,1,xdim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,28, &
                 'normalized momentum-per-mass',istatus)
      call check_err(istatus)

      vid=ncvdef0(ncid,'vnorm',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,33, &
                 'velocity (momentum-per-mass) norm',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7, &
                           'cms/sec',istatus)
      call check_err(istatus)

      vid=ncvdef0(ncid,'enorm',NCDOUBLE,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,20, &
                           'Energy normalization',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,3, &
                           'keV',istatus)
      call check_err(istatus)

      vid=ncvdef0(ncid,'iy',NCLONG,0,0,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,25, &
                  'Max pitch angle dimension',istatus)
      call check_err(istatus)

      vid=ncvdef2(ncid,'y',NCDOUBLE,2,y_dims,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,11, &
                 'pitch angle',istatus)
      call ncaptc2(ncid,vid,'units',NCCHAR,7, &
                 'radians',istatus)
      call check_err(istatus)

      vid=ncvdef0(ncid,'iy_',NCLONG,1,rdim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,36, &
                  'Pitch angle dimension at each radius',istatus)
      call check_err(istatus)

      vid=ncvdef0(ncid,'itl',NCLONG,1,rdim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,26, &
                 'lower trapped-passing bndy',istatus)

      vid=ncvdef0(ncid,'itu',NCLONG,1,rdim,istatus)
      call ncaptc2(ncid,vid,'long_name',NCCHAR,26, &
                 'upper trapped-passing bndy',istatus)
      call check_err(istatus)

!  d_rr diffusion coeff, and d_r pinch (kopt.eq.2)
!  Could use ngen.gt.1 coding for ngen=1 case, if turns out preferable.

!  n_d_rr is the number of diffused general species

      if (n_d_rr.eq.1) then
         vid=ncvdef2(ncid,'d_rr',NCDOUBLE,3,dims,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,27, &
              'radial diffusion coefficient',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,8, &
              'cm**2/sec',istatus)
         !Define d_r pinch velocity if difus_io(1).eq."drrdrout"
         if (difus_io(1).eq."drrdrout") then
            vid=ncvdef2(ncid,'d_r',NCDOUBLE,3,dims,istatus)
            call ncaptc2(ncid,vid,'long_name',NCCHAR,38, &
                 'radial vel (pinch) term (pos, outwards)',istatus)
            call ncaptc2(ncid,vid,'units',NCCHAR,6, &
                 'cm/sec',istatus)
         endif
      endif

      if (n_d_rr.gt.1) then
         vid=ncvdef2(ncid,'d_rr',NCDOUBLE,4,dimsg,istatus)
         call ncaptc2(ncid,vid,'long_name',NCCHAR,27, &
              'radial diffusion coefficient',istatus)
         call ncaptc2(ncid,vid,'units',NCCHAR,8, &
              'cm**2/sec',istatus)
         !Define d_r pinch velocity if difus_io(1).eq."drrdrout"
         if (difus_io(1).eq."drrdrout") then
            vid=ncvdef2(ncid,'d_r',NCDOUBLE,4,dimsg,istatus)
            call ncaptc2(ncid,vid,'long_name',NCCHAR,38, &
                 'radial vel (pinch) term (pos, outwards)',istatus)
            call ncaptc2(ncid,vid,'units',NCCHAR,6, &
                 'cm/sec',istatus)
         endif
      endif

!  End define mode
      istatus= NF_ENDDEF(ncid)

      endif  !On kopt.eq.0

!  Close file

      istatus= NF_CLOSE(ncid)



!.......................................................................
!  Write data  (assuming prior call with kopt=0)
!.......................................................................

      if (kopt.eq.1 .or. kopt.eq.2) then   !Down to 818

      istatus=NF_OPEN(t_, NF_WRITE, ncid)

!  Write version, setup0%mnemonic,and grid related quantities

      istatus= NF_INQ_VARID(ncid,'version',vid)
      ll=length_char(version)
      call ncvptc0(ncid,vid,1,ll,version,ll,istatus)

      istatus= NF_INQ_VARID(ncid,'setup0%mnemonic',vid)
      ll=length_char(setup0%mnemonic)
      call ncvptc0(ncid,vid,1,ll,setup0%mnemonic,ll,istatus)


      istatus= NF_INQ_VARID(ncid,'setup0%lrzmax',vid)
      istatus = NF_PUT_VARA_INT(ncid,vid,1,1,setup0%lrzmax)

      istatus= NF_INQ_VARID(ncid,'rya',vid)
      istatus = NF_PUT_VARA_DOUBLE(ncid,vid,(1),setup0%lrzmax,rya(1))

      istatus= NF_INQ_VARID(ncid,'rpconz',vid)
      istatus = NF_PUT_VARA_DOUBLE(ncid,vid,(1),setup0%lrzmax,rpconz(1))

      istatus= NF_INQ_VARID(ncid,'rhomax',vid)
      istatus = NF_PUT_VARA_DOUBLE(ncid,vid,(1),1,(rhomax))

      istatus= NF_INQ_VARID(ncid,'setup0%lrz',vid)
      istatus = NF_PUT_VARA_INT(ncid,vid,1,1,setup0%lrz)

      istatus= NF_INQ_VARID(ncid,'setup0%lrindx',vid)
      istatus = NF_PUT_VARA_INT(ncid,vid,1,setup0%lrz,setup0%lrindx(1))

      istatus= NF_INQ_VARID(ncid,'jx',vid)
      istatus = NF_PUT_VARA_INT(ncid,vid,1,1,(jx))

      istatus= NF_INQ_VARID(ncid,'x',vid)
      istatus = NF_PUT_VARA_DOUBLE(ncid,vid,(1),jx,x)

      istatus= NF_INQ_VARID(ncid,'vnorm',vid)
      istatus = NF_PUT_VARA_DOUBLE(ncid,vid,(1),1,(vnorm))

      istatus= NF_INQ_VARID(ncid,'enorm',vid)
      istatus = NF_PUT_VARA_DOUBLE(ncid,vid,(1),1,(enorm))

      istatus= NF_INQ_VARID(ncid,'iy',vid)
      istatus = NF_PUT_VARA_INT(ncid,vid,1,1,(iymax))

      if (iy*lrors.gt.iyjx2) stop 'netcdfrf:  Need to set jx>lrza'
      call pack21(y,1,iy,1,lrors,tem1,iymax,lrors)
      istatus= NF_INQ_VARID(ncid,'y',vid)
      istatus = NF_PUT_VARA_DOUBLE(ncid,vid,start1,y_count,tem1)

      istatus= NF_INQ_VARID(ncid,'iy_',vid)
      istatus = NF_PUT_VARA_INT(ncid,vid,1,setup0%lrz,iy_)

      istatus= NF_INQ_VARID(ncid,'itl',vid)
      istatus = NF_PUT_VARA_INT(ncid,vid,1,setup0%lrz,itl_)

      istatus= NF_INQ_VARID(ncid,'itu',vid)
      istatus = NF_PUT_VARA_INT(ncid,vid,1,setup0%lrz,itu_)

!  n_d_rr is the number of diffused general species
      if (n_d_rr.eq.1) then

      istatus= NF_INQ_VARID(ncid,'d_rr',vid)
         do ll=1,setup0%lrz
            do j=1,jx
               do i=1,iy
                  temp1(i,j)=d_rr(i,j,1,setup0%lrindx(ll))
               enddo
            enddo
            call pack21(temp1,0,iyp1,0,jxp1,wkpack,iy,jx)
            start1(3)=ll
            istatus = NF_PUT_VARA_DOUBLE(ncid,vid,start1,count1,wkpack)
         enddo
         if (kopt.eq.2) then
         istatus= NF_INQ_VARID(ncid,'d_r',vid)
         do ll=1,setup0%lrz
            do j=1,jx
               do i=1,iy
                  temp1(i,j)=d_r(i,j,1,setup0%lrindx(ll))
               enddo
            enddo
            call pack21(temp1,0,iyp1,0,jxp1,wkpack,iy,jx)
            start1(3)=ll
            istatus = NF_PUT_VARA_DOUBLE(ncid,vid,start1,count1,wkpack)
         enddo
         endif  !On kopt.eq.2

      else                      !n_d_rr.ge.2

         istatus= NF_INQ_VARID(ncid,'d_rr',vid)
         do k=1,n_d_rr
            do ll=1,setup0%lrz
               do j=1,jx
                  do i=1,iy
                     temp1(i,j)=d_rr(i,j,k,setup0%lrindx(ll))
                  enddo
               enddo
               call pack21(temp1,0,iyp1,0,jxp1,wkpack,iy,jx)
               startg(3)=ll
               startg(4)=k
               istatus = NF_PUT_VARA_DOUBLE(ncid,vid,startg,countg,wkpack)
            enddo               !  On ll
         enddo                  !  On k
         if (kopt.eq.2) then
         istatus= NF_INQ_VARID(ncid,'d_r',vid)
         do k=1,n_d_rr
            do ll=1,setup0%lrz
               do j=1,jx
                  do i=1,iy
                     temp1(i,j)=d_r(i,j,k,setup0%lrindx(ll))
                  enddo
               enddo
               call pack21(temp1,0,iyp1,0,jxp1,wkpack,iy,jx)
               startg(3)=ll
               startg(4)=k
               istatus = NF_PUT_VARA_DOUBLE(ncid,vid,startg,countg,wkpack)
            enddo               !  On ll
         enddo                  !  On k
         endif  !On kopt.eq.2

      endif                     ! on n_d_rr

      endif  !On kopt.eq.1 .or. kopt.eq.2

!.......................................................................
!  Read diffusion/pinch terms
!  It is not necessary to have prior call to diffus_io
!.......................................................................

      if (kopt.eq.3 .or. kopt.eq.4) then

      if (difus_io_file.eq."setup0%mnemonic") then
         t_=setup0%mnemonic//"_difus_io.nc"
      else
         t_=difus_io_file
      endif

      istatus=NF_OPEN(t_,NF_NOWRITE,ncid)
      if (istatus .NE. NF_NOERR) then
         WRITE(*,*)'   ***   Problem opening d_rr .nc data file   ***'
         Stop
      endif

!  Checking dimensions
      istatus= NF_INQ_DIMID(ncid,'xdim',xdim)
      istatus= NF_INQ_DIMID(ncid,'ydim',ydim)
      istatus= NF_INQ_DIMID(ncid,'rdim',rdim)
      istatus= NF_INQ_DIMID(ncid,'drr_gen_species_dim',gdim)
      istatus= NF_INQ_DIM(ncid,xdim,name,jx_file) !name should ='xdim', etc.
      istatus= NF_INQ_DIM(ncid,ydim,name,iy_file)
      istatus= NF_INQ_DIM(ncid,rdim,name,lrz_file)
      istatus= NF_INQ_DIM(ncid,gdim,name,n_d_rr_file)

!     Number of general species treated by the _difus_io.nc file
      n_d_rr=0
      do k=1,ngen
         if (difus_io(k).ne."disabled") n_d_rr=n_d_rr+1
      enddo

!     Check the data dimensions in the input file

      if (jx*iy*setup0%lrz*n_d_rr.ne.jx_file*iy_file*lrz_file*n_d_rr_file) then
         write(*,*)
         write(*,*)'  WRONG DIMENSIONS IN _difus_io.nc file: STOP'
         STOP
      endif

!     Could check rya, etc., in the _difus_io.nc file (but not yet done).

!     Read d_rr and (kopt.eq.4) d_r

      if (n_d_rr.eq.1) then

!     start1=10  !Testing that sets all start1(1:3).  (Yes, 4 gfortran.)
      start1(1)=1
      start1(2)=1
      start1(3)=1
      count1(1)=iy
      count1(2)=jx
      count1(3)=1

      istatus= NF_INQ_VARID(ncid,'d_rr',vid)
      do ll=1,setup0%lrz
         start1(3)=ll
         istatus= NF_GET_VARA_DOUBLE(ncid,vid,start1,count1,wkpack)
         call unpack21(temp1,0,iyp1,0,jxp1,wkpack,iy,jx)
         do j=1,jx
            do i=1,iy
               d_rr(i,j,1,setup0%lrindx(ll))=temp1(i,j)
            enddo
         enddo
      enddo  !On ll

      if (kopt.eq.4) then
      istatus= NF_INQ_VARID(ncid,'d_r',vid)
      do ll=1,setup0%lrz
         start1(3)=ll
         istatus= NF_GET_VARA_DOUBLE(ncid,vid,start1,count1,wkpack)
         call unpack21(temp1,0,iyp1,0,jxp1,wkpack,iy,jx)
         do j=1,jx
            do i=1,iy
               d_r(i,j,1,setup0%lrindx(ll))=temp1(i,j)
            enddo
         enddo
      enddo  !On ll
      endif  !On kopt.eq.4

      else  !n_d_rr.gt.1

      startg=1  !should give startg(1:4)=1
      countg(1)=iy
      countg(2)=jx
      countg(3)=1
      countg(4)=1

      istatus= NF_INQ_VARID(ncid,'d_rr',vid)
      do k=1,n_d_rr
         startg(4)=k
         do ll=1,setup0%lrz
            startg(3)=ll
            istatus= NF_GET_VARA_DOUBLE(ncid,vid,startg,countg,wkpack)
            call unpack21(temp1,0,iyp1,0,jxp1,wkpack,iy,jx)
            do j=1,jx
               do i=1,iy
                  d_rr(i,j,1,setup0%lrindx(ll))=temp1(i,j)
               enddo
            enddo
         enddo
      enddo  !On k

      if (kopt.eq.4) then
      istatus= NF_INQ_VARID(ncid,'d_r',vid)
      do k=1,n_d_rr
         startg(4)=k
         do ll=1,setup0%lrz
            startg(3)=ll
            istatus= NF_GET_VARA_DOUBLE(ncid,vid,startg,countg,wkpack)
            call unpack21(temp1,0,iyp1,0,jxp1,wkpack,iy,jx)
            do j=1,jx
               do i=1,iy
                  d_r(i,j,1,setup0%lrindx(ll))=temp1(i,j)
               enddo
            enddo
         enddo
      enddo  !On k
      endif  !On kopt.eq.4

      endif  !On n_d_rr

      istatus = NF_CLOSE(ncid)
      call check_err(istatus)

      start1(3)=1  !For safety
      startg(3)=1
      startg(4)=1

      endif  !On kopt.eq.3 .or. kopt.eq.4

      return
      end subroutine diffus_io
!
!
      real(c_double) function difus_io_scale(k,iopt)
      use param_mod
      use cqlcomm_mod, only : difus_io_drrscale, difus_io_drscale
      use cqlcomm_mod, only : difus_io_t
      implicit integer (i-n), real(c_double) (a-h,o-z)

!......................................................................
!     Returns t-dependent scale factor for radial diffusion drr
!     or pinch vel dr, for each species k.
!     iopt=1,   for drr
!         =2,   for dr
!......................................................................

      if (iopt.eq.1) then

      if (ndifus_io_t.le.0) then
         difus_io_scale=one
      else
         itme=1
         do jtm=1,ndifus_io_t
            if (timet.ge.difus_io_t(jtm)) itme=jtm
         enddo
         itme1=itme+1
         if (itme.lt.ndifus_io_t) then
            difus_io_scale=difus_io_drrscale(itme,k) &
                 +(difus_io_drrscale(itme1,k)-difus_io_drrscale(itme,k)) &
                 /(difus_io_t(itme1)-difus_io_t(itme)) &
                 *(timet-difus_io_t(itme))
         else
            difus_io_scale=difus_io_drrscale(ndifus_io_t,k)
         endif
      endif  !  On ndifus_io_t

      elseif (iopt.eq.2) then

      if (ndifus_io_t.le.0) then
         difus_io_scale=one
      else
         itme=1
         do jtm=1,ndifus_io_t
            if (timet.ge.difus_io_t(jtm)) itme=jtm
         enddo
         itme1=itme+1
         if (itme.lt.ndifus_io_t) then
            difus_io_scale=difus_io_drscale(itme,k) &
                 +(difus_io_drscale(itme1,k)-difus_io_drscale(itme,k)) &
                 /(difus_io_t(itme1)-difus_io_t(itme)) &
                 *(timet-difus_io_t(itme))
         else
            difus_io_scale=difus_io_drscale(ndifus_io_t,k)
         endif
      endif  !  On ndifus_io_t

      endif  !  On iopt

      return

      end function difus_io_scale


end module tdtrdfus_mod
