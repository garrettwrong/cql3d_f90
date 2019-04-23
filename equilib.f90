!
!
module equilib_mod

  !---BEGIN USE

  use eqtopeol_mod, only : eqtopeol
  use eqwrng_mod, only : eqwrng
  use r8subs_mod, only : dcopy
  use zcunix_mod, only : coeff1
  use zcunix_mod, only : coeff2
  use zcunix_mod, only : terp1e
  use zcunix_mod, only : t2

  !---END USE
  use iso_c_binding, only : c_double
  use r8subs_mod, only : dcopy
  ! these are used outside the module in tdeqdisk
  integer, parameter, public ::  ncoila=5, nccoila=50
  real(c_double), public :: pcvac(0:9)
  real(c_double), public :: ccoil(nccoila,ncoila)
  integer,  public :: ncoil(ncoila)
  save
  
contains
  
  subroutine equilib(Requil,Zequil,index,PSIequil,BReq,BPHIeq,BZeq)
    use param_mod
    use comm_mod
    implicit integer (i-n), real*8 (a-h,o-z)
    
    !  Call with index=0 for setup,
    !            index=1 for PSI,BR,BPHI,BZ return (cgs units)
    !                    for a given Requil,Zequil coordinate.
    
    
    parameter(nworka=3*nnra+1)
    
    character*8 ntitle,dat   ! Added for g95 compiler, Urban 110708
    dimension ntitle(5),workk(nworka)
    !..................................................................
    !     CQL3D subroutine:  Oct. 27, 2004, RWH.
    !     This routine will read in the eqdsk file if eqsource="eqdsk".
    !     The eqdsk data is up-down symmetrized if so speciefied
    !     by the input variable eqsym.
    !     The cql3d code has been generalized for non-up-down symmertic
    !     equilibria (BH: Oct'09).
    !     The routine assumes that the data coming from the equilibrium
    !     code is in MKS. Conversion is made to cgs. Finally, the setup
    !     for the routine which generates f(psi) (called eqfpsi) is
    !     done. This involves a call to a spline setup routine.
    !     If eqsource="topeol", eqtopeol is called.
    !     If eqsource="tsc", data has already been read in from file
    !     tscinp. It is a matter only of changing data to cgs.
    !..................................................................
    
    
    !..................................................................
    !   EQDSK Description
    !..................................................................
    ! In the eqdsk, poloidal flux function psi [poloidal flux/(2*pi)] is
    ! given on a rectagular grid, equispaced in the horizontal (r)
    ! and vertical direction.  Other parameters of the equilibrium
    ! are given, as described in the following.  Units are entirely MKSA.
    !
    ! A standard convention of EFIT (a source of eqdsk files) is that
    ! the psi-funtion is minimum at the magnetic axis, regardless of
    ! the direction of the toroidal current.  This convention is
    ! not followed by all eqdsk generators.
    !
    ! nnr,nnz,nnv      are the numbers of points in the horizontal (r)
    !                  and vertical (z) directions.  In some cases,
    !                  nnv, a number of radial flux surfaces different
    !                  from nnr, MAY be given.
    ! rbox,zbox        are the horizontal and vertical full-widths
    !                  of the  rectangle on which poloidal flux values
    !                  are given.
    ! radmaj           is the nominal major radius of the torus.
    ! rboxdst          is the major radius of the inner edge of the box.
    !                  of rectangular grid.
    ! ymideqd          is the vertical shift of the rectangular box
    !                  midplane.
    ! raxis,zaxis      are the major and vertical height of magnetic axis.
    ! psimag,psilim    are the poloidal flux function values at the
    !                  magnetic axis and the last closed flux surface
    !                   (touching the limiter or the separatrix).
    ! btor             is the vacuum toroidal magnetic field at radmaj.
    ! toteqd           is the toroidal current.
    ! psimx1,psixm2,xax1,xax2,zax1,zax2,psisep,xsep,ysep - OBSOLETE.
    !   [These quantities have been used in the description of doublet
    !    plasmas, with two magnetic axes.]
    !  ATTENTION:
    ! psimx1=psimag
    ! xax1=xma
    ! zax1=zma
    ! psisep=psilim
    !
    ! feqd(nnr),pres(nnr), usually, although dimension is
    !                      specified by nnv if it is present.
    ! fpsiar = r * B_phi  at nnr equispaced points in psi
    !          from psimag to psilim.
    ! prar   = are the pressure values at the same points.
    ! ffpar  = fpsiar_prime, [i.e., d f/d (psi)] at the same points.
    ! ppeqd  =  p_prime [i.e., d prar/d (psi)] at the same points.
    ! epsi(nnr,nnz) are the psi values on the nnr * nnz
    !               equispaced grid.
    ! qar(nnr)  gives safety factor q on the equispaced psi grid.
    
    ! The following quantities are given in some eqdsks, but are not
    ! presently used in cql3d:
    ! nlimit,nves  are the numbers of point at limiters and
    !              vacuum vessel wall.
    ! rlimit(nlimit),zlimit(nlimit):
    ! rlimit,zlimit      is the r,z location of the limiters wall.
    !       rves(nves),zves(nves):
    ! rves,zves          is the r,z location of the limiters
    !                    vacuum vessel wall.
    !
    ! Additional plasma data and coil data appears after the above
    ! data is some eqdsks.  In general, eqdsks need to contain
    ! at least all the above data down through epsi.  Most also
    ! contain the qar(), although is can be derived from the
    ! preceding data.
    !..................................................................
    

    if (index .eq. 0) then !==========================================
       
       if(lr_.ne.lrzmax) return ! equilib is called from aingeom in lr loop
       ! starting at lr_=lrzmax. Only calc once.
       
       if (eqsource.eq."tsc") then
          do j=1,nnz
             do i=1,nnr
                epsi(i,j)=epsi(i,j)-psilim
             end do
          end do
        do 170 i=1,nnv
          psiar(i)=-(psiar_(nnv+1-i)-psilim)*1.e+8
          fpsiar(i)=fpsiar_(i)
170    end do
       psimag=psimag-psilim
       psilim=0.
       rbox=rright-rleft
       zbox=ztop_-zbot_
       rboxdst=rleft
       nfp=nnv
    else if (eqsource.eq."eqdsk") then

       !BH020822  Adding nveqd (.le. nxeqd) for different number of
       !BH020822  flux surfaces on which p,feqd,p',feqd',q are tabulated.
       !BH020822  The standard EQDSK does not incorporate this feature,
       !BH020822  although is is available in cql3d.
       !BH020822  Bonoli uses it for ACCOME eqdsk output.
       
       open(unit=10,file=eqdskin,status='old')
       read(10,110,iostat=io1)(ntitle(i),i=1,5),dat,ipestg,nnr,nnz,nnv
       if (io1.ne.0) then
          rewind 30
          read(10,110,iostat=io1)(ntitle(i),i=1,5),dat,ipestg,nnr,nnz
          nnv=0
       endif
       if (nnv.eq.0) nnv=nnr
       if (nnv.gt.nnr) call eqwrng(11)
       if(nnr.gt.nnra .or. nnz.gt.nnza) call eqwrng(12)
       read(10,120) rbox,zbox,radmaj,rboxdst,ymideqd
       read(10,120) raxis,zaxis,psimag,psilim,btor
       read(10,120) toteqd,psimx1,psimx2,xax1,xax2
       read(10,120) zax1,zax2,psisep,xsep,zsep
       read(10,120) (fpsiar(i),i=1,nnv)
       
       read(10,120) (prar(i),i=1,nnv)
       read(10,120) (ffpar(i),i=1,nnv)
       read(10,120) (ppar(i),i=1,nnv)
       read(10,120) ((epsi(i,j),i=1,nnr),j=1,nnz)
       read(10,120) (qar(i),i=1,nnv)

       read (10, 8210)   ncontr, nlimiter
8210   format (2i5)
       !write(*,*)'ncontr, nconmax ',ncontr, nlimiter
       
       if (ncontr.gt.0) then  !-> LCFS
          allocate(rcontr(ncontr),STAT=istat) 
          allocate(zcontr(ncontr),STAT=istat) 
          read (10, 8200) (rcontr(i), zcontr(i),i = 1,ncontr) ! [m]
       else 
          !-YuP: Set LCFS = outer boundary of equilibrium grid
          ncontr=5
          allocate(rcontr(ncontr),STAT=istat) 
          allocate(zcontr(ncontr),STAT=istat) 
          rcontr(1) = rboxdst
          zcontr(1) = -0.5*zbox
          rcontr(2) = rboxdst
          zcontr(2) = +0.5*zbox
          rcontr(3) = rboxdst+rbox
          zcontr(3) = +0.5*zbox
          rcontr(4) = rboxdst+rbox
          zcontr(4) = -0.5*zbox
          rcontr(5) = rcontr(1)
          zcontr(5) = zcontr(1)
       endif
         
       if (nlimiter.gt.0) then !-> Limiter
          allocate(rlimiter(nlimiter),STAT=istat) 
          allocate(zlimiter(nlimiter),STAT=istat) 
          read (10, 8200) (rlimiter(i), zlimiter(i), i = 1,nlimiter) ![m]
       else 
          !-YuP: Set limiter = outer boundary of equilibrium grid
          nlimiter=5
          allocate(rlimiter(nlimiter),STAT=istat) 
          allocate(zlimiter(nlimiter),STAT=istat) 
          rlimiter(1) = rboxdst
          zlimiter(1) = -0.5*zbox
          rlimiter(2) = rboxdst
          zlimiter(2) = +0.5*zbox
          rlimiter(3) = rboxdst+rbox
          zlimiter(3) = +0.5*zbox
          rlimiter(4) = rboxdst+rbox
          zlimiter(4) = -0.5*zbox
          rlimiter(5) = rlimiter(1)
          zlimiter(5) = zlimiter(1)
       endif
         
8200   format(5e16.9)
110    format(6a8,4i4)
120    format(5e16.9)
250    format( (5(e21.14)) )
251    format(5i5)
       
       write(*,*)
       write(*,98) raxis,zaxis,psimag,psilim,qar(1),qar(nnv)
       write(*,99) rbox,zbox,rboxdst,ymideqd
98     format('ORIGINAL Eqdsk values [mks units]:  &
       raxis,zaxis,psimag,psilim,qar(1),qar(nnv)=',  6(1pe12.4))
99     format('Eqdsk values: rbox,zbox,rboxdst,ymideqd=',4(1pe12.4))
       write(*,*)

       !.......................................................................
       !     CQL cannot work correctly with B<0 (Lin-Liu)
       !     A namelist variable, bsign, has been added to
       !     adjust the toroidal field equilibrium data, and
       !     to keep track of the sign for RF QL diffusion coeff purposes.
       !.......................................................................
       if (btor.ge.0. .and. bsign.lt.0.) then
          stop "BTOR .gt.0, but bsign.lt.0. Check BSIGN"
       endif
       if (btor.lt.0.) then
          if (rdcmod.ne."disabled") then
             write(*,*) " " 
             write(*,*) " " 
             write(*,*) "   *******************************************"
             write(*,*) "   **BTOR .lt.0, Check RDC_UPAR_SIGN,      ***"
             write(*,*) "   **when using rdcmod.  DC coeffs use  -1.***"
             write(*,*) "   *******************************************"
             write(*,*) " " 
             write(*,*) " "
          endif
          if ( bsign.lt.0. ) then
             btor=bsign * btor
             do ii=1,nnv
                fpsiar(ii)=bsign * fpsiar(ii)
             enddo
             write(*,*)'equilib: Sign of btor and fpsiar is reversed'
          else
             write(*,*) " " 
             write(*,*) " " 
             write(*,*) "   *******************************************"
             write(*,*) "   ****BTOR .lt.0, Check BSIGN if using urf***"
             write(*,*) "   *******************************************"
             write(*,*) " " 
             write(*,*) " " 
             stop "BTOR .lt.0, Check BSIGN"
          endif
       endif
       !
       !$$$C%OS  
       !$$$C%OS  not used and creates error when ncoil not defined properly
       !$$$        read(10,251,END=299) (ncoil(i),i=1,5)
       !$$$        do 230 i=1,5
       !$$$          if(ncoil(i).le.0) go to 230
       !$$$C%OS  if(ncoil(i).gt.nccoila) call eqwrng(100)
       !$$$          if(ncoil(i).gt.nccoila) go to 299
       !$$$          read(10,250) (ccoil(nn,i),nn=1,ncoil(i))
       !$$$ 230    continue
       !$$$        read(10,250,END=299) (cvac(i),i=1,9)
       !$$$ 299    continue
       close(unit=10)


       !.......................................................................
       !     The convention we expect for the read in epsi from eqdsk
       !     within the plasma
       !     is that it will have a minimum at the magnetic axis.
       !     This is in conformity with usual eqdsk convention. 
       !     If psimag.gt.psilim (from an off-brand eqdsk), 
       !         we reverse the sign of epsi,psilim,psimag. 
       !     The sign of the current, toteqd, will be used
       !     to specify the current  and poloidal B-field directions.
       !     (Positive is CCW, viewed from above.)
       !     Keep sign of current in cursign.
       !
       !     HOWEVER, at the end of this subroutine, and throughout
       !     cql3d, we use epsi with reversed sign, i.e., epsi is maximum
       !     at the magnetic axis.  The bi-cubic splines of epsi use
       !     this latter convention.
       !
       !     Also, absolute value of fpsi and btor is effectively taken,
       !     and the sign of Bphi is maintained through bsign.
       !
       !     We have:
       !        vector-B = 
       !        bsign*fpsi_cql*grad(phi) - cursign*grad(phi) X grad(epsi_cql)
       !
       !.......................................................................
       
       if (psimag.gt.psilim) then ! reverse to make minimum at m.axis
          do j=1,nnz
             do i=1,nnr
                epsi(i,j)=-epsi(i,j)
             enddo
          enddo
          psimag=-psimag
          psilim=-psilim ! now psimag < psilim
          
          write(*,1000)
1000      format(//,1x,'WARNING: Sign of epsi,psilim,psimag reversed')
           
       endif

       cursign=sign(one,toteqd)
         
       !..................................................................
       !     Create the r,z meshes..
       !     Assumes ymideqd=0.
       !..................................................................
       
       !if(ymideqd.ne.zero) stop 'equilib: Consider effect ymideqd.ne.0'
       !YuP[04-2017] moved the above stop to few lines below.
       ! In case of eqsym=avg_zmag (or top, or bottom) the value of ymideqd
       ! will be redefined, so the code should not be stopped here.
       
       dzz=zbox/(nnz-1)
       drr=rbox/(nnr-1)
       er(1)=rboxdst
       ez(1)=-zbox*.5
       do nn=2,nnr
          er(nn)=er(nn-1)+drr
       enddo
       do nn=2,nnz
          ez(nn)=ez(nn-1)+dzz
       enddo
       
       !.......................................................................
       !     Up-down symmetrize the eqdsk about z=0, as specified by  eqsym. 
       !     (this has no effect if the equilibrium is initially 
       !     up-down symmetric about zmag=0.).
       !     New options beyond eqsym.eq."average" added (BobH: 020606).
       !.......................................................................
       
       if (eqsym.eq."none") then
          if(ymideqd.ne.zero)  stop 'equilib/eqsym=none: Consider effect ymideqd.ne.0'
          zshift=0.d0
          !          Do nothing.
       else if (eqsym.eq."average") then
          if(ymideqd.ne.zero) stop 'equilib/eqsym=average: Consider effect ymideqd.ne.0'
          zaxis=0.d0 ! YuP[2014] added
          !          Average psi values above and below z=0. This is nmlst deflt.
          do j=1,nnz/2
             do i=1,nnr
                epsi(i,j)=0.5*(epsi(i,j)+epsi(i,nnz-(j-1)))
             enddo
          enddo
          do j=nnz-nnz/2+1,nnz
             do i=1,nnr
                epsi(i,j)=epsi(i,nnz+1-j)
             enddo
          enddo
       else                 ! eqsym.eq.(avg_zmag, top, or bottom)
          !          Symmetrizing step:
          !
          !          Expand vertical height of eqdsk by amount 2*(zaxis-ymideqd),
          !          making midplane of computational grid at the magnetic axis.
          !    
          !          Interpolate psi onto new expanded grid.
          
          zshift=zaxis-ymideqd
          zaxis=0.d0
          ymideqd=0.d0
          if(ncontr.gt.1) then ! YuP[04-2017] Added: adjust LCFS, if available
             ncontr_= min(ncontr,lfielda)
             do ilim=1,ncontr_
                zcontr(ilim)=zcontr(ilim)-zshift
             enddo
          endif
          if(nlimiter.gt.1) then ! YuP[04-2017] Added: adjust plasma border
             nline= min(nlimiter,lfielda)
             do ilim=1,nline
                zlimiter(ilim)=zlimiter(ilim)-zshift
             enddo
          endif
          zdim_ex=zbox+2.*abs(zshift)
          dr_ex=drr
          dz_ex=zdim_ex/(nnz-1)
          if (zshift.ge.0.d0) then
             dummyaz(1)=ez(1)
          else
             dummyaz(1)=ez(1)-2.*abs(zshift)
          endif
          do ix=1,nnr
             dummyar(ix)=er(ix)
          enddo
          !          Gives a shifted, expanded z-array:
          do iz=2,nnz
             dummyaz(iz)=dummyaz(iz-1)+dz_ex
          enddo
           
          !          Interpolating:
          kz1=1
          kz2=2
          do 350  j=1,nnz
             zval=dummyaz(j)
360          if (zval.le.ez(kz2).or.kz2.eq.nnz) go to 365
             
             kz1=kz1+1
             kz2=kz1+1
             go to 360
             !     
365       continue
          kr1=1
          kr2=2
          do 370  i=1,nnr
             rval=dummyar(i)
380          if (rval.le.er(kr2).or.kr2.eq.nnr) go to 385
             
             kr1=kr1+1
             kr2=kr1+1
             go to 380

385       continue
          f1eq=epsi(kr1,kz1)+(rval-er(kr1))*(epsi(kr2,kz1)-epsi(kr1,kz1))/(er(kr2)-er(kr1))
          f2eq=epsi(kr1,kz2)+(rval-er(kr1))*(epsi(kr2,kz2)-epsi(kr1,kz2))/(er(kr2)-er(kr1))
          !     
          val=f1eq+(zval-ez(kz1))*(f2eq-f1eq)/(ez(kz2)-ez(kz1))
          dummypsi(i,j)=val
370    continue
350    continue
                
          !BH031010:  Following changes do nothing, since xma_ex,rdim_ex
          !BH031010:    not used.
          !BH031010           xma_ex=xma
          xma_ex=raxis
          !BH031010           rdim_ex=xdimeqd
          rdim_ex=rbox
          zmid_ex=0.d0
          zma_ex=0.d0
          
          !     Redefine the dummyaz-grid so midplane at z=0.:
          dummyaz(1)=-zdim_ex/2.
          do i=2,nnz
             dummyaz(i)=dummyaz(i-1)+dz_ex
          enddo
          
          !          Redefine the eqdsk quantities:
          zbox=zdim_ex
          zaxis=0.d0
          !BH031010           dzz=dz_ez  !No effect, since not used here
          !BH031010                      !and it is redefined below.
          do i=1,nnz
             ez(i)=dummyaz(i)
          enddo
          
          if (eqsym.eq."avg_zmag") then
             !          Up-down symmetrize
             do i=1,nnr
                j=nnz/2+1
                epsi(i,j)=dummypsi(i,j)
             enddo
             do j=1,nnz/2
                do i=1,nnr
                   epsi(i,j)=0.5*(dummypsi(i,j)+dummypsi(i,nnz-(j-1)))
                enddo
             enddo
             do j=nnz-nnz/2+1,nnz
                do i=1,nnr
                   epsi(i,j)=epsi(i,nnz+1-j)
                enddo
             enddo
             
          else if (eqsym.eq."top") then
             STOP 'Need to fix eqsym=top'
             !          Reflect top to bottom    
             do i=1,nnr
                j=nnz/2+1
                epsi(i,j)=dummypsi(i,j)
             enddo
             do j=nnz-nnz/2+1,nnz
                do i=1,nnr
                   epsi(i,j)=dummypsi(i,j)
                enddo
             enddo
             do j=1,nnz/2
                do i=1,nnr
                   epsi(i,j)=dummypsi(i,nnz-(j-1))
                enddo
             enddo
             
          else if (eqsym.eq."bottom") then
             STOP 'Need to fix eqsym=bottom'
             !          Reflect bottom to top
             do i=1,nnr
                j=nnz/2+1
                epsi(i,j)=dummypsi(i,j)
             enddo
             do j=1,nnz/2
                do i=1,nnr
                   epsi(i,j)=dummypsi(i,j)
                enddo
             enddo
             do j=nnz-nnz/2+1,nnz
                do i=1,nnr
                   epsi(i,j)=dummypsi(i,nnz-(j-1))
                enddo
             enddo
             
          endif
          
       endif  !On eqsym
        
        
                 
       !.......................................................................
       !     Determine the equally spaced psi array.
       !     Note: psilim, psimag and epsi will have their sign changed below.
       !     psiar is ordered from edge to magnetic axis.
       !.......................................................................
       
       nfp=nnv 
       ! Up to now, psimag < psilim;  epsi has min. at m.axis
       delpsi=1.e+8*(psilim-psimag)/(nnv-1)  ! positive
       do 10 ix=1,nnv
          !%OS  bug?!
          !%OS  psiar(ix)=1.e+8*psilim+(ix-1)*delpsi
          psiar(ix)=-1.e+8*psilim+(ix-1)*delpsi ! reversed, and mks->cgs
10        continue
          
       endif  !On eqsource
       
       if (eqsource.eq."eqdsk" .or. eqsource.eq."tsc") then
          
          !.......................................................................
          !    Convert to cgs, re-order fpsiar,prar,ppar,qar from edge to mag axis
          !    Reverse sign of epsi,psilim,psimag in the code, so max psi will
          !    occur on axis.
          !.......................................................................
          
          do i=1,nnv
             fpsiar(i)=fpsiar(i)*1.e+6
             ffpar(i)=ffpar(i)/1.e+2
             prar(i)=prar(i)*10.
             ppar(i)=ppar(i)/1.e+7
          enddo
          call dcopy(nnv,fpsiar,1,dummyar,1)
          do i=1,nnv
             fpsiar(i)=dummyar(nnv+1-i)
          enddo
          call dcopy(nnv,ffpar,1,dummyar,1)
          do i=1,nnv
             ffpar(i)=dummyar(nnv+1-i)
          enddo
          call dcopy(nnv,prar,1,dummyar,1)
          do i=1,nnv
             prar(i)=dummyar(nnv+1-i)
          enddo
          call dcopy(nnv,ppar,1,dummyar,1)
          do i=1,nnv
             ppar(i)=dummyar(nnv+1-i)
          enddo
          call dcopy(nnv,qar,1,dummyar,1)
          do i=1,nnv
             qar(i)=dummyar(nnv+1-i)
          enddo
          psilim=-psilim*1.e8  ! reversed !
          psimag=-psimag*1.e8  ! reversed !
          do 5 j=1,nnz
             do 6 i=1,nnr
                epsi(i,j)=-epsi(i,j)*1.e8  ! reversed !
6            continue ! now psimag > psilim;  epsi has min. at edge !
5         continue
          write(*,1000)
          ! 1000      format(//,1x,'WARNING: Sign of epsi,psilim,psimag reversed')
          
          btor=btor*1.e+4
          rbox=rbox*1.e+2
          zbox=zbox*1.e+2
          rboxdst=rboxdst*1.e+2
          radmaj=radmaj*1.e+2
          toteqd=toteqd*3.e9
          zshift=zshift*1.e+2
          
          psi_lim = psilim ! cgs
          psi_mag = psimag ! cgs
          R_axis  = raxis*1.e+2 ! cgs
          Z_axis  = zaxis*1.e+2 ! cgs
          
          do ilim = 1,nlimiter
             rlimiter(ilim)= rlimiter(ilim)*1.d2
             zlimiter(ilim)= zlimiter(ilim)*1.d2
          enddo
          do ilim = 1,ncontr
             rcontr(ilim)= rcontr(ilim)*1.d2
             zcontr(ilim)= zcontr(ilim)*1.d2
          enddo
          
          !..................................................................
          !        Set initial values of rmag,zmag 
          !        (to be refined in eqrhopsi for bicubic splines of epsi).
          !..................................................................
          
          rmag=raxis*1.e+2
          zmag=zaxis*1.e+2
          
          !..................................................................
          !     Re-Create the z,r meshes..
          !     Assumes ymideqd=0.
          !..................................................................
          
          if(ymideqd.ne.zero) stop 'equilib:Consider effect ymideqd.ne.0'
          dzz=zbox/(nnz-1) ! cgs
          drr=rbox/(nnr-1) ! cgs
          er(1)=rboxdst    ! cgs
          ez(1)=-zbox*.5
          do nn=2,nnr
             er(nn)=er(nn-1)+drr
          enddo
          do nn=2,nnz
             ez(nn)=ez(nn-1)+dzz
          enddo
          
          ezmin=ez(1)   ! cgs
          ezmax=ez(nnz)
          ermin=er(1)
          ermax=er(nnr)
          
          !..................................................................
          !     Set up spline array for the f =(R*BTOR) subroutine,etc.
          !..................................................................
          
          i1p(1)=4
          i1p(2)=4
          call coeff1(nnv,psiar,fpsiar,d2fpsiar,i1p,1,workk)
          call coeff1(nnv,psiar,ffpar,d2ffpar,i1p,1,workk)
          call coeff1(nnv,psiar,prar,d2prar,i1p,1,workk)
          call coeff1(nnv,psiar,ppar,d2ppar,i1p,1,workk)
          call coeff1(nnv,psiar,qar,d2qar,i1p,1,workk)
          
          
          ibd(1)=4
          ibd(2)=4
          ibd(3)=4
          ibd(4)=4
          call coeff2(nnr,er,nnz,ez,epsi,epsirr,epsizz,epsirz,nnra,ibd,wkepsi)
          
          !..................................................................
          !     Adjust z-height of diagnostics, if have shifted eqdsk array
          !..................................................................
          
          if (zshift.ne.zero) then
             !           X-ray detector height
             do 500 nn=1,nv
                if (z_sxr(nn).ne.0.) z_sxr(nn)=z_sxr(nn)-zshift
500          continue
             !           Neutral beam pivot point height adjusted in frnfreya.
             !           Must also shift URF ray data, after it is read in.
          endif

       elseif (eqsource.eq."topeol") then
          call eqtopeol
      else
         call eqwrng(10)
      endif
      
      write(*,*)
      write(*,598) rmag,zmag,psimag,psilim,qar(1),qar(nnv)
598   format('ADJUSTED (possible sign change; rescaled) [cgs]: &
           rmag,zmag,psimag,psilim,qar(1),qar(nnv)=',  6(1pe12.4))
      write(*,*)' btor,bsign,fpsiar(nnv)',btor,bsign,fpsiar(nnv)
      write(*,*)'equilib/setup: done'
      !pause


       
   else                      ! index.ne.0 !==========================
      
      !..................................................................
      !     Find  BR, BPHI, BZ for a given Requil,Zequil on equil. grid
      !..................................................................
      der=er(2)-er(1) ! cgs
      dez=ez(2)-ez(1) ! cgs
      
      PSIequil=t2(der,dez,Requil,Zequil,nnr,er,nnz,ez,epsi, &
           epsirr,epsizz,epsirz,nnra,0,0)           ! cgs
      dpsidr=  t2(der,dez,Requil,Zequil,nnr,er,nnz,ez,epsi, &
           epsirr,epsizz,epsirz,nnra,1,0)
      dpsidz=  t2(der,dez,Requil,Zequil,nnr,er,nnz,ez,epsi, &
           epsirr,epsizz,epsirz,nnra,0,1)
      !..................................................................
      !     Determine f(psi)
      !..................................................................
      
      itab(1)=1
      itab(2)=0
      itab(3)=0
      dpsiar=psiar(2)-psiar(1) ! cgs
      call terp1e(nfp,psiar,fpsiar,d2fpsiar,PSIequil, &
           1,tab,itab,dpsiar)
      ! YuP Note: The sign of fpsiar was adjusted
      ! during setup by using bsign (see line ~237 above) 
      fpsieq=tab(1)
      
      !..................................................................
      !     Determine components of B.
      !..................................................................
         
      BReq=  -cursign*dpsidz/Requil     ! cgs
      BPHIeq= bsign*fpsieq/Requil ! change to original sign of Bphi
      BZeq=   cursign*dpsidr/Requil
      
      !         if(abs(Zequil+0.8).le.1.e-2) then
      !          write(*,'(a,4e12.4)')'equilib',Requil,Zequil,PSIequil,BZeq
      !         endif
         
   endif                     ! endif on index !======================

   return
 end subroutine

end module equilib_mod
