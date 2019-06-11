module micxiniz_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use diagwrng_mod, only : diagwrng
  use micgetr_mod, only : micgetr
  use tdxin13d_mod, only : tdxin13d
  use psif_mod, only : psif
  use psif_mod, only : psifp
  use zcunix_mod, only : coeff1
  use zcunix_mod, only : terp1

  !---END USE

!
!

contains

      subroutine micxiniz
      use param_mod
      use cqlcomm_mod
      use cqlconf_mod, only : setup0
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!..................................................................
!     This routine determines certain constants and the mesh along
!     the magnetic field line z and sz, as well as related quantities.
!     This was part of micxinit in CQL3D.
!     These meshes are tailored to suit requirements specified by
!     the user in the namelist input.
!..................................................................

      dimension zd2bpol(lfielda),zd2solr(lfielda),zd2solz(lfielda), &
        znormsh(0:lza)

!%OS  compute psi_prime as df/ds
!%OS  fpsismi(l)=psis(l)*dls(2,2,1,l)+psis(l+1)*(1-dls(2,2,1,l))
      fpsismi(l)=(0.5-0.5*(1/(l+1))+0.5*(1/(setup0%ls+1-l)))*psis(l)+ &
        (0.5+0.5*(1/(l+1))-0.5*(1/(setup0%ls+1-l)))*psis(l+1)

!.......................................................................
!l    1. Determine the field line orbit mesh z(1:lz,lr_).
!     tfacz is a mesh squeeze (nonuniformity) factor:
!     tfacz=1. means uniform mesh; tfacz < 1. means geometric mesh with
!     points closer near zero. tfacz > 1. means points closer near zmax
!     For eqsym.eq."none", tfacz appears best, since z(lz,)=z(1,).
!......................................................................

      if (eqsym.ne."none") then
         z(1,lr_)=0.
         z(lz,lr_)=zmax(lr_)
         htfg=tfacz*(zmax(lr_)/(lz-1))
         z(2,lr_)=htfg
         call micgetr(lz,zmax(lr_),htfg,ram,ksingul)
         if (ksingul .eq. 1) call diagwrng(8)
         do 100 l=3,lz-1
            z(l,lr_)=z(l-1,lr_)+ram*(z(l-1,lr_)-z(l-2,lr_))
 100     continue
         do l=1,lz-1
            zmid(l,lr_)=0.5*(z(l,lr_)+z(l+1,lr_))
         enddo
         zmid(lz,lr_)=zmax(lr_)
         z_bmax(lr_)=zmax(lr_)
         lz_bmax(lr_)=lz
         bpsi_z_bmax(lr_)=psif(zmax(lr_))
!         do l=1,lz
!         write(*,*) l,z(l,lr_),zmid(l,lr_)
!         enddo
!         write(*,*)'micxiniz: l,z,zmid for lr_=',lr_
!         pause
      endif

!......................................................................
!     If eqsym.eq."none", use equispaced mesh on either side of
!     the inboard max-|B| point, with step size slightly adjusted
!     so the max-|B| point is a mesh point.

      if (eqsym.eq."none") then
         z_bmax(lr_)=es_bmax(lr_)
         !-YuP Was: lz_bmax(lr_)=nint(z_bmax(lr_)/zmax(lr_)*lz) !old version
         lz_bmax(lr_)= 1 + nint(z_bmax(lr_)/zmax(lr_)*(lz-1)) !-YuP Try (lz+1)/2
         dlz1=z_bmax(lr_)/(lz_bmax(lr_)-1)
         !-YuP Was: dlz2=(zmax(lr_)-z_bmax(lr_))/(lz-lz_bmax(lr_)-1) !old version
         dlz2=(zmax(lr_)-z_bmax(lr_))/(lz-lz_bmax(lr_)) !-YuP
         do l=1,lz_bmax(lr_)
            z(l,lr_)=(l-1)*dlz1                 !range: [0., z_bmax]
         enddo
         do l=lz_bmax(lr_)+1,lz
            z(l,lr_)=z_bmax(lr_)+(l-lz_bmax(lr_))*dlz2 ![z_bmax+dlz2, zmax]
         enddo
         do l=1,lz-1
            zmid(l,lr_)=0.5*(z(l,lr_)+z(l+1,lr_))
         enddo
         zmid(lz,lr_)=zmax(lr_) !-YuP: this is an extra ending point
         !-YuP Note: zmid is defined at lz-1 points;
         !-YuP The mesh z is defined at lz points,
         !-YuP  but last point = 1st point in R,Z coords.
!         do l=1,lz
!         write(*,*) l,z(l,lr_),zmid(l,lr_)
!         enddo
!         write(*,*)'micxiniz: l,z,zmid for lr_=',lr_
!         pause
         bpsi_z_bmax(lr_)=bpsi_max(lr_)  !New variable name for z-mesh.
      endif


!..................................................................
!l    2. Determine the bpsi array (ratio of B(z)/B(0)).
!     coeff1, terp1 are spline routines frequently utilized by CQL.
!..................................................................

      if (psimodel.eq."spline") then
        if (eqmod.eq."disabled") then
          lfield=lfielda
          lorbit(lr_)=lfield
          delz=zmax(lr_)/(lfield-1)
          do 24 lu=1,lorbit(lr_)
            es(lu,lr_)=(lu-1)*delz
            thtpol(lu,lr_)=pi*es(lu,lr_)*zmaxi(lr_)
            bpsi(lu,lr_)=(1.+eps(lr_))/(1.+eps(lr_)*cos(thtpol(lu,lr_)))
            eqbpol(lu,lr_)=bthr0(lr_)*bpsi(lu,lr_)
            solr(lu,lr_)=radmaj*(1.+eps(lr_)*cos(thtpol(lu,lr_)))
 24       continue
        endif
        if (eqsym.eq."none") then
           i1p(1)=3
           i1p(2)=3
        else
           i1p(1)=2
           i1p(2)=2
           d2bpsi(lorbit(lr_),lr_)=0.
           d2bpsi(1,lr_)=0.
        endif

!-YuP Print this:
!         step2=55./nstps
!         pt=0.
!         do il=1,nstps
!            write(*,*)'micxiniz before coeff1 pt,psif=',pt,psif(pt)
!            pt=pt+step2
!         enddo

!-YuP: Here psif(pt) is ok; monotonically increasing with pt.

!-YuP        call coeff1(lorbit(lr_),es(1,lr_),bpsi(1,lr_),d2bpsi(1,lr_)
!-YuP     1    ,i1p,1,work) !-> to get d2bpsi (2nd derivatives)

!-YuP: Here psif(pt) is bad - non-monotonic (after coeff1 interpolation)

!-YuP From print out of d2bpsi: The 2nd derivative is quite bad.
!-YuP It results in an error with spline interpolation
!-YuP in psif(pt), near B=Bmin point.
!-YuP However, if these two lines calling coeff1 are commented
!-YuP (which makes all d2bpsi remaining zeroed),
!-YuP the subsequent call of terp1 in psif(pt) is
!-YuP equivalent to a linear interpolation.
!-YuP In this case there is no problem around B=Bmin point.


        i1p(1)=4
        i1p(2)=4
        call coeff1(lorbit(lr_),es(1:lorbit(lr_),lr_),thtpol(1:lorbit(lr_),lr_),d2thtpol(1:lorbit(lr_),lr_),i1p,1,work)
        itab(1)=1
        itab(2)=0
        itab(3)=0



!..................................................................
!     Determine an evenly spaced psi array.
!..................................................................
!%OS
!%OS
!%OS  dels=es(lorbit(lr_),lr_)/float(incza-1)
!%OS  psiesfi(1,lr_)=1.
!%OS  esfi(1,lr_)=0.
!%OS  do 29 iv=2,incza-1
!%OS  esfi(iv,lr_)=(iv-1)*dels
!%OS  call terp1(lorbit(lr_),es(1,lr_),bpsi(1,lr_),
!%OS  1     d2bpsi(1,lr_),esfi(iv,lr_),1,tab,itab)
!%OS  psiesfi(iv,lr_)=tab(1)
!%OS29continue
!%OS  psiesfi(incza,lr_)=bpsi(lorbit(lr_),lr_)
!%OS  esfi(incza,lr_)=es(lorbit(lr_),lr_)
!%OS
!%OSc..................................................................
!%OSc Determine an evenly spaced (in argument) s array (arg is psi)
!%OSc..................................................................
!%OS
!%OS  call coeff1(lorbit(lr_),bpsi(1,lr_),es(1,lr_),d2es(1,lr_),
!%OS  1     i1p,1,work)
!%OS  deltapsi(lr_)=(bpsi(lorbit(lr_),lr_)-1.)/(incza-1)
!%OS  psifi(1,lr_)=1.
!%OS  espsifi(1,lr_)=0.
!%OS  do 28 iv=2,incza-1
!%OS  psifi(iv,lr_)=psifi(iv-1,lr_)+deltapsi(lr_)
!%OS  call terp1(lorbit(lr_),bpsi(1,lr_),es(1,lr_),d2es(1,lr_),
!%OS  1      psifi(iv,lr_),1,tab,itab)
!%OS  espsifi(iv,lr_)=tab(1)
!%OS28continue
!%OS  espsifi(incza,lr_)=zmax(lr_)
!%OS  espsifi(inczpa,lr_)=zmax(lr_)
!%OS  psifi(inczpa,lr_)=bpsi(lorbit(lr_),lr_)
!%OS  psifi(incza,lr_)=bpsi(lorbit(lr_),lr_)
        pol(1,lr_)=0.
        if (eqsym.eq."none") then
           pol(lz,lr_)=2*pi
        else
           pol(lz,lr_)=pi
        endif
        do 22 l=2,lz-1
          call terp1(lorbit(lr_),es(1:lorbit(lr_),lr_),thtpol(1:lorbit(lr_),lr_),d2thtpol(1:lorbit(lr_),lr_) &
            ,z(l,lr_),1,tab,itab)
          pol(l,lr_)=tab(1)
 22     continue

!.......................................................................
!     solrz(1:lz,lr_),solzz(1:lz,lr_)
!.......................................................................

        if (eqsym.eq."none") then
           i1p(1)=3
           i1p(2)=3
        else
           i1p(1)=2
           i1p(2)=2
           d2solrz(1,lr_)=0.0
           d2solrz(lorbit(lr_),lr_)=0.0
        endif
        call coeff1(lorbit(lr_),es(1:lorbit(lr_),lr_),solr(1:lorbit(lr_),lr_),d2solrz(1:lorbit(lr_),lr_), &
          i1p,1,work)
        do 221 l=1,lz
           call terp1(lorbit(lr_),es(1:lorbit(lr_),lr_),solr(1:lorbit(lr_),lr_),d2solrz(1:lorbit(lr_),lr_), &
                z(l,lr_),1,tab,itab)
           solrz(l,lr_)=tab(1)
 221    continue

        if (eqsym.eq."none") then
           i1p(1)=3
           i1p(2)=3
        else
           i1p(1)=2
           i1p(2)=2
           d2solzz(1,lr_)=0.0
           d2solzz(lorbit(lr_),lr_)=0.0
        endif
        call coeff1(lorbit(lr_),es(1:lorbit(lr_),lr_),solz(1:lorbit(lr_),lr_),d2solzz(1:lorbit(lr_),lr_), &
             i1p,1,work)
        do 223 l=1,lz
           call terp1(lorbit(lr_),es(1:lorbit(lr_),lr_),solz(1:lorbit(lr_),lr_),d2solzz(1:lorbit(lr_),lr_), &
                z(l,lr_),1,tab,itab)
           solzz(l,lr_)=tab(1)
 223    continue



      else if(psimodel.eq."axitorus") then

        do 25 l=1,lz
          pol(l,lr_)=pi*z(l,lr_)*zmaxi(lr_)
 25     continue

!BH000416        if (vlfmod.eq."enabled") then
          do 26 l=1,lz
            solrz(l,lr_)=radmaj*(1.+eps(lr_)*cos(pol(l,lr_)))
           !solzz(l,lr_)=radmaj*(1.+eps(lr_)*sin(pol(l,lr_))) ! YuP[07-2017] BUG
            solzz(l,lr_)=radmaj*( zmag/rmag +eps(lr_)*sin(pol(l,lr_)) ) !correct
            !Note: rmag=radmaj in this model, and zmag=0.
!            write(*,'(a,2i4,3f13.2)')
!     +         'micxiniz: lr_,l, solrz(),solzz(),pol(l,lr_)',
!     +          lr_,l, solrz(l,lr_), solzz(l,lr_), pol(l,lr_)
 26       continue
!BH000416        endif

      endif  ! On psimodel

      bbpsi(1,lr_)=1.+em14
      do 220 l=2,lz
        bbpsi(l,lr_)=psif(z(l,lr_))
 220  continue
      if (eqsym.eq."none") bbpsi(lz,lr_)=bbpsi(1,lr_)

!.......................................................................
!l    3. Construct mesh and profiles in CQLP mode on given flux surface
!     Note: If setup0%cqlpmod=enabled, lz=setup0%lsmax, full mesh along B
!     setup0%ls gives mesh on which equations are solved
!     (if periodic: sz(1)=sz(setup0%ls))
!     psis=B(s)/B(midplane), psisp=d(psis)/ds
!.......................................................................

!.......................................................................
!l    3.1 s, sz, psis, psisp and related meshes
!     Note: only one flux surface in CQLP so far (setup0%lrz=1)
!.......................................................................

      if (setup0%cqlpmod.ne."enabled" .or. lr_.ne.setup0%lrindx(1)) go to 999  !To END

!     Following, not checked for eqsym.eq."none". Possibly OK. BH090923.

!     Assumes midplane at l=1. If not, bmidplne, z mesh, bpsi etc should be
!     redefine
      do 310 l=1,setup0%ls
        sz(l)=z(setup0%lsindx(l),lr_)
        psis(l)=bbpsi(setup0%lsindx(l),lr_)/bbpsi(lmdpln(indxlr(lr_)),lr_)
        psisp(l)=psifp(z(setup0%lsindx(l),lr_))/bbpsi(lmdpln(indxlr(lr_)),lr_)
 310  continue

!     psipols, solrs, solzs
      i1p(1)=2
      i1p(2)=2
      zd2bpol(1)=0.0
      zd2bpol(lorbit(lr_))=0.0
      call coeff1(lorbit(lr_),es(1:lorbit(lr_),lr_),eqbpol(1:lorbit(lr_),lr_),zd2bpol(1:lorbit(lr_)), &
        i1p,1,work)
      i1p(1)=2
      i1p(2)=2
      zd2solr(1)=0.0
      zd2solr(lorbit(lr_))=0.0
      call coeff1(lorbit(lr_),es(1:lorbit(lr_),lr_),solr(1:lorbit(lr_),lr_),zd2solr(1:lorbit(lr_)), &
        i1p,1,work)
      i1p(1)=4
      i1p(2)=4
      call coeff1(lorbit(lr_),es(1:lorbit(lr_),lr_),solz(1:lorbit(lr_),lr_),zd2solz(1:lorbit(lr_)), &
        i1p,1,work)

      itab(1)=1
      itab(2)=0
      itab(3)=0

      do 311 l=1,setup0%ls
        call terp1(lorbit(lr_),es(1:lorbit(lr_),lr_),eqbpol(1:lorbit(lr_),lr_),zd2bpol(1:lorbit(lr_)) &
          ,sz(l),1,tab,itab)
        psipols(l)=tab(1)/bmidplne(lr_)
        call terp1(lorbit(lr_),es(1:lorbit(lr_),lr_),solr(1:lorbit(lr_),lr_),zd2solr(1:lorbit(lr_)) &
          ,sz(l),1,tab,itab)
        solrs(l)=tab(1)
        call terp1(lorbit(lr_),es(1:lorbit(lr_),lr_),solz(1:lorbit(lr_),lr_),zd2solz(1:lorbit(lr_)) &
          ,sz(l),1,tab,itab)
        solzs(l)=tab(1)
 311  continue
      do 312 l=2,setup0%ls-1
        dsz(l)=0.5*(sz(l+1)-sz(l-1))
        dszp5(l)=sz(l+1)-sz(l)
        eszp5(l)=1./dszp5(l)
        dszm5(l)=sz(l)-sz(l-1)
        eszm5(l)=1./dszm5(l)
 312  continue
      dsz(1)=(sz(2)-sz(1))*0.5
      dsz(setup0%ls)=(sz(setup0%ls)-sz(setup0%ls-1))*0.5
      dszm5(1)=0.0
      eszm5(1)=0.0
      dszp5(1)=sz(2)-sz(1)
      eszp5(1)=1./dszp5(1)
      dszm5(setup0%ls)=sz(setup0%ls)-sz(setup0%ls-1)
      eszm5(setup0%ls)=1./dszm5(setup0%ls)
      dszp5(setup0%ls)=0.0
      eszp5(setup0%ls)=0.0

!.......................................................................
!l    3.2 Construct parallel profile n(s) and T(s)
!.......................................................................

!     tdxin13d assumes a normalized mesh
      znormsh(0)=0.0
      do 320 l=1,lz
        znormsh(l)=z(l,lr_)/z(lz,lr_)
 320  continue
!     only parabola option for n(s), T(s) so far
      do 321 ik=1,ntotal
        call tdxin13d(denpar,znormsh,setup0%lsmax,ntotala,ik,npwr(0),mpwr(0))
        call tdxin13d(temppar,znormsh,setup0%lsmax,ntotala,ik,npwr(ik), &
          mpwr(ik))
 321  continue

!.......................................................................
!l    4. Define end points accordingly whether mesh is periodic or not
!.......................................................................

      if (transp .ne. "enabled") go to 999

!     Note: if transp=enabled then setup0%ls=setup0%lsmax
      if (sbdry.ne."periodic" .and. numclas.ne.1) then
        do 400 ik=1,ntotal
          do 401 il=0,setup0%lsmax+1,setup0%lsmax+1
            denpar(ik,il)=0.0
            temppar(ik,il)=0.0
 401      continue
 400    continue
        do 410 il=0,setup0%lsmax+1,setup0%lsmax+1
          sz(il)=0.0
          psisp(il)=0.0
          psipols(il)=0.0
          solrs(il)=0.0
          solzs(il)=0.0
          dsz(il)=0.0
          dszm5(il)=0.0
          dszp5(il)=0.0
          eszm5(il)=0.0
          eszp5(il)=0.0
 410    continue
      else
!     0 <=> setup0%lsmax and setup0%lsmax+1 <=> 1
        do 420 ik=1,ntotal
          do 421 il=0,setup0%lsmax+1,setup0%lsmax+1
            iequiv=il/(setup0%lsmax+1)+setup0%lsmax*((setup0%lsmax+1-il)/(setup0%lsmax+1))
            denpar(ik,il)=denpar(ik,iequiv)
            temppar(ik,il)=temppar(ik,iequiv)
 421      continue
 420    continue
        do 430 il=0,setup0%lsmax+1,setup0%lsmax+1
          iequiv=il/(setup0%lsmax+1)+setup0%lsmax*((setup0%lsmax+1-il)/(setup0%lsmax+1))
          sz(il)=sz(iequiv)
          psisp(il)=psisp(iequiv)
          psipols(il)=psipols(iequiv)
          solrs(il)=solrs(iequiv)
          solzs(il)=solzs(iequiv)
          dsz(il)=dsz(iequiv)
          dszm5(il)=dszm5(iequiv)
          dszp5(il)=dszp5(iequiv)
          eszm5(il)=eszm5(iequiv)
          eszp5(il)=eszp5(iequiv)
 430    continue
      endif

!.......................................................................
!l    5. Define dpsis/ds in the same way as df/ds in transport equation
!.......................................................................

!%OS
!%OS  if (noffso(1,2) .eq. 999) go to 999
!%OS
!%OS  do 500 l=1,setup0%ls
!%OS  psisp(l)=(fpsismi(l)-fpsismi(l-1))/dsz(l)
!%OS  500  continue
!%OS  psisp(0)=0.0
!%OS  psisp(setup0%ls+1)=0.0
!%OS  if (sbdry .eq. "periodic") then
!%OS  psisp(0)=psisp(setup0%ls)
!%OS  psisp(setup0%ls+1)=psisp(1)
!%OS  endif
!%OS

!.......................................................................

 999  continue

      return
      end subroutine micxiniz


end module micxiniz_mod
