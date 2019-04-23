module tdoutput_mod

  !---BEGIN USE

  use bcast_mod, only : bcast
  use cfpgamma_mod, only : cfpgamma
  use restcon_mod, only : restcon
  use resthks_mod, only : resthks
  use tdnflxs_mod, only : tdnflxs
  use tdtrflx_mod, only : tdtrflx

  !---END USE

!
!

contains

      subroutine tdoutput(kopt)
      use param_mod
      use comm_mod
      implicit integer (i-n), real*8 (a-h,o-z)

!..................................................................
!     diagnostics on paper
!     kopt=1 ==> printout before first time step, at end of ainitial
!                (if lrzmax.eq.1) or end of tdinitl
!     kopt=2 ==> printout after specified number of time steps,
!                including at least the last time step.
!..................................................................

!MPIINSERT_INCLUDE

      dimension ztr(lrza)
      dimension xj(lrza),xp(lrza),xe(lrza),xc(lrza)
      dimension rban_vth(lrzmax),dn_scale(lrzmax),dt_scale(lrzmax) !print-out
      character(len=8) :: ztext


!MPIINSERT_IF_RANK_NE_0_RETURN
 ! PRINT on mpirank.eq.0 only


!.......................................................................
!l    0. initialize "radial" mesh for print out
!.......................................................................
      if (pltvs .ne. "psi") then
        ztext=" rovera "
        do 10 l=1,lrzmax
          ztr(l)=rovera(l)
 10     continue
      else
        ztext=pltvs
        do 11 l=1,lrzmax
          ztr(l)=(equilpsi(0)-equilpsi(l))/equilpsi(0)
 11     continue
      endif

      call bcast(xj,zero,lrza)
      call bcast(xp,zero,lrza)
      call bcast(xe,zero,lrza)
      call bcast(xc,zero,lrza)

!.......................................................................
!     kopt=1 ==> printout before first time step, at end of ainitial
!                (if lrzmax.eq.1) or end of tdinitl
!     kopt=2 ==> printout after specified number of time steps,
!                including at least the last time step.

      go to (100,200) kopt

!.......................................................................
!l    1.0   Label the version and prepare versus which variable
!     the output values are listed
!
!     OS stands for Olivier Sauter, presently (Nov 23/93) at CRPP
!     Lausanne, Switzerland.
!     bh is Bob Harvey, CompX, Del Mar, CA.
!     Major change of code implies change in letter (A, B, ...).
!     Minor changes implies change in final date.
!     Final synthesis of OS and BH version is first B version: B_940425
!     bh_940222 is a synthesis of OS_A_940128 and the bh changes below.
!     The version OS_F_931123 used version OS_E_930217, and has added
!     into it capacity for cyclotron damping on the ions (in urf routines)
!     (previously only electrons were treated), and for damping by
!     multiple cyclotron harmonics. (Bob Harvey: (858)509-2131, 11/23/93.
!     The version OS_E_930217 used version OS_D_930209 and added the option
!     of solving CQL on a few flux surfaces, while having many radial
!     surfaces as main mesh (lrz<lrzmax, lrzdiff="enabled")
!
!.......................................................................
 100  continue

      WRITE(t_,1000) version
 1000 format(15x,"CQL3D/CQLP VERSION: ",a)
      WRITE(6,'(//16x, &
                "=================================================")')
      WRITE(6,*) t_(1:length_char(t_))
      WRITE(6,'(16x, &
                "=================================================",/)')

!.......................................................................
!l    1.1 some equilibrium quantities
!.......................................................................

      WRITE(6,9110) radmaj, rhomax, &
        r0geomp, rgeomp, eps0, &
        zgeomp, rgeom1, rgeom2, zgeomp/rgeomp, &
        rmag, rmag-r0geomp,zshift, &
        psimag,psilim, &
        btor, toteqd/3.e9 &
        ,volmid(lrzmax), areamid(lrzmax)
      if (cqlpmod .eq. "enabled") WRITE(6,9114) z(1,lrindx(1)), &
        z(lsmax,lrindx(1))

      WRITE(6,9111)
      do l=1,lrzmax
      if( (psilim-psimag).ne.0.d0 ) then !YuP[2018-01-03] added
        ! Sometimes psilim=psimag (in lrz=1 runs)
        sqrt_psinorm=sqrt((equilpsi(l)-psimag)/(psilim-psimag))
      else
        sqrt_psinorm=0.d0
      endif
      WRITE(6,91111) l,rovera(l),eps(l),btor0(l),bthr0(l), &
        bmod0(l),btor0(l)/bmod0(l),qsafety(l),psimx(l), &
        equilpsi(l),  sqrt_psinorm
      enddo


      WRITE(6,9109) (l,rmcon(l),zmcon(l),rpconz(l),zpcon(l),l=1,lrzmax)
!
      WRITE(6,9112)
!     The last number in the following write is the ratio
!     of curr/(<j.B>/<B>), where curr is the toroidal-area
!     averaged parallel current plotted out from cql3d.
      do 110 ll=1,lrzmax
        zr0=rpcon(ll)
        WRITE(6,9113) ll,rgeom(ll),rovera(ll)*rhomax &
          ,zr0*onovrp(1,ll),psiavg(1,ll) &
          ,zr0**2*onovrp(2,ll),psiovr(ll)*zr0 &
          ,psiavg(2,ll),darea(ll),dvol(ll) &
          ,psiovr(ll)/onovrp(1,ll)*psiavg(1,ll)/psiavg(2,ll)
 110  continue
!
!l    1.1.2 parallel quantities
!     (assumes lrz=1)
      if (cqlpmod .eq. "enabled") then
        ilr=lrindx(1)
        WRITE(6,9115)
        do 112 ll=1,lrors
          zbfiob=fpsi(ilr)/solrs(ll)/psis(ll)/bmidplne(ilr)
          WRITE(6,9116) lsindx(ll),sz(ll),solrs(ll),solzs(ll), &
            zbfiob,psis(ll),psisp(ll),psipols(ll)
 112    continue
      endif
!
 9110 format(/," equilibrium parameters",/,1x,22("="),// &
        ,7x,"Nominal major radius=",1pe14.4,10x,"rhomax      =",e14.4/ &
        ,7x,"Geometrical center  =", e14.4,10x,"minor radius=",e14.4 &
        ,10x,"Inv Aspect ratio=",e10.4,/,7x, &
        "Maximum height      =",e14.4,10x,"Inner Maj R =",e14.4 &
        ,10x,"Outer Maj R =",e14.4,10x,"elongation    =",e12.4/ &
        ,7x,"Magnetic axis       =",  e14.4,10x,"radial shift=",e14.4 &
        ,10x,"vertical zshift=",e11.4/ &
        ,7x,"Pol flux psimag     =",  e14.4,10x,"Pol flux psilim="e11.4, &
       /,7x,"Nominal tor mag fld =",  e14.4,10x,"Equ.curr.[A]=",e14.4,/ &
        ,7x,"Total volume        =",  e14.4,10x,"Total area  =",e14.4)

 9111 format(/"  l",5x,"rovera",2x,"inv aspct ratio",3x,"btor(0)",6x, &
        "bpol(0)",8x,"b(0)",6x,"b_phi/b(0)",7x,"q",8x,"Bmax/Bmin" &
        ,4x,"equilpsi",  4x, "sqrt((psi-psimag)/(psilim-psimag))")

91111  format(i3,1p10e13.4)
 9109 format(/"  l",5x,"Min R",8x,"Z@Min R",6x,"Max R",8x, &
        "Z@Max R"/,(i3,1p4e13.4))
 9112 format(/,"  l",3x,"[Rp-Rm](l)",5x,"rho(l)",6x,"<R(0)/R>",1x, &
      "<B/B(0)=bpsi>",1x,"<(R(0)/R)^2>",2x,"R(0)<bpsi/R>",3x,"<bpsi**2>" &
        ,5x,"darea",9x,"dvol",1x,"(<bpsi/R>/<1/R>)/(<bpsi**2>/<bpsi>)")
 9113 format(i3,1p10e13.4)
 9114 format(7x,"smin",16x,"=",1pe14.4,10x,"smax",8x,"=",e14.4)
 9115 format(/," Equilibrium quantities along the magnetic field:",/ &
        " ===============================================",//, &
        " ls_",6x,"s",10x,"R(s)",10x,"Z(s)",7x,"B_phi/B",6x, &
        "B(s)/B(0)",2x,"1/B(0)*dB/ds",1x,"B_pol(s)/B(0)")
 9116 format(i3,1p7e13.4)



!.......................................................................
!l    1.2.0 species numbers
!.......................................................................

      WRITE(6,8110)
 8110 format(//,'Species numbers, from ainspec.f:',/, &
      'kelecg = the index of the general species electrons ',/, &
      'kelecm = the background Maxwlln species electrons (0,if none)',/, &
      'kionn =  1st ionic species (general or background species)',/, &
      'kelec = kelecg (if .ne.0),  kelecm (else)',/, &
      'niong = number of general species ions',/, &
      'kiong(1:niong) = general ion specie (otherwise 0 for 1:ngen)',/, &
      'nionm = number of background Maxwellian species ions',/, &
      'kionm(1:nionm) = Maxwellian ion specie indices (otherwise 0)')

      WRITE(6,8111) kelecg,kelecm,kionn,kelec,niong
      if (niong.ne.0) then
         WRITE(6,8112) (kiong(i),i=1,niong)
      else
         WRITE(6,*) 'kiong(1) =  0'
      endif
      WRITE(6,8113) nionm
      if (nionm.ne.0) then
         WRITE(6,8114) (kionm(i),i=1,nionm)
      else
         WRITE(6,*) 'kionm(1) =  0'
      endif

      WRITE(6,*)
! 9118 format(30i3)
 8111 format(' kelecg =',i3,/,' kelecm =',i3,/,' kionn =',i3,/, &
           ' kelec =',i3,/,' niong =',i3)
 8112 format(' kiong(1:niong) =',10i3)

 8113 format(' nionm =',i3)
 8114 format(' kionm(1:nionm) =',10i3)


!.......................................................................
!l    1.2.1  Coulomb logrithm for general species at before initial step
!.......................................................................

      do igen=1,ngen
         WRITE(6,*)
         WRITE(6,8115) igen, n
         do ll=1,lrz
            call tdnflxs(ll)
            call cfpgamma
            WRITE(6,8116) lr_,(gama(igen,ispec),ispec=1,ntotal)
         enddo
      enddo
      WRITE(6,*)
 8115 format('Coulomb logrithm for general species number igen =',i2, &
           '   time step =', i4)
 8116 format(' lr_ =',i3, &
              '  ln_lambda(igen,ispec=1,ntotal) =',10(1pe10.3))

!.......................................................................
!l    1.2 density,temperature,magnetic field,etc.
!.......................................................................

!dir$ nextscalar
      do 123 jk=1,ntotal
         do ll=1,lrzmax
           qb_mc=bnumb(jk)*charge*bthr(ll)/(fmass(jk)*clight)
           ! Banana width (for v=vthermal, at t-p bndry):
           ! BH171230: For lrzdiff=enabled, lrz<lrzmax, the default
           ! values of itl_(ll>lrz) can cause out-of-bounds coss.
           if (ll.le.lrz) then
           rban_vth(ll)= abs(vth(jk,ll)*coss(itl_(ll),ll)/qb_mc) ! cm
           else
              rban_vth(ll)= zero
           endif
           llp=ll+1
           llm=ll-1
           llp=min(llp,lrzmax) ! not to exceeed lrzmax
           llm=max(llm,1)
           ll0=1 ! reference density or temp.: at plasma center.
           ! Density Gradient scale; defined as n(center)/(dn/dr) :
           dndr=(reden(jk,llp)-reden(jk,llm))/(rpcon(llp)-rpcon(llm))
           ! Temperature Gradient scale; defined as T(center)/(dT/dr) :
           dtdr=(temp(jk,llp)-temp(jk,llm))/(rpcon(llp)-rpcon(llm))
           if(dndr.ne.0.d0)then
              dn_scale(ll)= abs(reden(jk,ll0)/dndr)  ! cm
           else ! dn/dr=0
              dn_scale(ll)= 1000*rmag ! just any large number [~Inf]
           endif
           if(dtdr.ne.0.d0)then
              dt_scale(ll)= abs(temp(jk,ll0)/dtdr)  ! cm
           else ! dT/dr=0
              dt_scale(ll)= 1000*rmag ! just any large number [~Inf]
           endif
         enddo ! ll
!        The NPA signal is derived by neutralization from ion general
!        species nnspec.
         if (ipronn.eq."disabled".or.jk.ne.nnspec) then
            if(kspeci(2,jk).eq.'general' .and. &
                        iprovphi.ne.'disabled') then
              WRITE(6,9126) jk,kspeci(1,jk),kspeci(2,jk),bnumb(jk),ztext
              WRITE(6,9127) (il,ztr(il),reden(jk,il),temp(jk,il), &
                 vth(jk,il),energy(jk,il),vphipl(il),il=1,lrzmax)
            else
              WRITE(6,9120) jk,kspeci(1,jk),kspeci(2,jk),bnumb(jk),ztext
              WRITE(6,9121) (il,ztr(il),reden(jk,il),temp(jk,il), &
                 vth(jk,il), energy(jk,il), &
                 rban_vth(il), dn_scale(il), dt_scale(il), il=1,lrzmax)
            endif
         elseif (nnspec.eq.jk) then
            if(kspeci(2,jk).eq.'general' .and. &
                        iprovphi.ne.'disabled') then
              WRITE(6,9128) jk,kspeci(1,jk),kspeci(2,jk),bnumb(jk),ztext
              WRITE(6,9129) (il,ztr(il),reden(jk,il),temp(jk,il), &
                 vth(jk,il),energy(jk,il),enn(il,1),vphipl(il), &
                 il=1,lrzmax)
            else
              WRITE(6,9122) jk,kspeci(1,jk),kspeci(2,jk),bnumb(jk),ztext
              WRITE(6,9123) (il,ztr(il),reden(jk,il),temp(jk,il), &
                 vth(jk,il),energy(jk,il),enn(il,1),il=1,lrzmax)
            endif
         endif
        if (cqlpmod .eq. "enabled") &
          WRITE(6,9125) lrindx(1),(il,z(il,lrindx(1)),denpar(jk,il) &
          ,temppar(jk,il),vthpar(jk,il),enrgypa(jk,il),il=1,lsmax)
 123  continue
      !pause
!
 9120 format(/" species no. ",i3,2x,a,a8,"    charge number: ",f6.2 &
        ,/,1x,15("="),//,"  l",4x,a8,5x,"density",4x, &
        "temperature",6x,"vth",9x,"energy", &
        4x,"rban_vth(cm)", 1x,"dn_scale(cm)", 1x,"dT_scale(cm)" )
 9121 format(i3,1p8e13.5)
 9122 format(/" species no. ",i3,2x,a,a8,"    charge number: ",f6.2 &
        ,/,1x,15("="),//,"  l",4x,a8,5x,"density",4x, &
        "temperature",6x,"vth",9x,"energy",4x,"neutral den")
 9123 format(i3,1p6e13.5)
 9125 format(/," along magnetic field line, at lrindx=",i3,":",/, &
        9x,"s",/,(i3,1p5e13.5))
 9126 format(/" species no. ",i3,2x,a,a8,"    charge number: ",f6.2 &
        ,/,1x,15("="),//,"  l",4x,a8,5x,"density",4x, &
        "temperature",6x,"vth",9x,"energy",3x,"tor vel(cm/s)")
 9127 format(i3,1p6e13.5)
 9128 format(/" species no. ",i3,2x,a,a8,"    charge number: ",f6.2 &
        ,/,1x,15("="),//,"  l",4x,a8,5x,"density",4x, &
        "temperature",6x,"vth",9x,"energy",4x,"neutral den",2x, &
        "tor vel(cm/s)")
 9129 format(i3,1p7e13.5)

!.......................................................................
!l    1.3 tauee, vth, etc.
!.......................................................................


      if (cqlpmod .ne. "enabled") then
         WRITE(6,9130) (il,tauee(il),zmax(il),zmax(il)/vth(1,il), &
              zmax(il)/vth(2,il),taueeh(il),starnue(il), &
          starnue(il)*eps(il)**1.5,zeff(il),il=1,lrzmax)
      else
        WRITE(6,9131) (il,tauee(il),zmax(il),zmax(il)/vth(1,il), &
              zmax(il)/vth(2,il),taueeh(il),starnue(il), &
          starnue(il)*eps(lrindx(1))**1.5,zeff(il),il=1,lsmax)
      endif

!
 9130 format(//,"  lr","   tauee(l)  ","   zmax(l)   "," zmax/vth(1) " &
        ," zmax/vth(2) ",4x,"taueeh",7x,"nuestar",2x,"nuest*eps**3/2", &
        4x,"zeff"/,(i3,1p8e13.4))
 9131 format(//," ls","   tauee(l)  ","   zmax(l)   "," zmax/vth(1) " &
        ," zmax/vth(2) ",4x,"taueeh",7x,"nuestar",2x,"nuest*eps**3/2", &
        4x,"zeff",/,(i3,1p8e13.4))

!.......................................................................
!l    1.3.1 Neoclassical related output.  See tdrmshst.f.
!.......................................................................

      if (cqlpmod .ne. "enabled") then
        WRITE(6,9132) (il,tauii(il),rhol(il),rhol_pol(il), &
           drr_gs(il),tau_neo(il),taubi(il),rhol_b(il),rhol_pol_b(il), &
           drr_gs_b(il),tau_neo_b(il),il=1,lrzmax)
      endif

 9132 format(/, &
         "Neoclassical related output:tauii,Larmor rad,pol Lar rad,",/ &
        ,"   Galeev-Sagdeev rad diffusion coeff, Neocl confmnt time",/, &
         "(First 5 real numbers are for thermal particles, then beam)" &
         ,/,"(Beam mass from first ion, Using beam energy= 80 keV.)", &
          /," lr",4x,"tauii(l)",7x,"rhol",7x,"rhol_pol", &
           6x,"drr_gs",6x,"tau_neo",7x,"taubi",8x,"rhol_b", &
           5     x,"rhol_pol_b",4x,"drr_gs_b",3x," tau_neo_b", &
           /,(i3,1p10e13.4))


!.......................................................................
!l    1.4 Meshes y and x, if nloutp1(4)=.T.
!.......................................................................

      if (.not. nlotp1(4)) go to 149

!     y mesh for l_=1 and l_=lrors
      WRITE(6,9140) (l,iy_(l),iyh_(l),itl_(l),itu_(l),l=1,lrors)
      WRITE(6,9141)
      do 140 i=1,iy_(1),10
        WRITE(6,9142) (y(ii,1),ii=i,min(i+9,iy_(1)))
        WRITE(6,9143) (dyp5(ii,1),ii=i,min(i+8,iy_(1)-1))
 140  continue
      if (lrors.gt.1) then
!        WRITE(6,9144) lrors
        do 141 i=1,iy_(lrors),10
          WRITE(6,9142) (y(ii,lrors),ii=i,min(i+9,iy_(lrors)))
          WRITE(6,9143) (dyp5(ii,lrors),ii=i,min(i+8,iy_(lrors)-1))
 141    continue
      endif

!     x mesh
      WRITE(6,9145) jx,vnorm
      do 142 j=1,jx,10
        WRITE(6,9142) (x(jj),jj=j,min(j+9,jx))
        WRITE(6,9143) (dxp5(jj),jj=j,min(j+8,jx-1))
        !Note: dxp5(j)=x(j+1)-x(j)
 142  continue

      do j=1,jx ! YuP[2018-01-05] added - print in columns: x,dx,u/c
       WRITE(*,'(a,i7,3e13.5)')'j,x,dx,u/c=',j,x(j),dxp5(j),uoc(j)
      enddo

 149  continue
 9140 format(/," y mesh:",/,1x,7("="),//,"  l","  iy"," iyh"," itl", &
        " itu",/,(i5,4i4))
 9141 format(/," y(i,1):",/)
 9142 format(1p10e12.3)
 9143 format(6x,1p9e12.3)
 9144 format(/," y(i,",i5,"):",/)
 9145 format(/," x(1:",i5,"), vnorm=",1pe12.4," :",/)

      return

!.......................................................................
!l    2. Diagnostics called during the loop over n
!.......................................................................

 200  continue
      WRITE(6,9200) n,timet
 9200 format(///,15x,55("*"),/,15x,"*    PERIODIC OUTPUT:   n =",i5, &
        " , time =",1pe11.4,"    *" ,/,15x,55("*"),//)

!.......................................................................
!l    2.0  Coulomb logrithm for general species
!.......................................................................

      do igen=1,ngen
         WRITE(6,*)
         WRITE(6,9201) igen
         do ll=1,lrz
            call tdnflxs(ll)
            call cfpgamma
            WRITE(6,9202) lr_,(gama(igen,ispec),ispec=1,ntotal)
         enddo
      enddo
      WRITE(6,*)
 9201 format('Coulomb logrithm for general species number igen =',i2)
 9202 format(' lr_ =',i3, &
              '  ln_lambda(igen,ispec=1,ntotal) =',10(1pe10.3))


!.......................................................................
!l    2.1  density and energy at each time-step
!     Same structure as in tdpltmne.f
!.......................................................................
      if (kelecg.ne.0) then
        WRITE (6,9210) kelecg,pltvs
        WRITE(6,9211) (ztr(lrindx(l)),reden(kelecg,lrindx(l)), &
          energy(kelecg,lrindx(l)),(powrf(lrindx(l),kk),kk=1,3), &
          sorpwt(lrindx(l)),sorpwti(lrindx(l)), &
          bdre(lrindx(l)),bdrep(lrindx(l)),vfluxz(lmdpln(l)),l=1,lrz)
        WRITE(6,9212) sorpwtza,(powurf(kk),kk=0,3), &
          powurfc(0),powurfl(0)
      endif

      if (niong.ne.0) then
        WRITE (6,9210) kiong(1),pltvs
        WRITE(6,9211) (ztr(lrindx(l)),reden(kiong(1),lrindx(l)), &
          energy(kiong(1),lrindx(l)),(powrf(lrindx(l),kk),kk=1,3), &
          sorpwt(lrindx(l)),sorpwti(lrindx(l)),bdre(lrindx(l)), &
          bdrep(lrindx(l)),vfluxz(lmdpln(l)),l=1,lrz)
        WRITE(6,9212) sorpwtza,(powurf(kk),kk=0,3), &
          powurfc(0),powurfl(0)
      endif

      if (niong.ne.0) then
         do kkk=1,mrfn
         WRITE(*,*) 'mode = ',kkk, 'nharm(kkk) = ',nharm(kkk), &
                 'powurf =',powurf(kkk)
         enddo
      endif

!     Normalized current and rf powers
      if (cqlpmod.ne."enabled") then
        do 221 k=1,ngen
          do 220 l=1,lrz
!BH120221:  Added following call.  Req'd at each flux surface.
            call tdnflxs(l)
            call cfpgamma
            if (entr(k,3,l).ge. 1.e-20 .and. kelecg.ne.0) then
              fnu0=2.0/tauee(lrindx(l))
              xj(lrindx(l))=curr(k,lrindx(l))/reden(k,lrindx(l)) &
                /charge/vth(kelec,lrindx(l))
              xp(lrindx(l))=entr(k,3,l)/reden(k,lrindx(l)) &
                /vth(kelec,lrindx(l))**2/fmass(k)/fnu0*1.e7

              xe(lrindx(l))=xj(lrindx(l))/xp(lrindx(l))

              fnuc=8.*pi*reden(k,lrindx(l))*charge**4 &
                *gama(kelec,kelec)/(fmass(kelec)**2*clight**3)
              xc(lrindx(l))=curr(k,lrindx(l))/(entr(k,3,l)*1.e7) &
                /(charge/(fnuc*fmass(kelec)*clight))
            endif
 220      continue

        WRITE(6,9270)
        WRITE(6,9271) (ztr(lrindx(l)),xj(lrindx(l)), &
          xp(lrindx(l)),xe(lrindx(l)),xc(lrindx(l)), &
          l=1,lrz)
 221    continue
      endif
 9270 format(/" Normalized RF Power and Current drive",/, &
       "   In order: rho,  current (units: charge*ne*vth(kelec,lr_),",/, &
       "     power (units: ne*vth(kelec,lr_)**2*me*nu0), ",/, &
       "     efficiency (j/p) (Fisch 1978 units),  ",/, &
       "     efficiency (j/p) (e/(m*c*nuc units)),  ")
 9271 format(1p5e14.4)

!     parallel profile
      if (cqlpmod .eq. "enabled") then
        WRITE(6,9213) kelecg,lrindx(1)
        WRITE(6,9214) (sz(l),denpar(kelecg,lsindx(l)), &
          enrgypa(kelecg,lsindx(l)),l=1,lrors)
        if (ls .ge. 3) then
          zdensto=0.0
          zenrgto=0.0
          zlsto=0.0
          do 210 l=1,ls
            zdensto=zdensto+denpar(kelecg,lsindx(l))*dsz(l)
            zenrgto=zenrgto+enrgypa(kelecg,lsindx(l))*dsz(l)
            zlsto=zlsto+dsz(l)
 210      continue
          WRITE(6,9215) zdensto/zlsto,zenrgto/zlsto
        endif
      endif

 9210 format(/" Misc.:  (density, energy, power density(W/cm**3)," &
        " power(W), j/p for  ", &
        "general species (k=",i2,"))",/,1x,14("="),/,4x,a8,2x, &
        "density",5x,"energy",4x,"powrf(1)",3x,"powrf(2)",3x, &
        "powrf(3)",1x,"tot pwr den",1x,"intgr pwr",2x, &
        "J/P(fi)",4x,"J/P(fi+e)",1x,"runaway rate")
 9211 format(1p11e11.3)
 9212 format(/," RF Power integrated over radius      [W] :",1pe14.4, &
        /,"       from internal ray diag. for rf [W] :",1pe14.4, &
        /,"                               mode 1 [W] :",1pe14.4, &
        /,"                               mode 2 [W] :",1pe14.4, &
        /,"                               mode 3 [W] :",1pe14.4, &
        /,"       by collisions (from ray data)  [W] :",1pe14.4, &
        /,"       by lin. damping (from ray data)[W] :",1pe14.4)
 9213 format(/," density and energy as a function of s for ", &
        " general species k=",i2," , at lrindx=",i3," :", &
        //,8x,"s",11x,"density",7x,"energy")
 9214 format(1p3e14.4)
 9215 format(/," total per ds:",1p2e14.4)

!.......................................................................
!l    2.2  flux surf. av. curnt
!.......................................................................
!l    2.2.1 radial profile
!.......................................................................

      WRITE(6,9220) pltvs
      WRITE(6,9221) (ztr(lrindx(l)),currtz(lrindx(l)), &
        currtpz(lrindx(l)),bscurm(lrindx(l),1,1),totcurz(lrindx(l)), &
        totcurzi(lrindx(l))-totcurzi(lrindx(l-1)), &
        totcurzi(lrindx(l)),l=1,lrz) !NOTE: bscurm(l,1,1) is for e_maxw only !
      WRITE(6,9222) currtza,currtpza,bscurma,totcurza
 9220 format(//" Flux surf. avg. current densities [Amps/cm**2] " &
        "and integrated currents [Amps]:",/,1x, &
        78("="),/,7x,a8,6x,"fi",9x,"fi+e",7x,"bootstr_e",6x,"total", &
        6x,"local int.",2x,"integrated")
 9221 format(1x,1p7e13.4)
 9222 format(/," tot curr[A]:",1p3e13.4,13x,e13.4)

!.......................................................................
!l    2.2.2 radial profile of toroidal and poloidal cmpnts of current
!.......................................................................

      WRITE(6,9320) pltvs
      WRITE(6,9321) (ztr(lrindx(l)), &
        curpol(lrindx(l)),ccurpol(lrindx(l)), &
        curtor(lrindx(l)),ccurtor(lrindx(l)),l=1,lrz)
      WRITE(6,9322) ccurpol(lrzmax),ccurtor(lrzmax)
 9320 format(//" Poloidal and toroidal  current densities [Amps/cm**2] " &
        "and integrated currents [Amps]:",/,1x, &
        78("="),/,25x,a8,5x,"pol",3x,"intg. pol.",13x,"tor", &
        1x,"intg. tor.")
 9321 format(19x,1p5e13.4)
 9322 format(/," total integrat. pol. current [Amps]:",1pe13.4,10x, &
          "tor. current [Amps]:",1pe13.4)

!.......................................................................
!l    2.2.3 parallel profile
!.......................................................................

      if (cqlpmod .ne. "enabled") go to 222

      WRITE(6,9223)
      WRITE(6,9224) (sz(l),currmtz(l),currmtpz(l),currmtpz(l)/psis(l), &
        l=1,lrors)
 9223 format(//" Flux surf. avg. current density [Amps/cm**2] ", &
        "(along magnetic field):"/,1x, &
        67("="),/,29x,"s",9x,"fi",9x,"fi+e",6x,"fi+e/psis")
 9224 format(19x,1p4e13.4)

 222  continue

!.......................................................................
!l    2.3   resistivity vs. radial and parallel mesh
!.......................................................................


!       R.D. Hazeltine, F.L. Hinton, and M.N. Rosenbluth,
!         Phys. Fluids 16, 1645 (1973), Eq. 69.
!       See also, Hinton and Hazeltine, Rev. Mod. Phys. 48,
!         239 (1976), Eq. 6.35, and comments therein.
!        rovsf=1. / (1.-1.95*zeps**.5+.95*zeps)

!       Connor formula for toroidal resistivity divided by Spitzer
!         (using a model e-e collision operator),
!         J.W. Connor, R.C. Grimm, R.J. Hastie, P.M. Keeping,
!         Nucl. Fus. 13, 211 (1973):
!         rovsc= Connor resistivity/Spitzer resistivity,
!         rovsc*0.706/0.93= Connor asymptotic result as eps==>1.
!         xconn= transiting particle fraction, provided by restcon.
!
!         From CQL3D, Connor formulas turn out to be most
!         accurate.
!
!        call restcon

!     onetwo formula for resistivity
!        zi33o95=4./3.*(1.-xconn)/1.95/sqrt(zeps)
!        zokap11=(0.29+0.46/(1.08+zeff(lr_)))/0.51
!        zk033=1.46 + 0.37 / zeff(lr_)
!        zft=(zk033*sqrt(zeps)-(zk033-1.)*zeps)*zi33o95
!        zre12=zokap11 / (1. - zft)

!     Kulsrud et al (PRL 31, 690 (1973)) expression for E_Driecer
!     with the 2. in the denominator "(he) wish(es) (he) had never
!     put there" (BobH).
!     The 300. below in zelecr converts from statvolts/cm to volts/cm.
!     Kulsrud plots runaway rates versus elecfld/(2.*zelecr).
!       zelecr=300.*fmass(kelec)*vth(kelec,lr_)/(2.*charge*tauee(lr_))
!
!
!
!     If efflag="toroidal", then code electric field is
!                           assumed to be purely toroidal,
!                           varying like elecfld*(Rmag/R).
!     If efflag="parallel", then code electric field is
!                           assumed to be purely parallel,
!                           varying like elecfld*(Rmag/R).
!     restp(nch(l_),lr_)=Resistivity, calc'd from distn fnctn results
!                  =<E_phi/R>/<j_phi/R>, toroidal resistivity.
!                  Except, if efswtchn.eq."neo_hh" .and.
!                    cqlpmod.ne."enabled" ==>
!                    restp=(pol cross-section-area avg of E)/currpar
!                    and currpar is sum of Hinton-Hazeltine neoclassical
!                    current + runaway current.
!     rovs(lr_)=restp/sptzr (with Zeff.ne.1 corrections in ONETWO manual)
!     rovsn(lr_)=<E_parall*B>/<j_parall*B>/sptzr
!     elecfld(lr_)=toroidal of parallel electric field (V/cm) at Rmag,
!                  depending on efflag.
!     rovsf=Hazeltine, Hinton and Rosenbluth (1973) formula.
!     rovsc=Connor formula
!     zre12=ONETWO (H&H, 1976) formula, with collisionality and Zeff.
!     rovsc*0.706/0.93=Connor asymptotic, as eps==>1.

      WRITE(*,99)
 99   format(//)
      WRITE(*,*)'If efflag="toroidal", then code electric field is'
      WRITE(*,*)'                      assumed to be purely toroidal,'
      WRITE(*,*)'                      varying like elecfld*(Rmag/R).'
      WRITE(*,*)'If efflag="parallel", then code electric field is'
      WRITE(*,*)'                      assumed to be purely parallel,'
      WRITE(*,*)'                      varying like elecfld*(Rmag/R).'
      WRITE(*,99)
      WRITE(*,*)'Explanation of following table, for each column:'
      WRITE(*,*)'================================================'
      WRITE(*,*)'epsilon=rho/R'
      WRITE(*,*)'resist_phi=Resistivity, calc''d from distn'
      WRITE(*,*)'                             fnctn results'
      WRITE(*,*)'          =<E_phi/R>/<j_phi/R>, toroidal resistivity'
      WRITE(*,*)'           Except, if efswtchn.eq."neo_hh" .and.'
      WRITE(*,*)'           cqlpmod.ne."enabled" ==>'
      WRITE(*,*)'           restp=(pol x-section-area avg of E)/currpar'
      WRITE(*,*)'               and currpar is sum of Hinton-Hazeltine '
      WRITE(*,*)'               neoclassical current + runaway current.'
      WRITE(*,*)'           Units: statV-cm/statA = seconds'
      WRITE(*,*)'restp/sptzr=resist_phi/sptzr (sptzr incls ONETWO Zeff)'
      WRITE(*,*)'resist_neo=<E_parall*B>/<j_parall*B> (seconds)'
      WRITE(*,*)'res_neo/sptzr=<E_parall*B>/<j_parall*B>/sptzr'
      WRITE(*,*)'elecfld(lr_)=toroidal or parallel electric fld (V/cm)'
      WRITE(*,*)'             at Rmag, depending on efflag.'
      WRITE(*,*)'E/E-Driecer=elecfld/E_Driecer'
      WRITE(*,*)'Connor=Connor formula, banana regime, over sptzr'
      WRITE(*,*)'rovsc*0.706/0.93=Connor asymptotic, as eps==>1.'
      WRITE(*,*)'Zeff'
      WRITE(*,*)'ONETWO low collisionality limit (H&H, 1976)/sptzr'
      WRITE(*,*)'Check tdoutput.f for references.'

      WRITE(6,9230)
 9230 format(//," Resistivity:",/,1x,12("="),/,"  l",3x, &
        "epsilon",2x,"resist_phi",1x,"resp/sptzr",1x,"resist_neo",1x, &
        "res_neo/sptzr", &
        1x,"elecfld",1x,"E/E-Dreicer",2x,"Connor",2x &
        ,"Con.larg eps",2x,"sgm_code/BK",2x,"One2/sptzr")

      ! Braams-Karney normalization resistivity [cgs]:
      res_BK = sqrt(fmass(kelec))*4.*pi*charge**2 &
       *gama(kelec,kelec)*zeff(1)/(temp(1,0)*ergtkev)**1.5

      il_old=l_
      do 230 il=1,lrz
        call tdnflxs(lmdpln(il))
        zeps=eps(lr_)
        rovsf=1. / (1.-1.95*zeps**.5+.95*zeps)
        call restcon
        zi33o95=4./3.*(1.-xconn)/1.95/sqrt(zeps)
        zokap11=(0.29+0.46/(1.08+zeff(lr_)))/0.51
        zk033=1.46 + 0.37 / zeff(lr_)
        zft=(zk033*sqrt(zeps)-(zk033-1.)*zeps)*zi33o95
        zre12=zokap11 / (1. - zft)
        zelecr=300.*fmass(kelec)*vth(kelec,lr_)/(2.*charge*tauee(lr_))
        rovsc_hi(lr_)=rovsc(lr_)*0.706/0.93
        if (restp(nch(l_),lr_).ne.zero) then
           res_BK_ratio=res_BK/restp(nch(l_),lr_)
        else
           res_BK_ratio=zero
        endif
        WRITE(6,9231) l_,zeps,restp(nch(l_),lr_),rovs(lr_) &
          ,restnp(nch(l_),lr_),rovsn(lr_),elecfld(lr_) &
          ,elecfld(lr_)/zelecr,rovsc(lr_) &
          ,rovsc_hi(lr_),res_BK_ratio,zre12
           !!! WRITE(*,*) 'l_,lr_,nch(l_)===',l_,lr_,nch(l_)   ! l_ or lr_== radial index; nch()== time step
 230  continue
 9231 format(i3,1p9e11.3,1p1e13.5,1p1e11.3)

      WRITE(*,*)'res_BK Braams-Karney normalization resistivity [cgs]=', &
       res_BK, restp(nch(l_),1)
      if(restp(nch(l_),1).ne.0.) WRITE(*,*) &
        'Normalized conductivity: sigma_code/sigma_BK=', &
       res_BK/restp(nch(l_),1)
      do ktot1=1,ntotal
      do ktot2=1,ntotal
        WRITE(*,*)'k1,k2, ln(Lambda)==gama(k1,k2)=', &
                   ktot1,ktot2,  gama(ktot1,ktot2)
      enddo
      enddo

!%OS
      WRITE(*,99)
      WRITE(*,*)'For explanation of following table, see tdoutput.f'
      WRITE(*,*)
      WRITE(6,9238)
 9238 format(/,"epsilon",2x,"res_neo/sp",4x,"Connor",3x &
        ,"Con.largep",2x,"Zeff",2x,"One2/spt(Z)",3x,"Iconnor",2x, &
        "I33/1.95/rep",3x,"Hinton",6x,"Kim",8x,"g",5x,"(eta-1)/g")
      if(machine.eq."toroidal")then
      do 2301 il=1,lrz
        call tdnflxs(lmdpln(il))
        zeps=eps(lr_)
!     R.D. Hazeltine, F.L. Hinton, and M.N. Rosenbluth (1973).
        rovsf=1. / (1.-1.95*zeps**.5+.95*zeps)
!     connor formula
        call restcon
!     Hinton, Kim and Sauter formulae for resistivity
        call resthks(l_,lr_,lmdpln_,zreshin(lr_), &
             zreskim(lr_),zressau1,zressau2)
!     onetwo formula for resistivity
        zi33o95=4./3.*(1.-xconn)/1.95/sqrt(zeps)
        zokap11=(0.29+0.46/(1.08+zeff(lr_)))/0.51
        zk033=1.46 + 0.37 / zeff(lr_)
        zft=(zk033*sqrt(zeps)-(zk033-1.)*zeps)*zi33o95
        zre12=zokap11 / (1. - zft)
        zelecr=300.*fmass(kelec)*vth(kelec,lr_)/(2.*charge*tauee(lr_))
        zg=(1.-xconn)/xconn
        rovsc_hi(lr_)=rovsc(lr_)*0.706/0.93
        WRITE(6,9239) zeps,rovsn(lr_),rovsc(lr_) &
          ,rovsc_hi(lr_),zeff(lr_),zre12, &
          xconn,zi33o95,zreshin(lr_),zreskim(lr_),zg,(rovsn(lr_)-1.)/zg
 2301 continue
      endif
 9239 format(0pf7.5,1p3e12.4,0pf6.2,1p7e12.4)
!%OS

      call tdnflxs(il_old)
      if (cqlpmod .eq. "enabled") &
        WRITE(6,9232) (l,sz(l),rovsloc(l),sptzr(l),l=1,lrors)
 9232 format(/," Local in s resistivity divided by spitzer:",//, &
        "  l",5x,"s",6x,"res/spitz",3x,"spitzer",/,(i3,1p3e11.3))

      if (n .lt. nstop) go to 233

      if (cqlpmod .ne. "enabled") WRITE(6,9233) lrindx(1)
 9233 format(/" res/spitz. per time step: for first flux surface:" &
        ," lr_= ",i3,":",/" ========================")
      if (cqlpmod .eq. "enabled") WRITE(6,9234)lrindx(1),lsindx(1)
 9234 format(/" res/spitz. per time step: at flux surface:" &
        ," lr_= ",i3," for first orbit position ls_=",i3,":", &
        /" ========================")
      illeff=lrindx(1)
      if (cqlpmod .eq. "enabled") illeff=lsindx(1)
      do 231 l=1,nch(1),8
        WRITE(6,9235) (il,rovsp(il,illeff),il=l,min(l+7,nch(1)))
        !-YuP: added: evolution of resistivity for first flux surface [cgs]
        WRITE(6,9235) (il,restp(il,1) ,il=l,min(l+7,nch(1)))
 231  continue
 9235 format(8(i5,":",1pe11.4))

!     Print out electric field iterations:
      if (efiter.eq."enabled") then
!BH171230:  Using efldn(,,) storage, not elecfldn(,,) for this case.
!BH171230:  This storage had been set up before, but not used.
         if (cqlpmod .ne. "enabled") WRITE(6,9247) lrindx(1)
 9247 format(/" iteration elecn per time step: for 1st flux surface:" &
        ," lr_= ",i3,":",/" ========================")
         do 234 niter=1,nefitera

!     But don't if beyond range of iteration:
            fiter=0.
            do 237 i=1,nch(1)
               !YuP: was a bug: fiter=fiter+abs(elecn(niter,i,illeff))
               !YuP[21-08-2017] from ampf: elecn(1:lrz,0:nstop,nefitera)
               fiter=fiter+abs(elecn(illeff,i,niter))
 237        continue
            if (fiter.eq.zero) go to 236

            do 235 l=1,nch(1),8
               !YuP[21-08-2017] from ampf: elecn(1:lrz,0:nstop,nefitera).
               !Corrected:
               WRITE(6,9235) (il,elecn(illeff,il,niter), &
                    il=l,min(l+7,nch(1)))
 235        continue
 234     continue
 236     continue
      endif

      if (lrors .eq. 1) go to 233

      ilprin=lrors/2+1
      if (cqlpmod .ne. "enabled") WRITE(6,9236) ilprin,lrindx(ilprin)
 9236 format(/" res/spitz. per time step: for surface: ", &
        "lrindx(",i3,")=",i3," :",/" ========================")
      if (cqlpmod .eq. "enabled") WRITE(6,9237) ilprin,lsindx(ilprin)
 9237 format(/" res/spitz. per time step: for orbit pos.: ", &
        "lsindx(",i3,")=",i3," :",/" ========================")
      illeff=lrindx(ilprin)
      if (cqlpmod .eq. "enabled") illeff=lsindx(ilprin)
      do 232 l=1,nch(ilprin),8
        WRITE(6,9235) (il,rovsp(il,illeff),il=l,min(l+7,nch(ilprin)))
 232  continue

 233  continue




!.......................................................................
!l    2.4   conservation diagnostic
!.......................................................................

      if (lrors .gt. 1) then
        WRITE(6,9240)
        do 240 l=1,lrors,8
          WRITE(6,9235) (il,consnp(nch(il),il),il=l,min(l+7,lrors))
 240    continue
      else if (nch(1) .ge. nstop) then
        WRITE(6,9241)
        do 241 l=1,nch(1),8
          WRITE(6,9235) (it,consnp(it,1),it=l,min(l+7,nch(1)))
 241    continue
      endif
 9240 format(//," Conservation diagnostic: per flux surface (should " &
        ,"yield machine accuracy in absence of sources/losses)",/, &
        1x,24("="),/)
 9241 format(//," Conservation diagnostic:  per time-step (should " &
        ,"yield machine accuracy in absence of sources/losses)",/, &
        1x,24("="),/)

!     transport conservation diagnostic

      if (transp.eq."enabled" .and. cqlpmod.ne."enabled") then
        WRITE(6,9242)
        do 242 l=1,lrors,8
          WRITE(6,9235) (il,constp(nch(il),il),il=l,min(l+7,lrors))
 242    continue
        call tdtrflx
        WRITE(6,9243) flxout
      endif
 9242 format(//," Transport conserv. diagnostic: per flux surface" &
        ,/,1x,9("="),/)
 9243 format(" outward flux of particles: ",1pe12.4/)

!.......................................................................
!l    2.5   component by component power flow
!     (from pltpower)
!.......................................................................

      do 250 k=1,ngen
        WRITE(6,9250) k,(l,entr(k,-1,l),entr(k,2,l) &
          ,entr(k,1,l),entr(k,3,l) &
          ,entr(k,5,l),entr(k,6,l),entr(k,7,l) &
          ,entr(k,8,l),entr(k,11,l),entr(k,10,l) &
          ,entr(k,12,l),entr(k,4,l),entr(k,0,l) &
          ,l=1,lrors)
 250  continue
 9250 format(//," Component by component power flow:",5x, &
        "general species nr: ",i2,/,1x,34("="),//, &
        " l  coll. to    Ohmic   coll. to  RF drive  ion part. ", &
        "loss-      losses    runaway synchrot. set negat. phenomen." &
        ,"    total",/ &
        "    Maxw e/ion  drive   gene. sp.            source   ", &
        "lossmode   torloss   losses  radiation f to zero enrgy loss" &
        ,/,(i2,1pe11.3,10e10.3,e11.3,/,1pe13.3))

!.......................................................................
!l    2.6 Soft Xray diagnostic (from tdsxrplt)
!.......................................................................

      if (softxry.eq."disabled" .or. kelecg.le.0) go to 260

      ienmax=inegsxr(1)
      do 261 iview=2,nv
        if (inegsxr(iview) .gt. ienmax) ienmax=inegsxr(iview)
 261  continue
!     at most 20 points printed
      ishift=ienmax/20
      if (ishift .le. 0) ishift=1
      if (nv .le. 9) WRITE(6,9261) nv
      if (nv .gt. 9) WRITE(6,92612) nv
      do 262 ien=1,ienmax,ishift
        WRITE(6,9262) ien,en_(ien),(eflux(ien,iview),iview=1,min(nv,9))
 262  continue
      if ((ienmax/ishift)*ishift .ne. ienmax) WRITE(6,9262) &
        ienmax,en_(ienmax),(eflux(ienmax,iview),iview=1,min(nv,9))
      WRITE(6,9264) (efluxt(iview),iview=1,min(nv,9))

      if (nv .gt. 9) then
        WRITE(6,9263) nv
        do 263 ien=1,ienmax,ishift
          WRITE(6,9262) ien,en_(ien),(eflux(ien,iview), &
            iview=10,min(nv,18))
 263    continue
        if ((ienmax/ishift)*ishift .ne. ienmax) WRITE(6,9262) &
          ienmax,en_(ienmax),(eflux(ienmax,iview),iview=10,min(nv,18))
      endif
      if (nv.ge.10) &
          WRITE(6,9264) (efluxt(iview),iview=10,min(nv,18))
 9261 format(//" SXR energy flux versus photon energy",/, &
        " ====================================",//, &
        "  i "," photon [keV]",i5," view chords energy", &
        " flux in [ergs/cm**2/sec/ster/eV]")
92612 format(//" SXR energy flux versus photon energy",/, &
        " ====================================",//, &
        "  i "," photon [keV]","   first 9 out of ",i2, &
        " view chords energy flux in [ergs/cm**2/sec/ster/eV]")
 9262 format(i3,1x,1p10e12.3)
 9263 format(/,"  i "," photon [keV]","   10th to ",i2,"th", &
        " view chords energy flux in [ergs/cm**2/sec/ster/eV]")
 9264 format(/,5x,"Total flux:",1p9e12.3)

!.......................................................................

 260  continue

      return
      end
end module tdoutput_mod
