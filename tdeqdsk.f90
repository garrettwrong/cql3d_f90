!
!
module tdeqdsk_mod
  
  use iso_c_binding, only : c_double
  !     used by tdtscinp
  real(c_double),public  :: psimago,psilimo
  save

contains

  subroutine tdeqdsk
    use param_mod
    use comm_mod
    use equilib_mod, only : ncoila, nccoila, pcvac, ccoil, ncoil
    use r8subs_mod, only : dcopy
    
    implicit integer (i-n), real*8 (a-h,o-z)
    save
    
    parameter (niterate=20)
    
    dimension fpsiareq(nnra),ffpareq(nnra),prareq(nnra), &
   ppareq(nnra),q(nnra),eqpsi_(nconteqa),d2q_(nconteqa)
      parameter(itlrza=3*lrza+1+3*nrada+3*nconteqa)
      dimension workk(itlrza),atw(ntotala),dpsiar(nrada),d2gpary(nrada), &
        d2pary(nrada),d2ppary(nrada),ene(lrza),te(lrza), &
        ti(lrza,ntotala),enp(lrza,ntotala),eni(lrza,ntotala), &
        z1(lrza,ntotala),dumm(400),xcontr(50),ycontr(50), &
        psiart(nrada),xlimiter(50),ylimiter(50)
      dimension zeqpsir(lrza)
      character(len=5) :: blanks
      data mode /0/
      data blanks/"     "/

!.......................................................................
!     Called from tdinitl, if partner.eq."bram",
!        and from tdtloop, if n.ge.nstop.
!.......................................................................


      if (eqmod.ne."enabled") return

!.......................................................................
!     pick equilpsp values according to lrindx mesh
!.......................................................................
      do 20 l=1,lrz
        zeqpsir(l)=equilpsp(lrindx(l))
20   enddo
!..................................................................
!     file (the file used to reinitialize the MHD equilibrium code).
!     The procedure is as follows (1) create the arrays on the flux
!     surface grid (2) interpolate to the MHD equilibrium grid (
!     required for eqdsk) and (3) rescale to get all quantities into
!     MKS.
!     The procedure differs depending on the state of variable
!     partner.
!..................................................................

      !ncontr=0    !-YuP: commented out. Why needed here? defined in equilib
      !nlimiter=0  !-YuP: commented out. Why needed here?

      if (partner.eq."selene") then

!..................................................................
!     Compute the derivative of Pressure  w.r.t. psi
!..................................................................

        ipestg=3
        call firstdrv(equilpsp(1),prest,prestp,lrzmax)

!..................................................................
!     Compute <JparB> (volume, or flux surface average).
!     Compute the bootstrap effect too, <JbsB>. The bootstrap
!     current (bscurm) is assumed to be area averaged.
!
!     Presently, only electron bootstrap contribution included 
!                for jhirsh=88/99. (BobH, 990823).
!
!     Input currents are in Amps/cm**2
!..................................................................

        do 5 ll=1,lrz
          ilr=lrindx(ll)
          jparb(ll)=0.0
          jparbp(ll)=0.0
          if (zmaxpsi(ilr).ne.zero) then
            jparb(ll)=currmtz(lmdpln(ll))*bmdplne(ilr)*psidz(ilr)/ &
              zmaxpsi(ilr)*3.e+9
            jparbp(ll)=currmtpz(lmdpln(ll))*bmdplne(ilr)*psidz(ilr)/ &
              zmaxpsi(ilr)*3.e+9
          endif
          jparbt(ll)=bscurm(ilr,1,1)*btor*3.e+9
5      enddo

!..................................................................
!     Procedure for computing f**2 is iterative. We presume f is
!     known at the magnetic axis (boundary condition) from the
!     previous call to the equilibrium code and integrate out.
!     equilpsp=-equilpsi
!..................................................................

        factr=-8.*pi/clight
        do 101 ll=1,lrz
          ilr=lrindx(ll)
          tr1(ll)=bmdplne(ilr)**2*psidz(ilr)/zmaxpsi(ilr)
 101    enddo
        if (mode.ne.0) then
!BH070116          fpsiz2(0)=fpsiar(nnv)**2
          fpsiz2(0)=fpsiar(nfp)**2
          do 100 nit=1,niterate
            do 120 ll=1,lrz
              ilr=lrindx(ll)
              fpsiz2(ll)=fpsiz2(ll-1)+(equilpsp(ilr)-equilpsp(ilr-1))* &
                factr*(1./tr1(ll))*(fpsiz(ll)*(jparbp(ll)+jparbt(ll))+ &
                clight*fpsiz(ll)**2*prestp(ilr))
120        end do

!..................................................................
!     f
!..................................................................


            do 130 ll=1,lrz
              fpsiz(ll)=sqrt(fpsiz2(ll))
130        end do
100     end do
        endif

!..................................................................
!     ff'
!..................................................................

        do 140 ll=1,lrz
          ilr=lrindx(ll)
          ffpsiz(ll)=(1./tr1(ll))*factr*(fpsiz(ll)*(jparbp(ll)+ &
            jparbt(ll))+clight*fpsiz(ll)**2*prestp(ilr))*.5
140    end do

!..................................................................
!     Compute the toroidal current, toteqd
!..................................................................

        toteqd=0.
        do 150 ll=1,lrz
          ilr=lrindx(ll)
          toteqd=toteqd-clight/(8.*pi**2)*ffpsiz(ll)*onovrpz(ilr,2) &
            *dvol(ilr)-clight/twopi*prestp(ilr)*dvol(ilr)
150    end do

!..................................................................
!     The eqdsk format requires that all the quantities which are
!     functions of psi be written on an evenly spaced grid in psi.
!     This will be done now. We have psilim=0
!     The psi array to which we wish to interpolate was determined in
!     subroutine equilib - psiar
!..................................................................

        i1p(1)=4
        i1p(2)=4

!..................................................................
!     ff' and f
!..................................................................

        call coeff1(lrz,zeqpsir(1),fpsiz,d2fpsiz,i1p,1,workk)
        call coeff1(lrz,zeqpsir(1),ffpsiz,d2ffpsiz,i1p,1,workk)
        itab(1)=1
        itab(2)=0
        itab(3)=0

!..................................................................
!     Pressure and derivative
!     psiar defined in equilib - rescale so that min is at mag axis;
!     max (=0 is at magnetic axis).
!..................................................................

!BH070116        call dcopy(nnv,psiar,1,dummyar,1)
        call dcopy(nfp,psiar,1,dummyar,1)
!BH070116        do 60 ne=1,nnv
        do 60 ne=1,nfp
!BH070116          psiar(ne)=-dummyar(nnv+1-ne)
          psiar(ne)=-dummyar(nfp+1-ne)
60     end do
        do 65 i=2,nconteqn
          eqpsi_(i)=-eqpsi(i)
65     end do
        call coeff1(nconteqn-1,eqpsi_(2),q_(2),d2q_(2),i1p,1,workk)
!BH070116        do 50 ne=1,nnv
        do 50 ne=1,nfp
          call terp1(lrzmax,equilpsp(1),prest,d2prest,psiar(ne),1,tab, &
            itab)
          prar(ne)=tab(1)
          if (prar(ne).lt.0.) prar(ne)=0.
          call terp1(lrz,zeqpsir(1),fpsiz,d2fpsiz,psiar(ne),1,tab, &
            itab)
          fpsiar(ne)=tab(1)
          call terp1(lrz,zeqpsir(1),ffpsiz,d2ffpsiz,psiar(ne),1,tab, &
            itab)
          ffpar(ne)=tab(1)
          call terp1(nconteqn-1,eqpsi_(2),q_(2),d2q_(2),psiar(ne),1,tab, &
            itab)
          q(ne)=tab(1)
50     end do
        ffpar(1)=ffpsiz(1)
!BH070116        call firstdrv(psiar(1),prar(1),ppar(1),nnv)
        call firstdrv(psiar(1),prar(1),ppar(1),nfp)

      else if (eqsource.eq."eqdsk" .or. eqsource.eq."topeol" ) then

!BH070116        nnv_=nnv
!BH070116        do 201 j=1,nnv
!BH070116          fpsiareq(j)=fpsiar(nnv+1-j)
!BH070116          prareq(j)=prar(nnv+1-j)
!BH070116          ppareq(j)=ppar(nnv+1-j)
       nfp_=nfp
        do 201 j=1,nfp
          fpsiareq(j)=fpsiar(nfp+1-j)
          prareq(j)=prar(nfp+1-j)
          ppareq(j)=ppar(nfp+1-j)
201    end do

      endif   ! endif on (partner.eq.'selene' .or. eqsource.eq.'eqdsk')

!..................................................................
!     Convert as appropriate to MKS
!..................................................................

!BH070116      do 200 i=1,nnv
      do 200 i=1,nfp
        fpsiareq(i)=fpsiareq(i)/1.e+6
        ffpareq(i)=ffpar(i)*1.e+2   !Previously 1e-4, so check if
                                    !using partner='selene', BH010226.
        prareq(i)=prareq(i)/10.   
        ppareq(i)=ppareq(i)*1.e+7
200  end do


!..................................................................
!     Now for the case partner="bramb"
!..................................................................

      if (eqsource.eq."tsc") then
!BH070116        nnv_=nnr
        nfp_=nnr
        dpsi_=(psilimo-psimago)/(nnr-1)
        do 300 i=1,nnr
          psiart(i)=psimago+(i-1)*dpsi_
300    end do
        do 307 i=2,nconteqn
          eqpsi_(i)=-eqpsi(i)*1.e-8+psilimo
307    end do

!..................................................................
!     The radial arrays from the TSC code are defined on an arbitrary
!     rho mesh. To create a GA style eqdsk file, we have to interpolate
!     onto the psiart psi mesh (defined above). We remain in MKS units.
!     The temperatures are in KeV.
!..................................................................

        i1p(1)=4
        i1p(2)=4
        call coeff1(npsitm,psiar_,fpsiar_,dpsiar,i1p,1,workk)
        call coeff1(npsitm,psiar_,gpary,d2gpary,i1p,1,workk)
        call coeff1(npsitm,psiar_,pary,d2pary,i1p,1,workk)
        call coeff1(npsitm,psiar_,ppary,d2ppary,i1p,1,workk)
        call coeff1(nconteqn-1,eqpsi_(2),q_(2),d2q_(2),i1p,1,workk)
        itab(1)=1
        itab(2)=0
        itab(3)=0
        do 301 i=1,nnr
          call terp1(npsitm,psiar_,fpsiar_,dpsiar,psiart(i),1,tab,itab)
          fpsiareq(i)=tab(1)
          call terp1(npsitm,psiar_,gpary,d2gpary,psiart(i),1,tab,itab)
          ffpareq(i)=tab(1)
          call terp1(npsitm,psiar_,pary,d2pary,psiart(i),1,tab,itab)
          prareq(i)=tab(1)
          call terp1(npsitm,psiar_,ppary,d2ppary,psiart(i),1,tab,itab)
          ppareq(i)=tab(1)
          call terp1(nconteqn-1,eqpsi_(2),q_(2),d2q_(2),psiart(i),1,tab &
            ,itab)
          q(i)=tab(1)
301    end do

!..................................................................
!     Now set up the plasma arrays that will also be output
!..................................................................

        nprim=nspc
        write(6,'(/A/)') ' WARNING: nimp not defined in comm.h'
         nion=nprim+nimp
      else if (eqsource.eq."eqdsk" .or. eqsource.eq."topeol") then
        ku=0
        if (kelecm.ne.0) ku=ku+1
        if (kelecg.ne.0) ku=ku+1
        nion=ntotal-ku
        nprim=nion
      endif   ! endif on (partner.eq.'tsc' .or. eqsource.eq.'eqdsk')

      nti=1
      nimp=0
      ku=0
      do 302 k=1,ntotal
        if (k.ne.kelecg .and. k.ne.kelecm) then
          ku=ku+1
          atw(ku)=fmass(k)/1.6724e-24
          do 303 l=1,lrzmax
            z1(l,ku)=bnumb(k)
303      end do
        endif
302  end do
      k=kelecm
      if (k.eq.0) k=kelecg
      do 304 l=1,lrzmax
        ene(l)=reden(k,l)*1.e6
        te(l)=temp(k,l)
304  end do
      ku=0
      do 305 k=1,ntotal
        if (k.ne.kelecg .and. k.ne.kelecm) then
          ku=ku+1
          do 306 l=1,lrzmax
            enp(l,ku)=reden(k,l)*1.e6
            ti(l,ku)=temp(k,l)
306      end do
        endif
305  end do

!..................................................................
     !     Continue conversion to end 
!..................................................................

      psilim=-psilim/1.e+8
      do 205 i=1,nnz
        do 206 j=1,nnr
          epsi(j,i)=-epsi(j,i)/1.e+8
206    end do
205 end do
      psimag=-psimag/1.e+8
      btor=btor/1.e+4
      rbox=rbox/1.e+2
      zbox=zbox/1.e+2
      rboxdst=rboxdst/1.e+2
      radmaj=radmaj/1.e+2
      raxis=rmag/1.e+2
      zaxis=0.
      toteqd=toteqd/3.e+9

!..................................................................
!     Disk writes
!..................................................................

      do i=1,5     ! for dummy data below
         tem1(i)=0.
      enddo

      if (eqsource.eq."eqdsk" .or. eqsource.eq."topeol") then
        if (eqdskalt.ne."enabled") go to 800
        open(unit=17,file="tdeqdsk",delim='apostrophe', &
             status='unknown')
      elseif (eqsource.eq."tsc") then
        open(unit=17,file="tdeqdsk",delim='apostrophe',status='unknown')
      endif
!BH:000906  Maybe problem with format of the following write(17,210)
!BH070116      write(17,210) mnemonic,ipestg,nnr,nnz,nnv_
!BH070116:  Swe ymideqd=0.  BE CAREFUL in future something else needed.
      ymideqd=0.
      write(17,210) mnemonic,ipestg,nnr,nnz,nfp_
      write(17,220) rbox,zbox,radmaj,rboxdst,ymideqd
      write(17,220) raxis,zaxis,psimag,psilim,btor
      write(17,220) toteqd,(tem1(i),i=1,4)
      write(17,220) (tem1(i),i=1,5)
!BH070116      write(17,220) (fpsiareq(i),i=1,nnv_)
!BH070116      write(17,220) (prareq(i),i=1,nnv_)
!BH070116      write(17,220) (ffpareq(i),i=1,nnv_)
!BH070116      write(17,220) (ppareq(i),i=1,nnv_)
      write(17,220) (fpsiareq(i),i=1,nfp_)
      write(17,220) (prareq(i),i=1,nfp_)
      write(17,220) (ffpareq(i),i=1,nfp_)
      write(17,220) (ppareq(i),i=1,nfp_)
      write(17,220) ((epsi(i,j),i=1,nnr),j=1,nnz)
      if (partner.ne."selene" .or. partner.eq."bram") then
!BH070116        write(17,220) (q(i),i=1,nnv_)
        write(17,220) (q(i),i=1,nfp_)
      endif
      if (partner.eq."selene") then
        write(17,221) (ncoil(i),i=1,5)
        do 230 i=1,5
          if(ncoil(i).le.0) go to 230
          write(17,250) (ccoil(nn,i),nn=1,ncoil(i))
230    end do
        write(17,250) (pcvac(i),i=1,9)
      endif
      if (partner.ne."selene") then
         write(*,*)  &
        'tdeqdsk: work reqd to write xcontr,ycontr,xlimiter,ylimiter'
!        write(17,221) ncontr,nlimiter
!        if (ncontr.ne.0) then
!          write(17,220)(xcontr(i),ycontr(i),i=1,ncontr)
!          write(17,220)(xlimiter(i),ylimiter(i),i=1,nlimiter)
!        endif

!..................................................................
!     The following output is not part of the "official" eqdsk. It is
!     what distinguishes eqdsk created in onetwo from those created by
!     EFIT.
!     we put out enough information so that the eqdsk could be used as
!     a restart file.
!
!     lrzmax  = number of radial mesh points
!     nprim= number of primary ion species
!     nimp = number of impurities
!     nti  = number of species for which ion temperatures given
!     (1, for lumped ion species).
!     (nion=nprim+nimp)
!     atw  = atomic mass numbers of ions
!     z1   = charge state of each ion (can vary as function or radius).
!     zeffz= radial variation of z_effective
!     1.e-2*radmin*rya(j) = radial rho mesh from mag. axis to edge.
!     -equilpsi(j)*1.e-8= corresponding psi values (MKS) from magnetic
!     axis to edge.
!     te   = electron temperature (keV) on radial mesh.
!     ti   = ion temperature (keV) on radial mesh.
!     ene  =  electron density (MKS)
!     enp  =  primary (and impurity if nimp.ne.0)ion species density (MKS)
!..................................................................

        nti=1
        write(17,221) lrzmax,nprim,nimp,nti,lrors,lrz
        write (17,220) (atw(k),k=1,nion)
        do 700  k=1,nion
 700    write(17,220) (z1(j,k),j=1,lrzmax)
        write(17,220) (zeff(j),j=1,lrzmax)
        do 705  j=1,lrzmax
 705    dumm(j)=1.e-2*radmin*rya(j)
        write(17,220) (dumm(j),j=1,lrzmax)
        do 706 j=1,lrzmax
          dumm(j)=-equilpsi(j)*1.e-8
706    end do
        write(17,220) (dumm(j),j=1,lrzmax)
        write(17,220) (te(j),j=1,lrzmax)
        do 710  k=1,nti
 710    write(17,220) (ti(j,k),j=1,lrzmax)
        write(17,220)  (ene(j),j=1,lrzmax)
        do 730 k=1,nprim
 730    write(17,220) (enp(j,k),j=1,lrzmax)
        if (nimp.ne.0)  then
          do 740 k=1,nimp
 740      write(17,220) (enp(j,k),j=1,lrzmax)
        endif
      endif
      close(unit=17)

 800  continue

!..................................................................
!     Convert back to cgs where necessary..
!..................................................................

      psilim=-psilim*1.e8
      psimag=-psimag*1.e+8
      do 809 j=1,nnz
        do 810 i=1,nnr
          epsi(i,j)=-epsi(i,j)*1.e+8
810    end do
809 end do
    btor=btor*1.e+4
    rbox=rbox*1.e+2
    zbox=zbox*1.e+2
    rboxdst=rboxdst*1.e+2
    radmaj=radmaj*1.e+2
    toteqd=toteqd*3.e+9
    raxis=rmag*1.e+2
    zaxis=0.
250 format( (5(e21.14)))
210 format(a48,4i4)
220 format(5e16.9)
221 format(6i5)
    return
  end subroutine

end module
