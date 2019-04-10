c
c
      subroutine eqfndpsi(psides,areades,volum)
      use param_mod
      use cqcomm_mod
      implicit integer (i-n), real*8 (a-h,o-z)
CMPIINSERT_INCLUDE     

c     Set parameter for trapped particle frac calc (150 in ONETWO).
      parameter (nlam=150)
      dimension suml(nlam)


c..................................................................
c     This subroutine is called from subroutine eqcoord.
c
c     If 0 .le. rovera(lr_) .le. 1. then:
c     This routine does a Newton's iteration to determine the
c     value of the equilibrium psi, psides, associated with radial
c     coordinate (passed in common) erhocon(lr_)=rovera(lr_)*rhomax:
c     rhomax=max[rya,1.e-8], rya() is namelist input.
c     rhomax is determined in eqrhopsi for the radial coord
c     choice specified by namelist variable radcoord.
cBH090811:  Actually, from subroutine eqrhopsi, rhomax is obtained
cBH090811:  from a linear extrapolation of eqrho(eqpsi) to the
cBH090811:  eqdsk value psilim (eqsource='eqdsk').  This gives
cBH090811:  a more accurate value of rhomax, evidently for increased
cBH090811:  accuracy.
c     
c     The arrays eqrho(j), eqpsi(j) and eqfopsi(j), j=1:nconteqn, 
c     have been calculated in eqrhopsi.
c     eqpsi,eqrho are corresponding radial psi, and coord values 
c     (in accord with radcoord), and  eqfopsi=f=R*B_phi.
c
c     If (rovera(lr_).lt.0.) then:
c     Set psides=povdelp*delp and find the contour such that
c     psi=psides directly (no iterations required).
c..................................................................

CMPIINSERT_IF_RANK_EQ_0      
c      write(*,*)'eqfndpsi: radcoord,rhomax,lr_,erhocon(lr_)= ',
c     +                     radcoord,rhomax,lr_,erhocon(lr_)
CMPIINSERT_ENDIF_RANK


      if (rovera(lr_).ge.0.) then
c..................................................................
c     Begin by finding the first index jval such that eqrho(jval) is
c     larger than rhodes.
c..................................................................

        rhodes=erhocon(lr_)
        if (rhodes.gt.rhomax) call eqwrng(8)
        do 10 j=2,nconteqn
          if (rhodes.le.eqrho(j)) go to 11
 10     continue
 11     continue
        jval=j
        
CMPIINSERT_IF_RANK_EQ_0      
        WRITE(*,*)
        WRITE(*,*)'eqfndpsi: lr_,rhodes,jval,eqrho(jval-1),eqrho(jval)',
     +                       lr_,rhodes,jval,eqrho(jval-1),eqrho(jval)
c        write(*,*)'eqfndpsi: eqrho(j),j=1,nconteqn',
c     +                     (eqrho(j),j=1,nconteqn)
CMPIINSERT_ENDIF_RANK

c..................................................................
c     Begin iteration loop
c..................................................................

        psi2=eqpsi(jval)
        rho2=eqrho(jval)
        psi1=eqpsi(jval-1)
	  rho1=eqrho(jval-1)
CMPIINSERT_IF_RANK_EQ_0      
        !WRITE(*,*)
        !WRITE(*,*)'eqfndpsi: (psi2-psimag)/psimag',(psi2-psimag)/psimag
CMPIINSERT_ENDIF_RANK
	  !psimag=0 in a mirror machine.
        if(jval.le.2) then !-YuP: for better convergence near m.axis 
           psi1=psimag
           rho1=0.d0
        endif
        
        iter=0
 20     continue !-> iteration loop (through Line~175) ------------------
        !-YuP: Bi-linear interpolation -> initial guess for psinew:
        psinew= psi1 + (psi2-psi1)*
     ~ ((rhodes-rho1)/(rho2-rho1))*((rhodes+rho1)/(rho2+rho1))
        epsicon(lr_)=psinew
        eqcall="disabled"
        !--------------------
        call eqorbit(psinew) ! Get (solr_(l),solz_(l)) tracing flux surface
        !--------------------
        do 70 l1=1,lorbit_
          ! The r.h.s. are values from eqorbit:
          solr(l1,lr_)=solr_(l1)
          solz(l1,lr_)=solz_(l1)
          es(l1,lr_)=es_(l1)
          eqbpol(l1,lr_)=eqbpol_(l1) !=Bpol along field line (pol.plane)
          bpsi(l1,lr_)=bpsi_(l1)
          thtpol(l1,lr_)=thtpol_(l1)
          eqdell(l1,lr_)=eqdell_(l1)
CMPIINSERT_IF_RANK_EQ_0      
          if(eqbpol_(l1).eq.0. )then
             WRITE(*,'(a,i5,2e17.10)')
     +       'eqfndpsi: l1,eqdell_(l1),eqbpol_(l1) ',
     +                  l1,eqdell_(l1),eqbpol_(l1)
          endif
CMPIINSERT_ENDIF_RANK
 70     continue
        ! The r.h.s. (*_) are values from eqorbit:
        eqdells(lr_)=eqdells_ ! =dl along field line (pol.plane)
        lorbit(lr_)=lorbit_
        rmcon(lr_)=rmcon_
        rpcon(lr_)=rpcon_
        zmcon(lr_)=zmcon_
        zpcon(lr_)=zpcon_
        es_bmax(lr_)=es_bmax_
        bpsi_max(lr_)=bpsi_max_
        bpsi_min(lr_)=bpsi_min_
        lbpsi_max(lr_)=lbpsi_max_
        lbpsi_min(lr_)=lbpsi_min_
        bthr(lr_)=bthr_
        btoru(lr_)=btoru_
        fpsi(lr_)=fpsi_
        zmax(lr_)=zmax_
        btor0(lr_)=btor0_
        bthr0(lr_)=bthr0_
        bmidplne(lr_)=bmidplne_  !At min bpsi_ point, not necessarily
                                 !the midplane, for eqsym.eq."none"
        !Get average values over surface epsicon(lr_) (corr.to rhodes)                      
        eqorb="disabled"
        call eqvolpsi(epsicon(lr_),volum,areac)
        call eqonovrp(epsicon(lr_),onovrp1,onovrp2) !get <1/R>, <1/R**2>
        
        ! Find flux-surf. averaged values of tlorb1(l),
        ! save into bpolsqa_,flxavgd_
        ! Note: epsicon_ is INPUT, but only needed when eqorb="enabled"
        ! which is not the case here.
        do 60 l=1,lorbit_
          tlorb1(l)=eqbpol_(l)**2 !=Bpol^2 along given surface (iterated)
 60     continue
        call eqflxavg(epsicon_,tlorb1,bpolsqa_,flxavgd_)
        !OUT= bpolsqa_= SUM(tlorb1*dl/Bpol)/SUM(dl/Bpol)
        
        ! Similarly, Find flux-surf. averaged values of tlorb1(l),
        ! save into psiovr_,flxavgd_
        do 40 l=1,lorbit_
          tlorb1(l)=bpsi_(l)/solr_(l) !=(B(l)/B0)/R
 40     continue
        call eqflxavg(epsicon_,tlorb1,psiovr_,flxavgd_)
        !OUT= psiovr_= SUM(tlorb1*dl/Bpol)/SUM(dl/Bpol)
        !(and again, epsicon_ is not used here)
        
        bpolsqa(lr_)=bpolsqa_
        psiovr(lr_)=psiovr_
        flxavgd(lr_)=flxavgd_ !=SUM(dl/Bpol)
        onovrp(1,lr_)=onovrp1 ! <1/R> over lr_ surface
        onovrp(2,lr_)=onovrp2 ! <1/R^2> over lr_ surface
        
        ! Get toroidal "f(psi)" (=R*B)
        ! for psi=epsicon(lr_) surface; based on interpolation:
        call eqfpsi(epsicon(lr_),fpsi_,fppsi_)

        ! Now the value of rhonew (iterated value) can be evaluated:
        if (radcoord.eq."sqtorflx") then
        !YuP[2019] Need to inspect this part in detail.
        !YuP: Try to use Bi-linear interpolation, 
        !similar to initial guess for psinew:
        ! psinew= psi1 + (psi2-psi1)*
        !         ((rhodes-rho1)/(rho2-rho1))*((rhodes+rho1)/(rho2+rho1))
        ! where psi1=eqpsi(jval-1), psi2=eqpsi(jval)
        ! and similarly for rho1,rho2
           dvolum=(volum-eqvol(jval-1))
           fpsih=(fpsi_+eqfopsi(jval-1))*.5 !similar to eqfh below
           onok=.5*(onovrp(1,lr_)+eqovrp(jval-1,1)) ! similar to eqrph
           tem=eqrho(jval-1)**2*pi*btor ! value at previous radial point
           onoh=(onovrp(2,lr_)+eqovrp(jval-1,2))*.5
           rhonew=tem+onoh*dvolum*fpsih/pi*0.5
           !Definition from eqrhopsi:
           ! dvolum=(eqvol(j)-eqvol(j-1))
           ! eqrph=(eqovrp(j,2)+eqovrp(j-1,2))*.5
           ! eqfh=(eqfopsi(j)+eqfopsi(j-1))*.5
           ! eqrho(j)=eqrho(j-1)+eqrph*dvolum*eqfh/pi*0.5
           !And then, take sqrt, as here below:
           areanew=areac
           rhonew=sqrt(rhonew/pi/btor)
        elseif (radcoord.eq."sqarea") then
           rhonew=sqrt(areac/pi)
           areanew=areac
        elseif (radcoord.eq."sqvol") then
           rhonew=sqrt(volum/(2.*pi**2*rmag))
           areanew=areac
        elseif (radcoord.eq."rminmax") then
           rhonew=0.5*(rpcon_-rmcon_)
           areanew=areac
        elseif (radcoord.eq."polflx") then
           rhonew=(psinew-psimag)/(psilim-psimag)
           areanew=areac
        elseif (radcoord.eq."sqpolflx") then
           rhonew=sqrt((psinew-psimag)/(psilim-psimag))
           areanew=areac
        endif
         
        !Compare the given rhonew (iterated) with target value of rhodes:    
        err=abs(rhonew-rhodes)/rhodes
        iter=iter+1 ! count iterations; usually 2-4 is sufficient

CMPIINSERT_IF_RANK_EQ_0      
        WRITE(*,'(a,2i5,2f16.10,3e17.7,e12.3)')
     +  'eqfndpsi: iter,lorbit_,rhonew,rhodes,psinew,R(1),R(l),err',
     +        iter,lorbit_,rhonew,rhodes,psinew,solr_(1),solr_(l1),err
CMPIINSERT_ENDIF_RANK

        if (err.gt.1.e-5 .and. iter.lt.35) then !max number of iter: was 25
          if (rhonew.gt.rhodes) then
            rho2=rhonew
            psi2=psinew
          else
            rho1=rhonew
            psi1=psinew
            !if(iter.gt.25)then
              !Poor convergence (usually near rho=0):
              ! try to reset rho2 and psi2 :
            !  rho2=0.d0
            !  psi2=psimag
            !endif
          endif
          go to 20  !  GO TO NEXT ITERATION
        endif
        
        if (err.gt.1.e-2) then

CMPIINSERT_IF_RANK_EQ_0      
           WRITE(*,'(a,i4,e12.3,i4,e12.3)')
     +      'eqfndpsi/WARNING: POOR CONVERG. lr,rhonew,iter,err=',
     +                                       lr_,rhonew,iter,err
CMPIINSERT_ENDIF_RANK

           !!!YuP call eqwrng(4) !-> will stop the job
           !Sometimes the convergence is poor near magn.axis,
           !presumably because PSI(rho) is nearly flat.
        endif
        psides=psinew
        areades=areanew

        fppsi(lr_)=fppsi_
        call eqppsi(epsicon(lr_),ppsi_,pppsi_)
        pppsi(lr_)=pppsi_
        

c.......................................................................
c     Now for the case that rovera(lr_)=.lt.0.
c.......................................................................

      else  !rovera(lr_).lt.0
        delp=(psimag-psilim)
        psides=psimag-povdelp*delp
        do 30 j=2,nconteqn
          if (eqpsi(j).lt.psides) go to 31
 30     continue
 31     continue
        jval=j
        epsicon(lr_)=psides
        call eqorbit(psides) ! Get (solr_(l),solz_(l)) tracing flux surface
        do 90 l1=1,lorbit_
          solr(l1,lr_)=solr_(l1)
          solz(l1,lr_)=solz_(l1)
          es(l1,lr_)=es_(l1)
          eqbpol(l1,lr_)=eqbpol_(l1)
          bpsi(l1,lr_)=bpsi_(l1)
          thtpol(l1,lr_)=thtpol_(l1)
          eqdell(l1,lr_)=eqdell_(l1)
          bpsi(l1,lr_)=bpsi_(l1)
 90     continue
        eqdells(lr_)=eqdells_
        lorbit(lr_)=lorbit_
        rmcon(lr_)=rmcon_
        rpcon(lr_)=rpcon_
        bthr(lr_)=bthr_
        btoru(lr_)=btoru_
        fpsi(lr_)=fpsi_
        zmax(lr_)=zmax_
        btor0(lr_)=btor0_
        bthr0(lr_)=bthr0_
        eqorb="disabled"
        call eqvolpsi(epsicon(lr_),volum,areac)
        call eqonovrp(epsicon(lr_),onovrp1,onovrp2)
        do 50 l=1,lorbit_
          tlorb1(l)=bpsi_(l)/solr_(l)
 50     continue
        call eqflxavg(epsicon_,tlorb1,psiovr_,flxavgd_)
        do 80 l=1,lorbit_
          tlorb1(l)=eqbpol_(l)**2
 80     continue
        call eqflxavg(epsicon_,tlorb1,bpolsqa_,flxavgd_)
        psiovr(lr_)=psiovr_
        flxavgd(lr_)=flxavgd_
        onovrp(1,lr_)=onovrp1
        onovrp(2,lr_)=onovrp2
        call eqfpsi(epsicon(lr_),fpsi_,fppsi_)

        if (radcoord.eq."sqtorflx") then
           fpsih=(fpsi_+eqfopsi(jval-1))*.5
           onok=.5*(onovrp(1,lr_)+eqovrp(jval-1,1))
           tem=eqrho(jval-1)**2*pi*btor
           onoh=(onovrp(2,lr_)+eqovrp(jval-1,2))*.5
           rhonew=tem+onoh*(volum-eqvol(jval-1))*fpsih/pi*0.5
           areanew=areac
           rhonew=sqrt(rhonew/pi/btor)
        elseif (radcoord.eq."sqarea") then
           rhonew=sqrt(areac/pi)
           areanew=areac
        elseif (radcoord.eq."sqvol") then
           rhonew=sqrt(volum/(2.*pi**2*rmag))
           areanew=areac
        elseif (radcoord.eq."rminmax") then
           rhonew=0.5*(rpcon_-rmcon_)
           areanew=areac
        elseif (radcoord.eq."polflx") then
           rhonew=psides
           areanew=areac
        elseif (radcoord.eq."sqpolflx") then
           rhonew=sqrt(psides)
           areanew=areac
        endif
              
        areades=areanew
        psinew=psides
        erhocon(lr_)=rhonew
        rovera(lr_)=erhocon(lr_)/rhomax

        fppsi(lr_)=fppsi_
        call eqppsi(epsicon(lr_),ppsi_,pppsi_)
        pppsi(lr_)=pppsi_

      endif  !on rovera(lr_).ge./.le. 0.
c
c.......................................................................
c     compute <bpsi> and <bpsi**2>,<1/(bpsi*R**3>
c.......................................................................
c
      do 41 l=1,lorbit_
         tlorb1(l)= bpsi_(l)
         tlorb2(l)= bpsi_(l)**2
 41   continue
      call eqflxavg(epsicon_,tlorb1,zpsiavg,flxavgd_)
      call eqflxavg(epsicon_,tlorb2,zpsi2av,flxavgd_)
      psiavg(1,lr_)=zpsiavg
      psiavg(2,lr_)=zpsi2av
      do l=1,lorbit_
         tlorb1(l)=1./(bpsi_(l)*solr_(l)**3)
      enddo
      call eqflxavg(epsicon_,tlorb1,zpsiavg,flxavgd_)
      onovpsir3(lr_)=zpsiavg

c.......................................................................
c     Calculate effective trapped particle fraction 
c     (see e.g., Hirshman and Sigmar, Nucl. Fus. 21, 1079 (1981), 
c      Eq. 4.54)
c     trapfrac=1.-0.75*<B**2>*integral[0,Bmax**-1]{lambda*dlambda/
c                                        <(1-lambda*B)**0.5>}
c     Manipulated (BH) this to following (with zeta**0.5=Bmax*lambda)
c     to be close to a ONETWO expression, and integrated as below:
c     trapfrac=1-0.75*<h**2>*integral[0,1]{d_zeta/
c                                         <(2.*sqrt(1-zeta**.5*h)>
c      where <...> is flux surface avg, h=B/Bmax.
c.......................................................................

      !bmaxbmin=bpsi(lorbit_,lr_) ! YuP: maybe bpsi(lbpsi_max(lr_),lr_) ?
      ! For a general case of eqsym:
      bmaxbmin=bpsi_max(lr_) ! YuP [July 2014] ==Bmax/Bmin
      
      h2fsa=psiavg(2,lr_)/bmaxbmin**2
      
      dlam=1./(nlam-1.)
      do ilam=1,nlam
         rtlam=sqrt((ilam-1)*dlam)
         do l=1,lorbit_
            hlam=bpsi_(l)/bmaxbmin ! = B(l)/Bmax
            val=abs(1.0-rtlam*hlam)
            tlorb1(l)=sqrt(val)
         enddo
         call eqflxavg(epsicon_,tlorb1,suml(ilam),flxavgd_)
      enddo
         
      xi0=0.
      do ilam=1,nlam-1
         xi0=xi0+0.25*(1.0/suml(ilam) + 1.0/suml(ilam+1))*dlam
      enddo

      trapfrac(lr_)=1.-0.75*h2fsa*xi0
      
CMPIINSERT_IF_RANK_EQ_0      
      WRITE(*,*)'eqfndpsi/END: lr_,iter,rhonew,rhodes=', 
     +  lr_,iter,rhonew,rhodes
      !pause
CMPIINSERT_ENDIF_RANK

      return
      end
