! YuP: should be removed from list of sources.
!BH171014: Commenting, since apparently not called, and gfortran
!BH171014: is having a problem with it.

!$$$!
!$$$!
!$$$!     ONETWO DIVERGENCE
!$$$      subroutine frhexdrv(mb,sgxn,sgxnmi,hxfrac)
!$$$      implicit integer (i-n), real(c_double) (a-h,o-z)
!$$$!     write statements commented out
!$$$!
!$$$!---hexdrv is a driver for routine hexnb.
!$$$!---hexnb calculates the total reaction rate,rerate(per second), and
!$$$!---the total mean free path rmfp.
!$$$!---hexdrv gets the cross section,sgxn, and inverse of the
!$$$!---pseudocross-section, sgxnmi. hexnb takes excited states of the
!$$$!---neutral beam into account if iexcit=2. if iexcit=1 hexnb should retu
!$$$!---essentially the same cross sections as crsecs does. note
!$$$!---however that the cross section parameterization is different
!$$$!---inhexnb and crsecs.
!$$$!     ONETWO DIVERGENCE
!$$$      include 'frcomm.h77'
!$$$!
!$$$!
!$$$!     ONETWO DIVERGENCE
!$$$      dimension sgxn(kz,ke,kb),sgxnmi(ke,kb),hxfrac(ke,kb)
!$$$!
!$$$!
!$$$!---setup variables for hexnb
!$$$!
!$$$      iexcm1=iexcit-1
!$$$      istart=1
!$$$      if(nouthx.ne.0)nouthx=37
!$$$!     if(iexcm1.eq.0)write(ncrt,10)
!$$$!     if(iexcm1.eq.1)write(ncrt,15)mstate+1
!$$$!     10 format(2x,'entering hexnb,no beam excitation ')
!$$$!     15 format(2x,'entering hexnb, with beam excitation accounted for',
!$$$!     &/,'max principal assumed quantum number is',i5,/,
!$$$!     &'lorentz ionization limit not implemented at this time')
!$$$!
!$$$!---separate primary and impurity ions for easier reading
!$$$!---note:nion=nimp+nprim
!$$$!
!$$$      do 100 j=1,nion
!$$$        if(j.gt.nprim)go to 120
!$$$        atwpm(j)=atw(j)
!$$$        go to 100
!$$$ 120    i=j-nprim
!$$$        atwim(i)=atw(j)
!$$$ 100  continue
!$$$!
!$$$!---get atomic number of impurities
!$$$!
!$$$      do 180 i=1,nimp
!$$$          if(trim(namei(i)).eq.'he' .or.
!$$$     +       trim(namei(i)).eq.'HE' .or.
!$$$     +       trim(namei(i)).eq.'He')  iz(i) =  2
!$$$
!$$$          if(trim(namei(i)).eq.'b' .or.
!$$$     +       trim(namei(i)).eq.'B' )  iz(i) =  5  ! YuP added [2015]
!$$$
!$$$          if(trim(namei(i)).eq.'c' .or.
!$$$     +       trim(namei(i)).eq.'C' )  iz(i) =  6
!$$$
!$$$          if(trim(namei(i)).eq.'o' .or.
!$$$     +       trim(namei(i)).eq.'O' )  iz(i) =  8
!$$$
!$$$          if(trim(namei(i)).eq.'si' .or.
!$$$     +       trim(namei(i)).eq.'SI' .or.
!$$$     +       trim(namei(i)).eq.'Si')  iz(i) = 14
!$$$
!$$$          if(trim(namei(i)).eq.'ar' .or.
!$$$     +       trim(namei(i)).eq.'AR' .or.
!$$$     +       trim(namei(i)).eq.'Ar')  iz(i) = 18
!$$$
!$$$          if(trim(namei(i)).eq.'cr' .or.
!$$$     +       trim(namei(i)).eq.'CR' .or.
!$$$     +       trim(namei(i)).eq.'Cr')  iz(i) = 24
!$$$
!$$$          if(trim(namei(i)).eq.'fe' .or.
!$$$     +       trim(namei(i)).eq.'FE' .or.
!$$$     +       trim(namei(i)).eq.'Fe')  iz(i) = 26
!$$$
!$$$          if(trim(namei(i)).eq.'ni' .or.
!$$$     +       trim(namei(i)).eq.'NI' .or.
!$$$     +       trim(namei(i)).eq.'Ni')  iz(i) = 28
!$$$
!$$$          if(trim(namei(i)).eq.'kr' .or.
!$$$     +       trim(namei(i)).eq.'KR' .or.
!$$$     +       trim(namei(i)).eq.'Kr')  iz(i) = 36
!$$$
!$$$          if(trim(namei(i)).eq.'mo' .or.
!$$$     +       trim(namei(i)).eq.'MO' .or.
!$$$     +       trim(namei(i)).eq.'Mo')  iz(i) = 42
!$$$
!$$$          if(trim(namei(i)).eq.'w' .or.
!$$$     +       trim(namei(i)).eq.'W' )  iz(i) = 74
!$$$ 180  continue
!$$$!
!$$$!---create diagnostic output file hexnbout
!$$$!---hexnbout contains detailed cross section info
!$$$!
!$$$      if(nouthx.gt.0)
!$$$     &  open(unit=nouthx,file='hexnbout',status='new')
!$$$!
!$$$!---loop over spatial zones,number of beams and beam components
!$$$!
!$$$      do 200 i=1,mfm1
!$$$        teev=1.e3*zte(i)
!$$$        tiev=1.e3*zti(i)
!$$$        do 220 j=1,nion
!$$$          if(j.gt.nprim)go to 230
!$$$          znipm(j)=zni(i,j)
!$$$          go to 220
!$$$ 230      k=j-nprim
!$$$          zniim(k)=zni(i,j)
!$$$ 220    continue
!$$$        do 240 ib=1,mb
!$$$          do 250 j=1,3
!$$$            er0=1.e3*ebkev(ib)/(j*atw(ibion))
!$$$!---no lorentz ionization limit for now
!$$$            bperp=0.0
!$$$            call hexnb(istart,iexcm1,ilorent,mstate,ncont,
!$$$     &        er0,teev,tiev,nprim,atwpm,znipm,nimp,
!$$$     &        iz,atwim,izstrp,zniim,bperp,kdene,kdeni,
!$$$     &        kdenz,ksvi,ksvz,ksve,krad,ngh,ngl,nouthx,
!$$$     &        ncorin,rerate,rmfp,hexfrac,ihxerr)
!$$$!     if(ihxerr.ne.0)write(ncrt,20)ihxerr
!$$$            if(ihxerr.eq.1)then
!$$$!     write(ncrt,21)
!$$$!     write(nout,21)
!$$$              stop 'frhexdrv: problem'
!$$$            endif
!$$$!     20 format(2x,'hexnb error = ',i5)
!$$$!     21 format(2x,'execution terminated - file coronb absent')
!$$$            sgxn(i,j,ib)=1./rmfp
!$$$            if(i.eq.1)hxfrac(j,ib)=hexfrac
!$$$ 250      continue
!$$$ 240    continue
!$$$ 200  continue
!$$$!---
!$$$!---get inverse of pseudo cross section
!$$$!---
!$$$      do 400 ib=1,mb
!$$$        if((ib.gt.1).and.(ebkev(ib).eq.ebkev(1)))go to 450
!$$$        do 420 j=1,3
!$$$          sgxnm=sgxn(1,j,ib)
!$$$          do 430 i=2,mfm1
!$$$c990131 430      sqxnm=amax1(sgxnm,sgxn(i,j,ib))
!$$$ 430      sqxnm=max(sgxnm,sgxn(i,j,ib))
!$$$ 420    sgxnmi(j,ib)=1./sgxnm
!$$$        go to 400
!$$$ 450    do 460 j=1,3
!$$$ 460    sgxnmi(j,ib)=sgxnmi(j,1)
!$$$ 400  continue
!$$$      write(*,*) 'sgxnmi',sgxnmi(1,1),sgxnmi(2,1),sgxnmi(3,1)
!$$$      write(*,*) 'sgxn',sgxn(1,1,1),sgxn(2,1,1)
!$$$      write(ncrt,16)
!$$$ 16   format(2x,'leaving hexnb ')
!$$$      return
!$$$      end
