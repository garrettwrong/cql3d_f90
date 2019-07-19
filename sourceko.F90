module sourceko_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use bcast_mod, only : ibcast
  use fle_mod, only : fle
  use lookup_mod, only : lookup
  use r8subs_mod, only : dcopy
  use r8subs_mod, only : dscal
  use soucrit_mod, only : soucrit
  use tdoutput_mod, only : tdoutput

  !---END USE

!
!

contains

      subroutine sourceko ! must be called for k=kelecg only
      use param_mod
      use cqlcomm_mod
      use r8subs_mod, only: lug, luf, dscal, dcopy
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save
!
!.......................................................................

!  Computes Marshall Rosenbluth's knock-on source
!    (which differs from Besedin-Pankratov expression.)
!    Source is obtained as the bounce-average of
!       n_e(s)*dn_r*v_r*(d(sigma)/d(gamma))*d(gamma)/d^3u,
!       where dn_r=runaway electrons at given energy gamma,
!       cross-section sigma on p. 240, of Heitler, Quantum Theory
!       of Radiation, 3rd edition,
!       s is poloidal distance, and d^3 is momentum-space-per-mass
!       volume element at given s.
!       Numerical method in Harvey et al, Phys. of Plasmas, 2000.
!
!.......................................................................
!      save
#ifdef __MPI
!MPI >>>
      include 'mpilib.h'
!MPI <<<
#endif

!      The ii0m1 method is an alternative translation from
!      pitch angle at given poloidal position to equatorial
!      plane pitch angle index (see below).
!      dimension ii0m1(jx,jfl,lz)

      character*8 ifirst
      data ifirst/"first"/
!
!     write(*,*) 'In SOURCEKO'
!
!     Should probably replace jfl with jx, etc (bh, 980501).
      !write(*,*)'sourceko jfl=', jfl
      !if (jfl.lt.jx) stop 'in sourceko'  ! YuP:  Why needed?

      if (n.lt.nonko .or. n.gt.noffko .or. &
             knockon.eq."disabled")  return

      k=kelecg
!
!     At the first call, set up quantities which need only be calculated
!       once: the p_par,p_perp intersection points between source
!       lines associated with each parallel velocity grid point and the
!       u-grid.  Also,  area elements, and source factors.
!       Momenta p_par and p_prp, are measured in units of mc.
!       Similarly, momenta area is in units of (mc)**2.

      if (ifirst.eq."first") then
         ifirst="notfirst"
         call bcast(ppars,zero,jx*jfl)
         call bcast(pprps,zero,jx*jfl)
!  pparea storage can be removed:   call bcast(pparea,zero,jx*jfl)
         call bcast(faci,zero,jx*jfl)

!     Set up parallel grid  for fl:
         call fle("setup",0)

         jflh=1

!        For avoiding divide errors (due to roundoff):
         abit=1.d-12
         onepabit=1.+abit
         onemabit=1.-abit

!    Determine index of u-contour intercepted by source line from each
!     parallel velocity such that gamma(j-1) satisfies
!     gamma(j-1)<= 0.5*(gamma_prime_1 +1),
!     but gamma(j) does not satisfy this relation.
         do jf=jflh,jfl
            gamap1=sqrt(xl(jf)*xl(jf)*cnorm2i+1.d0) !cnorm2i=0 when relativ.eq."disabled"
!BH180719            xmx1=5.e-1*sqrt(xl(jf)*xl(jf)+(2.*gamap1-2.e0)*cnorm2)
            xmx1=5.d-1*sqrt(xl(jf)*xl(jf)+(2.*gamap1-2.d0)*cnorm2)
            jmaxxl(jf)=luf(xmx1*onemabit,x,jx)
         enddo
         jmaxxl(jflh)=1

         do jf=jflh+1,jfl
            xpar1p=xl(jf)
!YuP            ppar1p=xpar1p*cnormi*onepabit !cnormi=0 when relativ.eq."disabled"
            ppar1p=xpar1p*onepabit/cnorm !YuP[07-2016] cannot have ppar1p or ppar1p2=0
            ppar1p2=ppar1p*ppar1p
            gam1p=sqrt(ppar1p2+1.d0) ! can be ~ 1.
            gam1pm1=gam1p-1.d0 ! can be ~ 0.
            gam1p2=gam1p*gam1p
            gam1p2m1=gam1p2-1.d0 ! can be ~0.
            alf1=2.*(gam1p-1.d0) ! can be ~0.
            if(ppar1p.eq.zero)then
             write(*,*)'sourceko: ppar1p=0=',ppar1p
             !pause
            endif
            if(alf1.eq.zero)then
             write(*,*)'sourceko: alf1=0=',alf1
             !pause
            endif
            aaa=(1.d0/alf1)-(1.d0/ppar1p2) !cannot have ppar1p2=0
            bbb=1.d0/ppar1p
            if(aaa.eq.zero)then
             write(*,*)'sourceko: aaa=0=',aaa ! see further: 1/(2.*aaa)
             !pause
            endif

            do j=2,jmaxxl(jf)-1
!YuP               p02=xsq(j)*cnorm2i*onepabit !cnorm2i=0 when relativ.eq."disabled"
               p02=(xsq(j)/cnorm2)*onepabit !YuP
               gam0=sqrt(p02+1.d0)
               gam0m1=gam0-1.d0
               gam1pmg=gam1p-gam0
               ccc=-p02/alf1

            if(bbb*bbb-4.d0*aaa*ccc.lt.zero)then
             write(*,*)'sourceko: bbb*bbb-4.d0*aaa*ccc=', &
              bbb*bbb-4.d0*aaa*ccc
             !pause
            endif
            if(p02-ppars(j,jf)*ppars(j,jf).lt.zero)then
             write(*,*)'sourceko: p02-ppars(j,jf)*ppars(j,jf)=', &
              p02-ppars(j,jf)*ppars(j,jf)
             !pause
            endif
            if(gam1p2m1*gam0m1*gam0m1*gam1pmg*gam1pmg.eq.zero)then
             write(*,*) &
              'sourceko: gam1p2m1*gam0m1*gam0m1*gam1pmg*gam1pmg=', &
              gam1p2m1*gam0m1*gam0m1*gam1pmg*gam1pmg
             !pause
            endif


               ppars(j,jf)=(-bbb+sqrt(bbb*bbb-4.d0*aaa*ccc))/(2.d0*aaa)
               pprps(j,jf)=sqrt(p02-ppars(j,jf)*ppars(j,jf))
!BH180719               pprps2=pprps(j,jf)*pprps(j,jf)
!BH180719               ppar1pmp=ppar1p-ppars(j,jf)
!BH180719               ppar1pm2=ppar1pmp*ppar1pmp
! Capital Sigma, Harvey etal, Checked BH180719
               faci(j,jf)=1.d0/(gam1p2m1*gam0m1*gam0m1*gam1pmg*gam1pmg) &
                    *(gam1pm1*gam1pm1*gam1p2 &
                    -gam0m1*gam1pmg*(2.*gam1p2+2.d0*gam1p-1.d0 &
                    -gam0m1*gam1pmg))

            enddo
         enddo



!     Setup equatorial plane pitch angle bin for
!     each source contribution.
!c     This is flux surface dependent, so must be recalculated
!c     if multiple flux surfaces are treated.
!c     This would be a large amount of storage for multiple
!c     flux surface runs.  (Could reduce storage by packing
!c     the ii0m1 into byte size words).
!c     It is more exact than the following method, but
!c     the replacement uses less storage.
!
!c     Starting value for lug() below.
!      i0=iyh/2
!
!         do l=1,lz
!            do jf=jflh+1,jfl-1
!               do j=1,jx
!c                 Positive pitch angles:
!                  sinth=pprps(j,jf)*cnorm*xi(j)
!                  i0=lug(sinth,sinz(1,l,lr_),iyh_(l_),i0)
!                  i0m1=i0-1
!                  if (i0m1.eq.0) then
!                     i0m1=1
!                     i0=2
!                  endif
!                  ii0m1(j,jf,l)=i0m1
!               enddo
!            enddo
!         enddo

!     The following method sets up a table of equatorial
!       plane pitch angles, corresponding to uniformly spaced
!       sin(theta) values for each poloidal position il,
!       and radial position ir.
!       This table is used later in the subroutine to find
!       equatorial pitch angle indices using simple divide
!       and multiply operators, and is much faster than
!       redoing the lug-searches.  (Differences (of +/- 1)
!       from the previous lug-method occur in typically 1-5
!       out of jx=350 points, for i0param=1000).

         call ibcast(i0tran,1,(i0param+1)*lz*setup0%lrz)
!000123BH         dsinth=0.5*pi/(i0param-1.) !Fix below increases accuracy abit
!BH170927:  But, this correction made a significant difference in the
!BH170927:  calculation of the wh80 test case, to the distribution.
!BH170927:  and the knockon rate. It apparently lead to the somewhat
!BH170927:  unphysical looking "ledge" on f near the runaway velocity.
!BH180416: Temporarily reverting this fix dsinth=1.0/(i0param-1.) !This is corrected value
!BH180416:         dsinth=0.5*pi/(i0param-1.)
         dsinth=1.0d0/(i0param-1)
!BH180416:       write(*,*)
!BH180416:       write(*,*) 'souceko: Temporary reversion of dsinth=. NEEDS'
!BH180416:       write(*,*) '         INVESTIGATION, 180416 !!!!'
!BH180416:       write(*,*)
         dsinthi=1./dsinth
         do ir=1,setup0%lrz
            do il=1,lz
               i0tran(1,il,ir)=1
               i0tran(i0param+1,il,ir)=iyh_(ir)
               i0=iyh_(ir)
               do ii=2,i0param
                  sinth=(ii-1)*dsinth
                  i0=lug(sinth,sinz(1:iyh_(ir),il, &
                   setup0%lrindx(ir)),iyh_(ir),i0)
                  i0m1=i0-1
                  if (i0m1.eq.0) then
                     i0m1=1
                  endif
                  i0tran(ii,il,ir)=i0m1
               enddo
            enddo
         enddo



!       second2=second()
!       write(*,*) 'Seconds for sourcko setup = ',second2-second1

      endif !On ifirst.eq."first"

!  end of setup for knock-on source integral.


!.......................................................................
!  Formation of knock-on source:
!.......................................................................

!      second3=second()


!  Constant factor:
!     tpir02=2*pi*r0**2 (cgs), r0=e**2/(mc**2), classical elec radius
      tpir02=2.*pi*(2.8179e-13)**2


!c     Source off below soffvte*vte (soffvte positive),
!                      abs(soffvte)*sqrt(E_c/E)*vte (soffvte negative).
!     E_c=2.*elecr (below, also in subroutine efield).
!     elecr is the Dreicer electric field (as in Kulsrud et al.),
!     converted to volts/cm.
!.......................................................................
!     call a routine to determine runaway critical velocity (used in
!       sourceko).
!.......................................................................

      call soucrit

!..................................................................
!     Note: the tauee and elecr definition below uses vth(),
!     vth is the thermal velocity =sqrt(T/m) (at t=0 defined in ainpla).
!     But, T==temp(k,lr) can be changed in profiles.f,
!     in case of iprote (or iproti) equal to "prbola-t" or "spline-t"
!..................................................................
      xvte3=0.
      if (isoucof.eq.0) then
         if (soffvte.gt.0.) then
            xvte3=soffvte*vthe(lr_)/vnorm
         else
            tauee(lr_)=vth(kelec,lr_)**3*fmass(kelec)**2 &
                 /(4.*pi*reden(kelec,lr_)*charge**4* &
                 gama(kelec,kelec))
            elecr(lr_)=300.*fmass(kelec)*vth(kelec,lr_)/ &
                 (2.*charge*tauee(lr_))
            if(elecfld(lr_).ne.zero) then
               xvte3=abs(soffvte)*sqrt(2.*elecr(lr_)/elecfld(lr_))* &
                    vthe(lr_)/vnorm
               !YuP[2018-09-19] Perhaps should use abs(elecfld(lr_)) in the above
               !to avoid sqrt(neg.value)
            endif
         endif
      else
         if (soffvte.gt.0.) then
            if (faccof.lt.0.7) faccof=0.7
!SC_did_not_like:        xvte3=amin1(soffvte,ucrit(kelec,lr_)*faccof)
!990131            xvte3=amin1(soffvte*vthe(lr_)/vnorm,ucrit(kelec,lr_)*faccof)
            xvte3=min(soffvte*vthe(lr_)/vnorm,ucrit(kelec,lr_)*faccof)
         else
            if (faccof.lt.0.5) faccof=0.5
            xvte3=ucrit(kelec,lr_)*faccof
         endif
      endif
      call lookup(xvte3,x,jx,wtu,wtl,lement)
      jvte3=lement

!  Distribution of primary knockon particles zeroed out (effectively)
!    below velocity soffpr*(above source cutoff velocity):
!    (default soffpr=0.0)

      xsoffpr=soffpr*xvte3
      ! XXX    YuP: was xl(jflh), where jflh=(jfl+1)/2, and I don't know why.
      call lookup(xsoffpr,xl,jfl,wtu,wtl,lement)
      jsoffpr=min(lement,jfl)


!  Obtain average source for positive and negative v_parallel:

      denfl(l_)=0.0  !FSA density, for diag purposes
      do l=1,lz  !Over poloidal angle

         call fle("calc",l)  !Reduced distn at each pol angle

         den_of_s(l)=0.  !Dens at each pol angle, save for diag purposes
         do jf=1,jfl-1
            den_of_s(l)=den_of_s(l)+dxl(jf)*(fl1(jf)+fl2(jf))
         enddo
         denfl(l_)=denfl(l_) &
                +(dz(l,lr_)/bbpsi(l,lr_))*den_of_s(l)/zmaxpsi(lr_)

         flmax=0.
         do jf=1,jfl-1
!990131            flmax=amax1(flmax,fl1(jf))
!990131            flmax=amax1(flmax,fl2(jf))
            flmax=max(flmax,fl1(jf))
            flmax=max(flmax,fl2(jf))
         enddo
         flmin=em100*flmax
!     do jf=jflh-jsoffpr,jflh+jsoffpr
         do jf=1,jsoffpr-1
            fl1(jf)=flmin
            fl2(jf)=flmin
         enddo




      do jf=jflh+1,jfl-1
!     Skip calc if fl(jf) near minimum value:
            if ((fl1(jf)+fl2(jf)).lt.em100*flmax) go to 99

            xpar1p=xl(jf)
            ppar1p=xpar1p*vnorm/clight*onepabit
            ppar1p2=ppar1p*ppar1p
            gam1p=sqrt(ppar1p2+1.)

            cnst0=tpir02*den_of_s(l)*dxl(jf)*ppar1p*vnorm/ &
                 (gam1p*cnorm)
            cnst1=cnst0*fl1(jf)
            cnst2=cnst0*fl2(jf)



      do j=jvte3,jmaxxl(jf)-2
         wc=x(j)/cnorm
         wc2=wc*wc

         src_mr1=cnst1*faci(j,jf)*gammi(j)*x(j)*dx(j)
         src_mr2=cnst2*faci(j,jf)*gammi(j)*x(j)*dx(j)

!         i0m1=ii0m1(j,jf,l)

         sinth=pprps(j,jf)*cnorm*xi(j)
         isinth=sinth*dsinthi+1.5
         i0m1=i0tran(isinth,l,l_)
!         if(i0m1.ne.i0m1o) then
!            write(*,*) 'l_,l,jf,j,i0m1o,i0m1  ', l_,l,jf,j,i0m1o,i0m1
!         endif

!BH091031         axl=dtau(i0m1,l,lr_)/tau(i0m1,lr_)
!BH091031         if(l.ne.lz .and. lmax(i0m1,lr_).eq.l) then
!BH091031            axl=axl+dtau(i0m1,l+1,lr_)/tau(i0m1,lr_)
!BH091031         endif
         if (eqsym.ne."none") then !i.e. up-down symm
            !if not bounce interval
            if(l.eq.lz .or. l.ne.lmax(i0m1,lr_)) then
               axl=dtau(i0m1,l,lr_)/tau(i0m1,lr_)
            else                !bounce interval: additional contribution
               axl=(dtau(i0m1,l,lr_)+dtau(i0m1,l+1,lr_))/tau(i0m1,lr_)
            endif
         else  !eqsym="none"
            if (l.lt.lz_bmax(lr_) .and. l.eq.lmax(i0m1,lr_))then
               !trapped, with tips between l and l+1 (above midplane)
               axl=(dtau(i0m1,l,lr_)+dtau(i0m1,l+1,lr_))/tau(i0m1,lr_)
               !-YuP  Note: dtau(i,l+1,lr_)=0
            elseif (l.gt.lz_bmax(lr_) .and. l.eq.lmax(i0m1+iyh,lr_))then
               !trapped, with tips between l and l-1 (below midplane)
               axl=(dtau(i0m1,l,lr_)+dtau(i0m1,l-1,lr_))/tau(i0m1,lr_) !NB:l-1
               !-YuP  Note: dtau(i,l-1,lr_)=0
            else
               axl=dtau(i0m1,l,lr_)/tau(i0m1,lr_)
               !passing (i<itl), or trapped but with tips at other l;
               !also, at l=lz_bmax, includes last trapped particle i=itl
               !(for such particle, lmax(itl)=lz_bmax; see micxinil)
            endif
         endif
         dtaudt=axl

         i0m1m=iy+1-i0m1

         source(i0m1,j,k,indxlr_)=source(i0m1,j,k,indxlr_)+ &
                     dtaudt*src_mr1*cosz(i0m1,l,lr_)/(coss(i0m1,l_) &
                     *cynt2(i0m1,l_)*cint2(j)*bbpsi(l,lr_))
         source(i0m1m,j,k,indxlr_)=source(i0m1m,j,k,indxlr_)+ &
                     dtaudt*src_mr2*cosz(i0m1m,l,lr_)/(coss(i0m1m,l_) &
                     *cynt2(i0m1m,l_)*cint2(j)*bbpsi(l,lr_))


      enddo  !On j=jvte3,jmaxxl(jf)-2
 99   continue
      enddo  !On jf=jflh+1,jfl-1
      enddo  !On l=1,lz

!     Symmetrize trapped particles
      do  j=jvte3,jx
         do  i=itl,iyh
            ii=iy+1-i
            source(i,j,k,indxlr_)=0.5*(source(i,j,k,indxlr_) &
                +source(ii,j,k,indxlr_)) ! for k=kelecg
            source(ii,j,k,indxlr_)=source(i,j,k,indxlr_)
         enddo
      enddo

!..................................................................
!     For primary distribution formed from the pitch angle averaged
!       reduced distribution, reverse the direction of source if
!       the electric field is positive.
!..................................................................

      if (elecfld(lr_).gt.0.) then
         do  j=jvte3,jx
            do  i=1,itl-1
               ii=iy+1-i
               source(ii,j,k,indxlr_)=source(i,j,k,indxlr_)
               source(i,j,k,indxlr_)=0.0
            enddo
         enddo
      endif

!..................................................................
!     Compute source rate in (particles/sec/cc)
!..................................................................

          s1=0.
          do j=jvte3,jx
          do i=1,iy
             s1=s1+ &
                source(i,j,k,indxlr_)*cynt2(i,l_)*cint2(j)*vptb(i,lr_)
          enddo
          enddo
          s1=s1/zmaxpsi(lr_)
          srckotot(lr_)=s1

!..................................................................
!     Compute pitch angle averaged source rate
!     in (particles/sec/cc) divided by du=vnorm*dx(j)
!     (Can examine this with CDBX).
!..................................................................

          call bcast(tam30,zero,jx)
          do i=1,iy
             do j=1,jx
                tam30(j)=tam30(j)+source(i,j,k,indxlr_)* &
                     cynt2(i,l_)*vptb(i,lr_)* &
                     cint2(j)/(vnorm*dx(j)*zmaxpsi(lr_))
             enddo
          enddo

!          second4=second()

!       write(*,*) 'Seconds for sourcko calc = ',second4-second3


!..................................................................
!     If knockon.eq."fpld_dsk":
!        write source out to disk file "fpld_disk", and stop,
!     If knockon.eq."fpld_ds1":
!        write source*dtr+(previous step f(since no preloading in it)
!        to disk file "fpld_dsk1", and previous step f to disk file
!        "fpld_dsk2",     and stop.
!..................................................................

      if (n.eq.nstop) then
      if (knockon.eq."fpld_dsk") then
         call dcopy(iyjx2,source(0:iy+1,0:jx+1,k,l_),1, &
                           temp1(0:iy+1,0:jx+1),1)
         call dscal(iyjx2,dtr,temp1(0:iy+1,0:jx+1),1)
         iunit=20
         open (unit=iunit,file='fpld_dsk',status='unknown')
         !do 80 k=1,1 !  k=kelecg
            write (iunit,1000) ((temp1(i,j),i=1,iy),j=1,jx)
 80      continue
         close(unit=iunit)
         call tdoutput(2)
!990307         call geglxx(0)
         call pgend
         stop 'Wrote disk file from subroutine souceko: fpld_dsk'
      endif
      if (knockon.eq."fpld_ds1") then
         call dcopy(iyjx2,source(0:iy+1,0:jx+1,k,l_),1, &
                           temp1(0:iy+1,0:jx+1),1)
         call dscal(iyjx2,dtr,temp1(0:iy+1,0:jx+1),1)
         do j=1,jx
            do i=1,iy
               temp1(i,j)=temp1(i,j)+f_(i,j,k,1)
            enddo
         enddo
         iunit=20
         open (unit=iunit,file='fpld_dsk1',status='unknown')
         !do 81 k=1,1  !  k=kelecg
            write (iunit,1000) ((temp1(i,j),i=1,iy),j=1,jx)
 81      continue
         close(unit=iunit)
         open (unit=iunit,file='fpld_dsk2',status='unknown')
         !do 82 k=1,1  !  k=kelecg
            write (iunit,1000) ((f_(i,j,k,1),i=1,iy),j=1,jx)
 82      continue
         close(unit=iunit)
         call tdoutput(2)
!990307         call geglxx(0)
         call pgend
         stop 'Wrote disk files from subroutine souceko: fpld_dskx'
      endif
      endif
 1000 format(5(1pe16.7))




!..................................................................
!     As passed to subroutine souplt:
!..................................................................
      xlncur(1,lr_)=s1*zmaxpsi(lr_)

 999  return
      end subroutine sourceko


end module sourceko_mod
