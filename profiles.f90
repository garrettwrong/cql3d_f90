module profiles_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use cfpmodbe_mod, only : cfpmodbe
  use profaxis_mod, only : profaxis
  use tdinterp_mod, only : tdinterp
  use tdxin13d_mod, only : tdxin13d

  !---END USE

!
!

contains

      subroutine profiles
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save
!
!.......................................................................
!     Obtain time-dependent
!     "parabolic" and "spline" profiles, for time-varying quantities.
!.......................................................................

!MPIINSERT_INCLUDE

      real(c_double):: tmpt(njene)  !Temporary array
      dimension ztr(lrza)   !For tdoutput-like printout
      dimension rban_vth(setup0%lrzmax) ! for print-out
      character*8 ztext

      if (nbctime.le.0) return

      itab(1)=1
      itab(2)=0
      itab(3)=0

      if (tmdmeth.eq."method1") then
         itme=0
         do jtm=1,nbctime
            if (timet.ge.bctime(jtm)) itme=jtm
         enddo
         itme1=itme+1
      endif

!     Temperatures



      if (iprote.eq."prbola-t") then

         do 5 k=1,ntotal
            if (tmdmeth.eq."method1") then

               if (itme.eq.0) then
                  temp(k,0)=tempc(1,k)
                  temp(k,1)=tempb(1,k)
               elseif (itme.lt.nbctime) then
                  temp(k,0)=tempc(itme,k)+(tempc(itme1,k)-tempc(itme,k)) &
                       /(bctime(itme1)-bctime(itme)) &
                       *(timet-bctime(itme))
                  temp(k,1)=tempb(itme,k)+(tempb(itme1,k)-tempb(itme,k)) &
                       /(bctime(itme1)-bctime(itme)) &
                       *(timet-bctime(itme))
               else
                  temp(k,0)=tempc(nbctime,k)
                  temp(k,1)=tempb(nbctime,k)
               endif

               !YuP: called later: call tdxin13d(temp,rya(1:ubound(rya,1)),setup0%lrzmax,ntotala,k,npwr(k),mpwr(k))

            elseif (tmdmeth.eq."method2")  then

               temp(k,0)=tdprof(timet,tempc(:,k),bctime)
               temp(k,1)=tdprof(timet,tempb(:,k),bctime)

            endif

            call tdxin13d(temp,rya,setup0%lrzmax,ntotala,k,npwr(k),mpwr(k))

 5       continue
      endif ! iprote=prbola-t


      if (iprote.eq."spline-t") then

!        Electron temperature

         if (itme.eq.0) then
            do l=1,njene
               tmpt(l)=tein_t(l,1)
            enddo
         elseif (itme.lt.nbctime) then
            do l=1,njene
               tmpt(l)=tein_t(l,itme)+(tein_t(l,itme1)-tein_t(l,itme)) &
                     /(bctime(itme1)-bctime(itme))*(timet-bctime(itme))
            enddo
         else
            do l=1,njene
               tmpt(l)=tein_t(l,nbctime)
            enddo
         endif

         do 16  k=1,ntotal
            if(bnumb(k).eq.-1.)  then
               call tdinterp("zero","linear",ryain,tmpt,njene,rya(1:ubound(rya,1)), &
                    tr(1),setup0%lrzmax)
               tr(0)=tmpt(1)
               do 19  ll=0,setup0%lrzmax
                  temp(k,ll)=tr(ll)
                  if(temp(k,ll).le.0.001)then
                   write(*,*)'temp<0.001.  tmpt(1:njene)=',tmpt
                   write(*,*)'temp<0.001.  temp(k,ll)=',k,ll,temp(k,ll)
                   !pause    ! temp<0.001
                  endif
 19            continue ! ll
            endif
 16      continue ! k
      endif ! iprote=spline-t

      if (iproti.eq."spline-t") then
!        Ion temperatures
         if (itme.eq.0) then
            do l=1,njene
               tmpt(l)=tiin_t(l,1)
            enddo
         elseif (itme.lt.nbctime) then
         do l=1,njene
            tmpt(l)=tiin_t(l,itme)+(tiin_t(l,itme1)-tiin_t(l,itme)) &
                     /(bctime(itme1)-bctime(itme))*(timet-bctime(itme))
         enddo
         else
            do l=1,njene
               tmpt(l)=tiin_t(l,nbctime)
            enddo
         endif

         do 26  k=1,ntotal
            if(bnumb(k).ne.-1.)  then
               call tdinterp("zero","linear",ryain,tmpt,njene,rya(1:ubound(rya,1)), &
                    tr(1),setup0%lrzmax)
               tr(0)=tmpt(1)
               do 29  ll=0,setup0%lrzmax
                  temp(k,ll)=tr(ll)
 29            continue
            endif
 26      continue
      endif  ! (iproti.eq."spline-t")

      do l=1,setup0%lrzmax
      if (tempc(l,kelec).le.zero .and. tein_t(l,1).le.zero) then
         WRITE(*,*) "Time-dependent Te profile input problem at l=",l
         STOP
      endif
      enddo

!     Assume any time dep in tempc is in first ion species.
      do l=1,setup0%lrzmax
      if (tempc(l,kionn).le.zero .and. tiin_t(l,1).le.zero) then
         WRITE(*,*) "Time-dependent Ti profile input problem at l=",l
         STOP
      endif
      enddo

!     Renormalize temperatures using tescal/tiscal
      do k=1,ntotal
         if (bnumb(k).eq.-1.) then
            do l=0,setup0%lrzmax
               temp(k,l)=tescal*temp(k,l)
            enddo
         else
            do l=0,setup0%lrzmax
               temp(k,l)=tiscal*temp(k,l)
            enddo
         endif
         do l=1,setup0%lrzmax
            if (temp(k,l).le.zero) then
              WRITE(*,*) "profiles.f: temp(k,l)<0 at k,l=",k,l
              STOP
            endif
         enddo ! l
      enddo ! k


!..................................................................
!     YuP[2018-01-05] Added resetting of vth() array in profiles.f.
!     vth is the thermal velocity =sqrt(T/m) (at t=0 defined in ainpla).
!     But, T==temp(k,lr) was changed above,
!     in case of iprote (or iproti) equal to "prbola-t" or "spline-t".
!     Then, reset the values of vth() [used for plots, and also
!     in lossegy, efield, restvty, sourceko, tdnpa, tdtrdfus,
!     tdoutput, and vlf*,vlh* subroutines].
      do k=1,ntotal
         do l=1,setup0%lrzmax
           vth(k,l)=(temp(k,l)*ergtkev/fmass(k))**.5
           if (k .eq. kelec) vthe(l)=vth(kelec,l)
         enddo
      enddo

      !YuP [2018-09-18] Added warning, similar to micxinit
      !Check that there are sufficient number of v-grid points
      !over thermal part of the distribution (at least 3 points).
      ! If not, print warning.
      if(n.gt.0)then ! skip it at n=0 (x(j) not defined yet)
!MPIINSERT_IF_RANK_EQ_0
      WRITE(*,'(a)')"==============================================="
      WRITE(*,*)'profiles.f: time step n, timet=',n,timet
!MPIINSERT_ENDIF_RANK
      do k=1,ngen
      do l=1,setup0%lrzmax
         j_thermal=1 ! to initialize
         do j=jx,1,-1 ! scan backward
           !write(*,*) j, vth(k,l)/vnorm, x(j)
           if(vth(k,l) .le. x(j)*vnorm)then
            j_thermal=j !it is such j that v(j_thermal) is just below vth()
           endif
         enddo
!MPIINSERT_IF_RANK_EQ_0
         WRITE(*,'(a,3i4,3f16.11)') &
         "profiles: k,lr, j_thermal, x(j_thermal), vth/vnorm, temp =", &
               k, l, j_thermal, x(j_thermal), vth(k,l)/vnorm, temp(k,l)
         if(j_thermal.lt.3)then
           WRITE(*,'(a)')" WARNING(profiles.f): V-grid is too coarse."
           WRITE(*,'(a)')" Thermal part of distrib. is covered by only"
           WRITE(*,'(a,i5,a)')"    j=", j_thermal," points."
           WRITE(*,'(a)')" The solution may become unstable."
           WRITE(*,'(a)')" Consider increasing jx or setting xfac<1."
           !pause !-------
         endif
!MPIINSERT_ENDIF_RANK
      enddo ! l=1,setup0%lrzmax
      enddo ! k=1,ngen
      endif ! n>0

!     The energy() array, initially set as 1.5*temp() [see ainpla.f]
!     is re-defined at each time step in subr. diaggnde
!     as the m*v^2/2 average over distr.function in vel. space.
!     BUT, it is done for k=1:ngen species only!
!     So, below, we re-define the energy() array for the Maxwellian
!     species only.
      do k=ngen+1,ntotal
         rstmss=fmass(k)*clite2/ergtkev
         do l=1,setup0%lrzmax
          thta=rstmss/temp(k,l)
          if (thta.gt.100. .or. relativ.eq."disabled") then
            energy(k,l)=1.5*temp(k,l)
          else
            call cfpmodbe(thta,bk1,bk2)
            energy(k,l)=rstmss*(bk1/bk2-1.+3./thta)
          endif
         enddo
      enddo
!     Another thought: Change the definition of vth() for k=1:ngen,
!     i.e. for the general species,
!     to be based on (2/3)*energy() rather than on temp() array.
!     [Not implemented for the now]
!..................................................................


!     Densities and Zeff  (cf. subroutine tdxinitl)

      if (iprone.eq."prbola-t") then

      do 11 k=1,ntotal

         if (iprozeff.ne."disabled" .and. &
                 (k.ne.kelecg .and. k.ne.kelecm)) go to 11

         if (redenc(1,k).ne.zero) then

            if (tmdmeth.eq."method1") then

               if (itme.eq.0) then
                  reden(k,0)=redenc(1,k)
                  reden(k,1)=redenb(1,k)
               elseif (itme.lt.nbctime) then
                  reden(k,0)=redenc(itme,k)+ &
                       (redenc(itme1,k)-redenc(itme,k)) &
                       /(bctime(itme1)-bctime(itme)) &
                       *(timet-bctime(itme))
                  reden(k,1)=redenb(itme,k)+ &
                       (redenb(itme1,k)-redenb(itme,k)) &
                       /(bctime(itme1)-bctime(itme)) &
                       *(timet-bctime(itme))
               else
                  reden(k,0)=redenc(nbctime,k)
                  reden(k,1)=redenb(nbctime,k)
               endif

            elseif (tmdmeth.eq."method2") then

               reden(k,0)=tdprof(timet,redenc(:,k),bctime)
               reden(k,1)=tdprof(timet,redenb(:,k),bctime)

            endif


            call tdxin13d(reden,rya,setup0%lrzmax,ntotala,k,npwr(0),mpwr(0))

         endif

 11   continue

      endif  !On iprone.eq.prbola-t


      if (iprone.eq."spline-t") then  !endif at line 245


      do  k=1,ntotal
         if (iprozeff.ne."disabled" .and. &
                 (k.ne.kelecg .and. k.ne.kelecm)) go to 12
         if (enein_t(1,k,1).ne.zero) then

            if (itme.eq.0) then
               do l=1,njene
                  tmpt(l)=enein_t(l,k,1)
               enddo
            elseif (itme.lt.nbctime) then
            do l=1,njene
               tmpt(l)=enein_t(l,k,itme)+ &
                    (enein_t(l,k,itme1)-enein_t(l,k,itme)) &
                    /(bctime(itme1)-bctime(itme))*(timet-bctime(itme))
            enddo
            else
            do l=1,njene
               tmpt(l)=enein_t(l,k,nbctime)
            enddo
            endif

            call tdinterp("zero","linear",ryain,tmpt,njene, &
                 rya(1:ubound(rya,1)),tr(1),setup0%lrzmax)
            tr(0)=tmpt(1)
            do 13 ll=0,setup0%lrzmax
               reden(k,ll)=tr(ll)
 13         continue
         else
            do 9  ll=0,setup0%lrzmax
               reden(k,ll)=tr(ll)/abs(bnumb(k))
 9          continue
         endif
 12      continue
      enddo  !On k

      endif  !On iprone.eq.spline-t

!     Finish up with zeff and ions if iprozeff.ne."disabled"
      if (iprozeff.ne."disabled") then  !endif at line 406
         if (iprozeff.eq."prbola-t") then

            if (zeffc(1).ne.zero) then

               if (tmdmeth.eq."method1") then

                  if (itme.eq.0) then
                     zeffin(0)=zeffc(1)
                     zeffin(1)=zeffb(1)
                  elseif (itme.lt.nbctime) then
                     zeffin(0)=zeffc(itme)+(zeffc(itme1)-zeffc(itme)) &
                          /(bctime(itme1)-bctime(itme)) &
                          *(timet-bctime(itme))
                     zeffin(1)=zeffb(itme)+(zeffb(itme1)-zeffb(itme)) &
                          /(bctime(itme1)-bctime(itme)) &
                          *(timet-bctime(itme))
                  else
                     zeffin(0)=zeffc(nbctime)
                     zeffin(1)=zeffb(nbctime)
                  endif

               elseif (tmdmeth.eq."method2") then

                  zeffin(0)=tdprof(timet,zeffc,bctime)
                  zeffin(1)=tdprof(timet,zeffb,bctime)

               endif  !On tmdmeth

            endif  !On zeffc(1)

            dratio=zeffin(1)/zeffin(0)
            do ll=1,setup0%lrzmax
               call profaxis(rn,npwrzeff,mpwrzeff,dratio,rya(ll))
               zeff(ll)=zeffin(0)*rn
            enddo

         elseif (iprozeff.eq."spline-t") then

               if (itme.eq.0) then
                  do l=1,njene
                     tmpt(l)=zeffin_t(l,1)
                  enddo
               elseif (itme.lt.nbctime) then
                  do l=1,njene
                     tmpt(l)=zeffin_t(l,itme) &
                          +(zeffin_t(l,itme1)-zeffin_t(l,itme)) &
                          /(bctime(itme1)-bctime(itme)) &
                          *(timet-bctime(itme))
                  enddo
               else
                  do l=1,njene
                     tmpt(l)=zeffin_t(l,nbctime)
                  enddo
               endif
               call tdinterp("zero","linear",ryain,tmpt,njene,rya(1:ubound(rya,1)), &
                    tr(1),setup0%lrzmax)
               do  ll=1,setup0%lrzmax
                  zeff(ll)=tr(ll)
               enddo
         endif

!     Scale zeff
         do ll=1,setup0%lrzmax
            zeff(ll)=zeffscal*zeff(ll)
         enddo

!     Check that range of bnumb for Maxl species brackets zeff
         fmaxx=0.
         fminn=0.
         do k=1,nionm
            fmaxx=max(fmaxx,bnumb(kionm(k)))
            fminn=min(fminn,bnumb(kionm(k)))
         enddo
         do 121 ll=1,setup0%lrzmax
            if(zeff(ll).gt.fmaxx .or. zeff(ll).lt.fminn) then
               WRITE(*,*)'profiles.f: ', &
                    'Adjust bnumb(kion) for compatibility with zeff'
               stop
            endif
 121     continue


!     Check number of ion Maxwl species with different bnumb.
!     (Need at least two in order to fit zeff. Check ainsetva.f.)
        if (nionm.lt.2) stop 'profiles: ion species problem'
        ndif_bnumb=1
        do k=2,nionm
           if (abs(bnumb(kionm(k))/bnumb(kionm(1))-1.).gt.0.01) &
                ndif_bnumb=ndif_bnumb+1
        enddo

!BH120627:  This appears to be already done above:
!BH120627:c    Interpolate input ion densities onto rya grid
!BH120627:     call tdxin13d(reden,rya,setup0%lrzmax,ntotala,k,npwr(0),mpwr(0))



!     Save ratios of densities for Maxwl ions with same charge in reden
!     to be used in following ion density calc from zeff.
!     (First ndif_bnumb ions have same change.  These are ion species
!     at the head of the ion species list.)

      nsame_bnumb=nionm-ndif_bnumb+1  !i.e., =1 if all diff bnumb ions
      if (nsame_bnumb.gt.1) then  !i.e., 2 or more ions have same bnumb.
                                  !Densities will contain density ratios.
                                  !Equal-bnumb() species at beginning of
                                  !Maxl ions.
!        Renormalizing the density ratios as fractions for equal-bnumb:
            do ll=1,setup0%lrzmax
               dsum=0.
               do k=1,nsame_bnumb
                  dsum=dsum+reden(kionm(k),ll)
               enddo
               do k=1,nsame_bnumb
                  reden(kionm(k),ll)=reden(kionm(k),ll)/dsum
               enddo
               reden(kionm(nionm),ll)=one
            enddo
      endif  ! on nsame_bnumb.gt.1

!     Set rest of ion densities to 1.
      do kk=nsame_bnumb,nionm
         k=kionm(kk)
         do l=0,setup0%lrzmax
            reden(k,l)=one
         enddo
      enddo

!     NOTE: reden(k, ) on rhs of following reden(,) formulas is
!     density ratio, .le.1. for diff species with same bnumb.
      do 14 k=kionm(1),kionm(nionm)
!     write(*,*)'tdxinitl, do 14 loop, k= ',k
!        For k pointing to equal-bnumb species (1 or more)
         if (k.le.(kionm(1)+nsame_bnumb-1)) then
            k1=k
            k2=kionm(nionm)
!        For k beyond equal-bnumb species
         else
            k1=kionm(nionm)
            k2=k-1
         endif
         reden(k,0)=reden(kelec,0)*reden(k,0)*(zeff(1) &
               -bnumb(k2))/(bnumb(k1)-bnumb(k2))/bnumb(k1)
         !YuP: From the above, reden can be a small negative value,
         !because of a rounding error. Add lower limit =0.d0
         reden(k,0)=max(reden(k,0),zero) !YuP[2018-09-18] added
         do 142 ll=1,setup0%lrzmax
               reden(k,ll)=reden(kelec,ll)*reden(k,ll)*(zeff(ll) &
                    -bnumb(k2))/(bnumb(k1)-bnumb(k2))/bnumb(k1)
            !YuP: From the above, reden can be a small negative value,
            !because of a rounding error. Add lower limit =0.d0
            reden(k,ll)=max(reden(k,ll),zero) !YuP[2018-09-18] added
 142     continue
         !write(*,*)'profiles:k1,k2=',k1,k2
         !write(*,*)'zeff=',zeff
         !write(*,*)'profiles:k,max(reden)=',k,MAXVAL(reden(k,:))
 14   continue



!     Copy Maxwellian ion densities to any ion general species
!     in order of the indexes.

         do k=1,niong
            if (nionm.ge.niong) then
               do l=0,setup0%lrzmax
                  reden(kiong(k),l)=reden(kionm(k),l)
               enddo
            endif
         enddo


      endif  ! on iprozeff.ne."disabled"

!     Renormalize densities using enescal
      do k=1,ntotal
         do l=0,setup0%lrzmax
            reden(k,l)=enescal*reden(k,l)
         enddo
      enddo



!--------------------------------------------------------------------
      if (ipronn.eq."spline-t") then ! time-dep. neutrals or impurities

         ! For NPA:
         do kkk=1,npaproc
            if (npa_process(kkk).ne.'notset'.and.kkk.ne.5) then
               if (itme.eq.0) then
                  do l=1,njene
                  tmpt(l)= ennin_t(l,1,kkk)
                  enddo
               elseif (itme.lt.nbctime) then
                  do l=1,njene
                  tmpt(l)= ennin_t(l,itme,kkk)+ &
                  (ennin_t(l,itme1,kkk)-ennin_t(l,itme,kkk)) &
                  /(bctime(itme1)-bctime(itme))*(timet-bctime(itme))
                  enddo
               else ! itme >= nbctime
                  do l=1,njene
                  tmpt(l)= ennin_t(l,nbctime,kkk)
                  enddo
               endif
               !------------
               call tdinterp("zero","linear",ryain,tmpt,njene, &
                    rya(1:ubound(rya,1)),tr(1),setup0%lrzmax)
               do ll=1,setup0%lrzmax
               enn(ll,kkk)= ennscal(kkk)*tr(ll) ! for tdnpa
               enddo
            endif ! npa_process(kkk).ne.'notset'.and.kkk.ne.5
         enddo ! kkk=1,npaproc
         ! For NPA, process#5 (recombination with electrons):
         if (npa_process(5).eq.'radrecom') then
             k=max(kelecg,kelecm)  !I.e., using bkgrnd distn, if avail.
             do ll=1,setup0%lrzmax
                enn(ll,5)= ennscal(5)*reden(k,ll)
                ! Note: reden is defined above (can be time-dependent)
             enddo
         endif

      endif ! ipronn.eq.'spline-t'
!--------------------------------------------------------------------



!     Electric field or target current
!     Skip, if ampfmod.eq.enabled .and. n+1.ge.nonampf

      if (ampfmod.eq.'enabled' .and. n+1.ge.nonampf) then
        continue ! Meaning: skip elecfld(0)=elecc(1),etc, setting.
         ! Here, profiles is called from tdchief when
         ! n is not updated yet. Here n=0,1,2,...,nstop-1.
         ! For example, if nstop=5 and nonampf=5
         ! (meaning: apply ampf calc. at the last step),
         ! we need to start skipping elecfld(0)=elecc(1),etc
         ! when n+1=nstop=nonampf
      else

      if (efswtch.eq."method1") then
         if (iproelec.eq."prbola-t") then

            if (tmdmeth.eq."method1".and.elecc(1).ne.zero) then

               if (itme.eq.0) then
                  elecfld(0)=elecc(1)
                  elecfld(1)=elecb(1)
                  elecfldc=elecfld(0)
                  elecfldb=elecfld(1)
               elseif (itme.lt.nbctime) then
                  elecfld(0)=elecc(itme)+(elecc(itme1)-elecc(itme)) &
                       /(bctime(itme1)-bctime(itme)) &
                       *(timet-bctime(itme))
                  elecfld(1)=elecb(itme)+(elecb(itme1)-elecb(itme)) &
                       /(bctime(itme1)-bctime(itme)) &
                       *(timet-bctime(itme))
                  elecfldc=elecfld(0)
                  elecfldb=elecfld(1)
               else
                  elecfld(0)=elecc(nbctime)
                  elecfld(1)=elecb(nbctime)
                  elecfldc=elecfld(0)
                  elecfldb=elecfld(1)
               endif

            elseif (tmdmeth.eq."method2") then

               elecfld(0)=tdprof(timet,elecc,bctime)
               elecfld(1)=tdprof(timet,elecb,bctime)
               elecfldc=elecfld(0)
               elecfldb=elecfld(1)

            endif

            dratio=elecfld(1)/elecfld(0)
            do 17 ll=1,setup0%lrzmax
               call profaxis(rn,npwrelec,mpwrelec,dratio,rya(ll))
               elecfld(ll)=elecfld(0)*rn
 17            continue
         endif  !On iproelec.eq.prbola-t

         if (iproelec.eq."spline-t") then
            if (tmdmeth.eq."method1") then

              if (itme.eq.0) then
                  do l=1,njene
                     tmpt(l)=elecin_t(l,1)
                  enddo
              elseif (itme.lt.nbctime) then
                  do l=1,njene
                     tmpt(l)=elecin_t(l,itme) &
                          +(elecin_t(l,itme1)-elecin_t(l,itme)) &
                          /(bctime(itme1)-bctime(itme)) &
                          *(timet-bctime(itme))
                  enddo
               else
                  do l=1,njene
                     tmpt(l)=elecin_t(l,nbctime)
                  enddo
               endif
               call tdinterp("zero","linear",ryain,tmpt,njene,rya(1:ubound(rya,1)), &
                    tr(1),setup0%lrzmax)
               do  ll=1,setup0%lrzmax
                  elecfld(ll)=tr(ll)
               enddo
!              Assume ryain(1)=0 (or close) to get central elecfld
!              needed, e.g., for ampfsoln
               elecfld(0)=tmpt(1)
            endif  !On tmdmeth
            elecfldc=tmpt(1)
            elecfldb=tmpt(njene)

         endif  !On iproelec.eq.spline-t

!100126  Scale elecfld
         do ll=1,setup0%lrzmax
            elecfld(ll)=elecscal*elecfld(ll)
         enddo
         elecfldc=elecscal*elecfldc
         elecfldb=elecscal*elecfldb

      elseif (efswtch.eq."method2" .or. efswtch.eq."method3" .or. &
              efswtch.eq."method4") then

         if (iprocur.eq."prbola-t") then

            if (tmdmeth.eq."method1") then

               if (itme.eq.0) then
                  currxj(0)=xjc(1)
                  currxj(1)=xjb(1)
               elseif (itme.lt.nbctime) then
                  currxj(0)=xjc(itme)+(xjc(itme1)-xjc(itme)) &
                       /(bctime(itme1)-bctime(itme)) &
                       *(timet-bctime(itme))
                  currxj(1)=xjb(itme)+(xjb(itme1)-xjb(itme)) &
                       /(bctime(itme1)-bctime(itme)) &
                       *(timet-bctime(itme))
               else
                  currxj(0)=xjc(nbctime)
                  currxj(1)=xjb(nbctime)
               endif

            elseif (tmdmeth.eq."method2") then

               currxj(0)=tdprof(timet,xjc,bctime)
               currxj(1)=tdprof(timet,xjb,bctime)

            endif



            dratio=currxj(1)/currxj(0)
            do 60 ll=1,setup0%lrzmax
               call profaxis(rn,npwrxj,mpwrxj,dratio,rya(ll))
               currxj(ll)=currxj(0)*rn
 60         continue

         elseif (iprocur.eq."spline-t") then
            if (tmdmeth.eq."method1") then
              if (itme.eq.0) then
                  do l=1,njene
                     tmpt(l)=xjin_t(l,1)
                  enddo
              elseif (itme.lt.nbctime) then
                  do l=1,njene
                     tmpt(l)=xjin_t(l,itme) &
                          +(xjin_t(l,itme1)-xjin_t(l,itme)) &
                          /(bctime(itme1)-bctime(itme)) &
                          *(timet-bctime(itme))
                  enddo
               else
                  do l=1,njene
                     tmpt(l)=xjin_t(l,nbctime)
                  enddo
               endif
               call tdinterp("zero","linear",ryain,tmpt,njene,rya(1:ubound(rya,1)), &
                    tr(1),setup0%lrzmax)
               do  ll=1,setup0%lrzmax
                  currxj(ll)=tr(ll)
               enddo
            endif  !On tmdmeth
         endif  !On iprocur .eq. prbola-t or spline-t

      elseif (efswtch.eq."method5") then

!          Parallel current from eqdsk:
!          Decided to set up currxj from eqdsk in
!          subroutine tdrmshst.  Also, have to
!          make this setup AFTER processing the eqdsk.

      endif  !On efswtch
      endif  !On ampfmod.eq.enabled .and. n.ge.nonampf

      if (ampfmod.eq.'enabled' .and. n.le.nonampf) then !Fill in for plots
         do it=0,nampfmax
            do ll=0,setup0%lrz
               elecfldn(ll,n,it)=elecfld(ll)/300.d0
            enddo
            elecfldn(setup0%lrz+1,n,it)=elecfldb/300.d0
         enddo
      endif

      if (totcrt(1).ne.zero) then

         if (tmdmeth.eq."method1") then

            if (itme.eq.0) then
               totcurtt=totcrt(1)
            elseif (itme.lt.nbctime) then
               totcurtt=totcrt(itme)+(totcrt(itme1)-totcrt(itme)) &
                    /(bctime(itme1)-bctime(itme)) &
                    *(timet-bctime(itme))
            else
               totcurtt=totcrt(nbctime)
            endif

         elseif (tmdmeth.eq."method2") then

            totcurrtt=tdprof(timet,totcrt,bctime)

         endif

!     Renormalize currxj
         currxjtot=0.
         do 80 ll=1,setup0%lrzmax
            currxjtot=currxjtot+darea(ll)*currxj(ll)
 80      continue
         do 81 ll=1,setup0%lrzmax
            currxj(ll)=totcurtt/currxjtot*currxj(ll)
 81      continue

      endif                     !On totcft(1)

!     Toroidal rotation velocity

      if (iprovphi.eq."prbola-t") then

         if (tmdmeth.eq."method1") then
            if (itme.eq.0) then
               vphiplin(0)=vphic(1)
               vphiplin(1)=vphib(1)
            elseif (itme.lt.nbctime) then
               vphiplin(0)=vphic(itme)+(vphic(itme1)-vphic(itme)) &
                    /(bctime(itme1)-bctime(itme)) &
                    *(timet-bctime(itme))
               vphiplin(1)=vphib(itme)+(vphib(itme1)-vphib(itme)) &
                    /(bctime(itme1)-bctime(itme)) &
                    *(timet-bctime(itme))
            else
               vphiplin(0)=vphic(nbctime)
               vphiplin(1)=vphib(nbctime)
            endif

         elseif (tmdmeth.eq."method2") then

            vphiplin(0)=tdprof(timet,vphic,bctime)
            vphiplin(1)=tdprof(timet,vphib,bctime)

         endif

         dratio=vphiplin(1)/vphiplin(0)
         do 90 ll=1,setup0%lrzmax
            call profaxis(rn,npwrvphi,mpwrvphi,dratio,rya(ll))
            vphipl(ll)=vphiscal*vphiplin(0)*rn

 90      continue

      elseif (iprovphi.eq."spline-t") then
         if (tmdmeth.eq."method1") then
            if (itme.eq.0) then
               do l=1,njene
                  tmpt(l)=vphiplin_t(l,1)
               enddo
            elseif (itme.lt.nbctime) then
               do l=1,njene
                  tmpt(l)=vphiplin_t(l,itme) &
                       +(vphiplin_t(l,itme1)-vphiplin_t(l,itme)) &
                       /(bctime(itme1)-bctime(itme)) &
                       *(timet-bctime(itme))
               enddo
            else
               do l=1,njene
                  tmpt(l)=vphiplin_t(l,nbctime)
               enddo
            endif
            call tdinterp("zero","linear",ryain,tmpt,njene,rya(1:ubound(rya,1)), &
                 tr(1),setup0%lrzmax)
            do  ll=1,setup0%lrzmax
               vphipl(ll)=vphiscal*tr(ll)
            enddo
         endif                  !On tmdmeth

      endif                     !On iprovphi


      write(*,*)
      write(*,*)'profiles: time step n, timet=',n,timet

!     From tdoutput.f:
!.......................................................................
!l    0. initialize "radial" mesh for print out
!.......................................................................
      if (pltvs .ne. "psi") then
        ztext=" rovera "
        do l=1,setup0%lrzmax
          ztr(l)=rovera(l)
        enddo
      else
        ztext=pltvs
        do l=1,setup0%lrzmax
          ztr(l)=(equilpsi(0)-equilpsi(l))/equilpsi(0)
        enddo
      endif

!.......................................................................
!l    1.2 density,temperature,magnetic field,etc.
!.......................................................................

!dir$ nextscalar
!
!MPIINSERT_IF_RANK_EQ_0
      do 123 jk=1,ntotal
         if(n.gt.0)then
           ! At n=0, profiles() is called before aingeom(),
           ! and so bthr(ll) is not defined yet.
           ! Then, rban_vth=NaN
           do ll=1,setup0%lrzmax
             qb_mc=bnumb(jk)*charge*bthr(ll)/(fmass(jk)*clight)
             ! Banana width (for v=vthermal, at t-p bndry):
             ! BH171230: For setup0%lrzdiff=enabled, setup0%lrz<setup0%lrzmax, the default
             ! values of itl_(ll>setup0%lrz) can cause out-of-bounds coss.
             if (ll.le.setup0%lrz) then
             rban_vth(ll)= abs(vth(jk,ll)*coss(itl_(ll),ll)/qb_mc) ! cm
             else
                rban_vth(ll)= zero
             endif
           enddo
         else ! n=0
           rban_vth=0.d0 ! not defined yet
           ! Not a big problem (just for a print-out,
           ! and also printed by tdoutput,
         endif
         if (ipronn.eq."disabled".or.jk.ne.nnspec) then
            if(kspeci(2,jk).eq.'general' .and. &
                        iprovphi.ne.'disabled') then
            WRITE(6,9126) jk,kspeci(1,jk),kspeci(2,jk),bnumb(jk),ztext
            WRITE(6,9127) (il,ztr(il),reden(jk,il),temp(jk,il), &
                 vth(jk,il),energy(jk,il),vphipl(il),il=1,setup0%lrzmax)
            else
            WRITE(6,9120) jk,kspeci(1,jk),kspeci(2,jk),bnumb(jk),ztext
            WRITE(6,9121) (il,ztr(il),reden(jk,il),temp(jk,il), &
                 vth(jk,il),energy(jk,il),rban_vth(il), il=1,setup0%lrzmax)
            endif
         elseif (nnspec.eq.jk) then
            if(kspeci(2,jk).eq.'general' .and. &
                        iprovphi.ne.'disabled') then
            WRITE(6,9128) jk,kspeci(1,jk),kspeci(2,jk),bnumb(jk),ztext
            WRITE(6,9129) (il,ztr(il),reden(jk,il),temp(jk,il), &
                 vth(jk,il),energy(jk,il),enn(il,1),vphipl(il), &
                 il=1,setup0%lrzmax)
            else
            WRITE(6,9122) jk,kspeci(1,jk),kspeci(2,jk),bnumb(jk),ztext
            WRITE(6,9123) (il,ztr(il),reden(jk,il),temp(jk,il), &
                 vth(jk,il),energy(jk,il),enn(il,1),il=1,setup0%lrzmax)
            endif
         endif
 123  continue
!MPIINSERT_ENDIF_RANK

 9120 format(/" species no. ",i3,2x,a,a8,"    charge number: ",f6.2 &
        ,/,1x,15("="),//,"  l",4x,a8,5x,"density",4x, &
        "temperature",6x,"vth", 9x,"energy", 4x,"rban_vth(cm)")
 9121 format(i3,1p6e13.5)
 9122 format(/" species no. ",i3,2x,a,a8,"    charge number: ",f6.2 &
        ,/,1x,15("="),//,"  l",4x,a8,5x,"density",4x, &
        "temperature",6x,"vth",9x,"energy",4x,"neutral den")
 9123 format(i3,1p6e13.5)
 9125 format(/," along magnetic field line, at setup0%lrindx=",i3,":",/, &
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

      write(*,*)
      write(*,*)'profiles: zeff(1:setup0%lrzmax)=',(zeff(il), il=1,setup0%lrzmax)
      write(*,*)

!

      return
      end subroutine profiles
!
!
      real(c_double) function tdprof(timet,y,t)
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save
!
!     computes value of tdprof according to:
!     tdprof=y(1), timet < t(1)
!     tdprof=y(2)+0.5*(y(1)-y(2))*(1.+cos((timet-t(1))/(t(2)-t(1))*pi)),
!        timet in (t(1),t(2).
!     tdprof=y(2), timet > t(2)
!
      dimension y(2),t(2)
      data pi/3.141592653589793/

      if (timet.le.t(1)) then
         tdprof=y(1)
      elseif (timet.gt.t(2)) then
         tdprof=y(2)
      else
         tdprof=y(2)+0.5*(y(1)-y(2))* &
              (1.+cos((timet-t(1))/(t(2)-t(1))*pi))
      endif

      return
      end function tdprof

end module profiles_mod
