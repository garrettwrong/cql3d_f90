module finit_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use cqlcomm_mod
  use lossorbm_mod, only : lossorbm
  use param_mod
  use tdreadf_mod, only : tdreadf

  !---END USE

contains

      subroutine finit
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!..................................................................
!     set up distributions
!     For n.eq.0, set up initial (relativistic) Maxwellians,
!       restart, or read in an initial distribution.
!     For n.eq.nfpld, add in further loading (a constant
!       plus a shifted "Maxwellian" or a read in distribution,
!       as specified through fpld).
!..................................................................

      data ireturn/0/
!.......................................................................

!.......................................................................
!     f from previous run (restart option), but only call at end
!     of l_=lrors,1,-1 sequence.
!.......................................................................

      if (nlrestrt.ne."disabled") then
         ireturn=1
         call tdreadf(2)
      endif

      if (ireturn .ne. 0) return

!.......................................................................
!     f from file "fpld_dsk"
!     In opening and closing "fpld_dsk", we use that finit
!       is called in the loop do ll=lrors,1,-1 (and ll => l_).
!.......................................................................

      if ((n.eq.0) .and. (fpld(1,1).eq.-1.0)) then
         iunit=20
         if (l_.eq.lrors) then
            open (unit=iunit,file='fpld_dsk',status='old')
         endif
         do 80 k=1,ngen
            read (iunit,1000) ((f(i,j,k,l_),i=1,iy),j=1,jx)
 80      continue
         if (l_.eq.1) then
            close(unit=iunit)
            ireturn=1
         endif
         go to 999
      endif
 1000 format(5(1pe16.7))

!.......................................................................
!     normal initialization
!.......................................................................

      do 70 k=1,ngen

!..................................................................
!     Given the temperature temp(k,lr_) in Kev, a Maxwellian distribution
!     will be generated. If relativ.ne."disabled", the Maxwellian will
!     be relativistic. Note that the distribution will be properly
!     normalized in subroutine diaggnde so that it reflects the desired
!     density.
!..................................................................

        if ((n.eq.0) .and. (fpld(1,1).ne.-1.0)) then
        thta=fmass(k)*clite2/(temp(k,lr_)*ergtkev)
        if (cqlpmod .eq. "enabled") &
          thta=fmass(k)*clite2/(temppar(k,ls_)*ergtkev)
        do 2 j=1,jx
           !XXXXXXXXX whole lotta underflow here. YuP: this is a standard Maxw.distr. Usually ok.
          swwtemp=exp(-gamm1(j)*thta)
          do 3 i=1,iy
            f(i,j,k,l_)=swwtemp
 3        continue
 2       continue
         endif

!..................................................................
!     Additional loading if fpld(1,k).ne.0.
!     fpld(1,k)= a real number.gt.0.:
!       fpld(1,k)=fraction of flux surface averaged density in
!               this additional loading of f (.ge.0., .le.1.).
!       fpld(2,k)=fraction of the additional loading which is
!               in the shifted "Maxwellian" (.ge.0., .le.1.).
!       The shifted "Maxwellian" is an  energy distribution (i.e.,
!       as above) except shifted in energy to fpld(3,k) with
!       energy spread fpld(4,k))*(guassian in theta about radian angle
!       fpld(5,k), with dispersion fpld(6,k)).
!       fpld(7:10,k) windows these functions in energy and theta.
!     fpld(1,k)=-2.0 ("fpld_ds1", previously on C90):
!      Then fpld(2,k) gives fraction of total density (reden) which
!        is distributed according the read-in distribution in
!        "fpld_dsk".  Windowing of this distribution still given
!        by fpld(7:10,k).
!..................................................................

! If fpld(1,k).ne.0, then determine quantities proportional to
! the flux surface averaged densities.
! (fpld(1,k).ne.-1.0 used to be fpld(1,k).ne."fpld_dsk" on C90)

      if (n.eq.nfpld) then
!990225      if (fpld(1,k).ne.0. .or. fpld(1,k).ne.-1.0) then
      if (fpld(1,k).ne.zero .and. fpld(1,k).ne.-1.0) then
         if (fpld(1,k).eq.-2.0) then
            fpld11=fpld(2,k)
            fpld21=1.
            ids1=1
         else
            fpld11=fpld(1,k)
            fpld21=fpld(2,k)
            ids1=0
         endif

         call bcast(tam1,zero,jx)
         call bcast(tam2,zero,jx)
         call bcast(tam3,zero,jx)
!
!  Put constant distribution in temp1, zero temp2
         call bcast(temp1(0:iy+1,0:jx+1),one,iyjx2)
         call bcast(temp2(0:iy+1,0:jx+1),zero,iyjx2)
!
! If fpld(2,k).ne.0., then set up the equatorial plane
! shifted "Maxwellian" in temp2.
         if (fpld21.ne.zero) then
            if (ids1.ne.1) then
               thta1=fmass(k)*clite2/(fpld(4,k)*ergtkev)
               gmpeak=fpld(3,k)*ergtkev/(fmass(k)*clite2)
               do 30 i=1,iy
                  temc1(i)=exp(-(y(i,l_)-fpld(5,k))**2*.5/fpld(6,k))
 30            continue
               do 31 j=1,jx
                  swwtmp=exp(-(gamm1(j)-gmpeak)*thta1)
                  do 32 i=1,iy
                     temp2(i,j)=swwtmp*temc1(i)
 32               continue
 31            continue
            else
               iunit=20
               open (unit=iunit,file='fpld_dsk',status='old')
               read (iunit,1000) ((temp2(i,j),i=1,iy),j=1,jx)
               close(unit=iunit)
            endif

         endif
!
!  Window down constant and shifted "Maxwellian" preloading:

         fmc2=fmass(k)*clite2/ergtkev
         iimin=1
         do 40 i=1,iy
            if (y(i,l_).ge.fpld(9,k)) then
               iimin=i
               go to 41
            endif
 40      continue
 41      continue
         iimax=iy
         do 42 i=iy,1,-1
            if (y(i,l_).le.fpld(10,k)) then
               iimax=i
               go to 43
            endif
 42      continue
 43      continue

         do 33 j=2,jx
            if (gamm1(j)*fmc2 .ge. fpld(7,k) .and. &
                 gamm1(j)*fmc2 .le. fpld(8,k)) then
               do 34 i=1,iimin-1
                  temp2(i,j)=0.
                  temp1(i,j)=0.
 34            continue
               do 35 i=iimax+1,iy
                  temp2(i,j)=0.
                  temp1(i,j)=0.
 35            continue
            else
               do 36 i=1,iy
                  temp2(i,j)=0.
                  temp1(i,j)=0.
 36            continue
            endif
 33      continue

!
!     Get flux surface average densities associated with
!     temp1, temp2, and f.
!     hn =density of relativistic Maxwellian
!     hn0=density of constant distribution
!     hn1=density of shifted "maxwellian"

         do 20 i=1,iy
          do 10 j=1,jx
            tam1(j)=tam1(j)+temp1(i,j)*cynt2(i,l_)* &
              abs(coss(i,lmdpln_))*tau(i,lr_)
            tam2(j)=tam2(j)+temp2(i,j)*cynt2(i,l_)* &
              abs(coss(i,lmdpln_))*tau(i,lr_)
            tam3(j)=tam3(j)+f(i,j,k,l_)*cynt2(i,l_)* &
              abs(coss(i,lmdpln_))*tau(i,lr_)
 10       continue
 20     continue

        hn=0.0
        hn0=0.0
        hn1=0.0
        do 39 j=1,jx
           hn0=hn0+tam1(j)*cint2(j)
           hn1=hn1+tam2(j)*cint2(j)
           hn=hn+tam3(j)*cint2(j)
 39     continue

! Using these densities, then scale distributions according
! to the definitions of fpld(1,k) and fpld(2,k):

        if (fpld21.eq.zero) then
           sp=0.
           htp=hn0
        elseif (fpld21.ne.1.) then
           sp=fpld21*hn0/((1.-fpld21)*hn1)
           htp=hn0/(1.-fpld21)
        else
           sp=1.
           htp=hn1
        endif

        if(fpld11.ne.1.) then
           stp=fpld11*hn/((1.-fpld11)*htp)
        else
           stp=1.
        endif

        if (fpld21.ne.1.) then
           do 51 j=1,jx
              do 50 i=1,iy
                 temp2(i,j)=stp*temp1(i,j)+stp*sp*temp2(i,j)
 50           continue
 51        continue
        else
           do 53 j=1,jx
              do 52 i=1,iy
                 temp2(i,j)=stp*sp*temp2(i,j)
 52           continue
 53        continue
        endif

        if (fpld11.ne.1.) then
           do 55 j=1,jx
              do 54 i=1,iy
                 f(i,j,k,l_)=f(i,j,k,l_)+temp2(i,j)
 54           continue
 55        continue
        else
           do 57 j=1,jx
              do 56 i=1,iy
                 f(i,j,k,l_)=temp2(i,j)
 56           continue
 57        continue
        endif

!     Symmetrize trapped particles
        do 61 j=1,jx
           do 60 i=itl+1,iyh
              iu=iy+1-i
              f(i,j,k,l_)=0.5*(f(i,j,k,l_)+f(iu,j,k,l_))
              f(iu,j,k,l_)=f(i,j,k,l_)
 60        continue
 61     continue

      endif  !on fpld(1,k)
      endif  !On n.eq.nfpld


!..................................................................
!     If this is a mirror  calculation with a loss region
!     delete the particles on loss orbits..
      if (lossmode(k) .eq. "mirrorcc") then
         ! YuP: how about other models, "mirrsnk" and "mirrsnk1" ?
         call lossorbm(ephicc,k)
         do 100 j=1,jx
            do 90 i=1,iy
               if(gone(i,j,k,indxlr_) .lt. -0.9) then
                 f(i,j,k,l_)=zero
               endif
 90         continue
 100     continue
      endif
!..................................................................


 70   continue ! k

 999  return
      end subroutine finit
      
      
end module finit_mod
