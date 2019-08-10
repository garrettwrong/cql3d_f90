module diagimpd_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use aminmx_mod, only : aminmx
  use bcast_mod, only : bcast
  use diagdens_mod, only : diagdens
  use diagentr_mod, only : diagentr
  use diagentr_mod, only : gfi
  use r8subs_mod, only : dcopy

  !---END USE

!
!

contains

      subroutine diagimpd(k)
      use param_mod
      use cqlcomm_mod
      use advnce_mod !here: in diagimpd(). To get gfu(), gfi()
      use r8subs_mod, only : dcopy
      use aminmx_mod, only : aminmx
      implicit integer (i-n), real(c_double) (a-h,o-z)
!..................................................................

!     This routine computes the density gain (loss) due to particle
!     sources (losses) from call to time advancement routine impavnc.
!     Gains are positive, losses negative.
!     The indices refer to the following processes:
!     1- density gained due to setting negative values of
!     distribution function to 0 during velocity sweep.
!     Will be zero for implicit time advancement.
!     Indices 1 and 2 used only with splitting scheme.
!     2- density gained due to setting negative values of
!     distribution function to 0 during theta sweep. (sub diagxswt)
!     Will be zero for implicit time advancement.
!     3- particle source contribution.
!     4- particle loss due to flux at high velocity end.
!     5- loss orbits (part of krook operator)
!     6- toroidal losses(part of krook operator)
!     7- fusion losses (part of krook operator)
!     8- density gained due to setting negative values of
!     distribution to 0 after call to impadv (implicit).
!BH180809: Should add calc of density change due to rescaling?
!BH180809: Done elsewhere? Yuri added plot of added particles at .ps end


!
!     Further the normalized velocity flux at mesh point x(j+.5) is put
!     into vflux(j,k,l_). In the absence of Krook or source terms
!     vflux(jx,k,l_)=1/n * dn/dt *tauee(lr_)  (n is line density).
!

!     Last the routine sets f(i,j,k,l_)=0 if f(i,j,k,l_) < 0
!     if ineg .eq. "enabled".or."trunc_d"
!..................................................................


#ifdef __MPI
      include 'mpilib.h'
#endif

      dimension lefcti(12) ! local working array-index
!..................................................................
!     Call routine to scale distribution if desired.
!..................................................................
!cc      call diagscal(k)  ! YuP: moved outside

!..................................................................
!     Initialize..
!..................................................................


      call dcopy(iyjx2,f(0:iy+1,0:jx+1,k,l_),1,temp2(0:iy+1,0:jx+1),1)
      sgain(1,k)=0.
      sgain(2,k)=0.
      sgain(3,k)=0.
      sgain(4,k)=0.
      sgain(5,k)=0.
      sgain(6,k)=0.
      sgain(7,k)=0.
      sgain(8,k)=0.
      lefcti(1)=-1
      lefcti(2)= 0
      lefcti(3)= 1
      lefcti(4)= 2
      lefcti(5)= 3
      lefcti(6)= 5
      lefcti(7)= 6
      lefcti(8)= 7
      lefcti(9)= 8
      lefcti(10)= 11
      lefcti(11)= 12
      call bcast(tam4,zero,jx)
      call bcast(tam5,zero,jx)
      call bcast(tam6,zero,jx)
!BH081125:  Added zeroing of vflux.  Otherwise, compiler dep results?
      call bcast(vflux(1:jx,k,l_),zero,jx)

!..................................................................
!     Krook operator terms are integrated over velocity space.
!     Simultaneously, compute new density.
!..................................................................

      do 10 i=1,iy
        do 2 j=1,jx
          tam4(j)=tam4(j)+vptb(i,lr_)*temp2(i,j)*cynt2(i,l_)
!
!     Bad orbits..
          tam5(j)=tam5(j)+vptb(i,lr_)*dtreff*gon(i,j)*temp2(i,j)* &
            cynt2(i,l_)
!
!     Toroidal losses
          tam6(j)=tam6(j)-vptb(i,lr_)*dtreff/taulos(i,j,indxlr_) &
            *temp2(i,j)*cynt2(i,l_)
 2      continue
 10   continue

!..................................................................
!     Integrate over x (speed) to obtain relevant diagnostics
!     and the density.
!..................................................................

      sden=0.
      do 3 j=1,jx
        sgain(5,k)=sgain(5,k)+tam5(j)*cint2(j)*one_
        sgain(6,k)=sgain(6,k)+tam6(j)*cint2(j)*one_
        sden=sden+tam4(j)*cint2(j)*one_
!     For diagnostic purposes to compare with tem5 at end of diaggnde:
!     Summing tam7 over j will give flux surface averaged density.
        tam7(j)=tam4(j)*cint2(j)/zmaxpsi(lr_)
 3    continue

!..................................................................
!     Ion source gain..
!..................................................................

      sgain(3,k)=xlncur(k,lr_)*dtreff

!..................................................................
!     Compute:
!     (1) losses at high velocity terminator
!     (2) normalized flux as a function of velocity
!..................................................................
      do 6 i=1,iy
        do 5 j=1,jx
          vflux(j,k,l_)=vflux(j,k,l_)+one_*gfi(i,j,k)*cynt2(i,l_)
 5      continue
 6    continue

!..................................................................
!     Compute density lost at upper terminator
!..................................................................

      sgain(4,k)=vflux(jx,k,l_)*dtreff
!cc       write(*,'(a,i5,e13.5)')'diagimpd sgain(4,k)=', l_,sgain(4,k)

!..................................................................
!     Compute normalized flux
!..................................................................

      do 7 j=1,jx
        vflux(j,k,l_)=vflux(j,k,l_)*tauee(lr_)/xlndn00(k,lr_)
 7    continue

      if (n .eq. 0 .or. n/nchec*nchec .ne. n)  go to 90

!..................................................................
!     determine energy transfer diagnostics
!..................................................................

      entr(k,4,l_)=0.

      do lfct=1,11 ! cpu-intensive calc. for lefct= -1:3; 11 and 12
        lefct=lefcti(lfct) ! == -1:3; 5:8; 11,12;
        !(4 is skipped because it's a sum;  9,10 calc-ed below)
#ifdef __MPI
         if(mpisize.gt.1) then
            mpiworker= MOD(lfct-1,mpisize-1)+1
         else
            PRINT*, '------- WARNING: mpisize=1 -------'
            mpiworker=0
         endif
#endif
#ifdef __MPI
      if(mpirank.eq.mpiworker) then

#endif
        call diagentr(lefct,k)
#ifdef __MPI
      endif  ! for if(mpirank.eq.***)
#endif
#ifdef __MPI
      if(mpirank.eq.0.or.mpirank.eq.mpiworker) then
         call send_entr(k,lefct)
         !send/recv entr(k,:,l_),pwrrf,pwrrfs(:,k,l_)
      endif

#endif
      enddo ! lfct

#ifdef __MPI
      call MPI_BARRIER(MPI_COMM_WORLD,mpiierr)
#endif
#ifdef __MPI
      call MPI_BCAST(entr(k,-1,l_),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(entr(k,0,l_),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(entr(k,1,l_),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(entr(k,2,l_),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(entr(k,3,l_),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(entr(k,4,l_),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(entr(k,5,l_),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(entr(k,6,l_),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(entr(k,7,l_),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(entr(k,8,l_),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(entr(k,11,l_),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(entr(k,12,l_),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(pwrrf(1,k,l_),jx,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)

#endif

      if (k.eq.kelecg) sorpw_rf(kelecg,lr_)=entr(k,3,l_) !-YuP
      if (k.eq.kiong(k)) sorpw_rf(k,lr_)=entr(k,3,l_) !-YuP RF source

 90   continue
      s=0.

!..................................................................
!     if ineg="enabled" or "enabled1" or "trunc_d", set f.le.0 to zero.
!..................................................................

      call dcopy(iyjx2,temp2(0:iy+1,0:jx+1),1,temp1(0:iy+1,0:jx+1),1)
      call bcast(temp4,zero,iyjx2)

      if (ineg .eq. "disabled") go to 141

      if (ineg.eq."enabled" .or. ineg.eq."enabled1" &
           .or.ineg.eq."trunc_d") then
      ! [YuP-101228: Was there an error? (ineg.ne."enabled")]
      !------------------------------------------------------------
      ! Before 03-2011: Set f to 0 for points where f(i,j)<0.
      ! After  03-2011: If a point f(i,j).le.0 is found,
      ! f(i,jj) is set to 0 for ALL jj.ge.jmin_nosource
      ! where jmin_nosource is the smallest j above which
      ! there is no source (for any i).
      ! For points below jmin_nosource, if f<0 is found,
      ! set f to zero for this point (as done before).
      ! As of 20110408, ineg.eq."enabled1" sets to zero
      ! as for "enabled", except don't check for source points.
      !------------------------------------------------------------
      jmin_nosource=1
      if (ineg.ne."enabled1") then
                !If source is present, jmin_nosource.gt.1 will be found.
      do j=jx,1,-1 ! Scan from largest v towards zero
        source_j=sum(abs(source(:,j,k,lr_))) ! Sum over i (all i-range)
        if (source_j .ne.zero) then !source is present for this j-level
          jmin_nosource= j+1 ! Smallest j above which there is no source
          goto 409
        endif
      enddo
 409  continue

      do 411 i=1,iy
      do 410 j=1,jx
         if ( temp2(i,j) .le. zero) then ! Found f(i,j)<=0.
            ! Set f to zero at this local point.
            temp1(i,j)=zero ! Set this f(i,j) to zero.
            temp4(i,j)=-temp2(i,j) ! For diagnostics.
            j0=max(j,jmin_nosource)
            do jj= j0,jx ! Scan all jj>=jmin_nosource
               temp1(i,jj)=zero ! Set f(i,jj) to zero.
               temp4(i,jj)=-temp2(i,jj) ! For diagnostics.
            enddo
         endif
 410  continue
 411  continue

      else ! on ineg

      f011=temp2(1,1)
      fmin=1.d-20*f011 !YuP[2019-04-23] was 1.e-20*f011
      do i=1,iy
      do j=1,jx
         if ( temp2(i,j) .le. fmin) then ! Found f(i,j)<=fmin
            ! Set f to zero at this local point.
            temp1(i,j)=zero ! Set this f(i,j) to zero.
            temp4(i,j)=-temp2(i,j) ! For diagnostics.
            j0=max(j,jmin_nosource)
            do jj= j0,jx ! Scan all jj>=jmin_nosource
               temp1(i,jj)=zero ! Set f(i,jj) to zero.
               temp4(i,jj)=-temp2(i,jj) ! For diagnostics.
            enddo
         endif
      enddo ! On i
      enddo ! On j

      endif  ! On ineg.ne."enabled1"
!
!     routine diagdens will compute density gained by setting negative
!     values of distribution to 0.
!
      call diagdens(xline,xmidp,eline)
      sgain(8,k)=xline*one_
      engain(k)=eline*one_

      elseif(ineg.eq."renorm") then

!     find min and max of f
        fmin=ep100
        fmax=-ep100
        do 500 j=1,jx
          call aminmx(temp1(1:jx,j),1,jx,1,ffmin,ffmax,jmin,jmax)
!990131          fmin=amin1(ffmin,fmin)
!990131          fmax=amax1(ffmax,fmax)
          fmin=min(ffmin,fmin)
          fmax=max(ffmax,fmax)
 500    continue

        if(fmin.lt.zero) then ! YuP[2019-04-23] was if(fmin.lt.0.)
!     calculate line average density of f
          call dcopy(iyjx2,temp1(0:iy+1,0:jx+1),1,temp4(0:iy+1,0:jx+1),1)
          call diagdens(xline0,xmidp0,eline0)
!     add 1.1*fmin to f
          fadd=1.1*abs(fmin)
          do 520 j=1,jx
            do 521 i=1,iy
              temp1(i,j)=fadd+temp1(i,j)
 521        continue
 520      continue
!     recalculate density
          call dcopy(iyjx2,temp1(0:iy+1,0:jx+1),1,temp4(0:iy+1,0:jx+1),1)
          call diagdens(xline1,xmidp1,eline1)
!         renormalize f to old density
          ffac=xline0/xline1
          do 530 j=1,jx
            do 531 i=1,iy
              temp1(i,j)=ffac*temp1(i,j)
 531        continue
 530      continue
!     change of line density (/2) and line energy density
          sgain(8,k)=(xline1-xline0)*one_
          engain(k)=(eline1-eline0)*one_
        endif ! on fmin

      endif  ! on ineg


 141  continue

#ifdef __MPI
      if(mpirank.eq.0) then
#endif
      WRITE(*,'(a,2i4,2e15.7)') &
       'diagimpd BEFORE f update. n,l_,MIN(temp2),MAX(temp2)', &
         n,l_,MINVAL(temp2),MAXVAL(temp2)
      WRITE(*,'(a,2i4,2e15.7)') &
       'diagimpd neg.val.adjusted n,l_,MIN(temp1),MAX(temp1)', &
         n,l_,MINVAL(temp1),MAXVAL(temp1)
#ifdef __MPI
      endif  ! for if(mpirank.eq.***)
#endif
      call dcopy(iyjx2,temp1(0:iy+1,0:jx+1),1,f(0:iy+1,0:jx+1,k,l_),1)
#ifdef __MPI
      if(mpirank.eq.0) then
#endif
      WRITE(*,'(a,2i4,3e15.7)') &
       'diagimpd AFTER  f update. n,l_,MIN(f),MAX(f),SUM(f)=', &
                n,l_,MINVAL(f),MAXVAL(f),SUM(f)
#ifdef __MPI
      endif  ! for if(mpirank.eq.***)
#endif
      !!!if(l_.eq.setup0%lrz)  pause

      if (n .gt. 0 .and. n/nchec*nchec .eq. n) then
        call diagentr(9,k)
        call diagentr(10,k)
      endif

      return
      end subroutine diagimpd

#ifdef __MPI
      subroutine send_entr(k,lefct)  !used in diagimpd.f90 (only)
      use param_mod
      use cqlcomm_mod
      !send/recv entr(k,lefct,l_),pwrrf(1:jx,k,l_),pwrrfs(1:jx,k,l_)
      implicit integer (i-n), real(c_double) (a-h,o-z)
      include 'mpilib.h'
      dimension buff(1+jx)
      if(mpirank.eq.0) then ! receive data from other ranks
         call MPI_RECV(buff, 1+jx, &
              MPI_DOUBLE_PRECISION, &
              MPI_ANY_SOURCE, MPI_ANY_TAG, &
              MPI_COMM_WORLD,mpistatus,mpiierr)
         mpitag=mpistatus(MPI_TAG)
         lefct_=mpitag-2 ! determine which lefct was sent
         entr(k,lefct_,l_)=buff(1) ! for a given lefct
         entr(k,4,l_)=entr(k,4,l_)+buff(1) ! sum
         if (lefct_.eq.3) then
           call dcopy(jx,buff(2:jx+1),1,pwrrf(1:jx,k,l_),1)
           pwrrfs(1,k,l_)=dx(1)*pwrrf(1,k,l_)
           do j=2,jx  ! sum over j
             pwrrfs(j,k,l_)=pwrrfs(j-1,k,l_)+dx(j)*pwrrf(j,k,l_)
           enddo
         endif
      else !-> all other ranks send data to rank 0
         buff(1)=entr(k,lefct,l_) ! for a given lefct
         mpisz=1
         if (lefct.eq.3) then
           call dcopy(jx,pwrrf(1:jx,k,l_),1,buff(2:jx+1),1)
           mpisz=1+jx
         endif
         mpitag=lefct+2 ! 2 added to make mpitag>0 (lefct can be -1)
         call MPI_SEND(buff, mpisz, &
              MPI_DOUBLE_PRECISION, &
              0, mpitag, &
              MPI_COMM_WORLD,mpiierr)
      endif
      return
      end subroutine send_entr
#endif

end module diagimpd_mod
