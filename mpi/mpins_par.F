!##########################################################################
!# This is just a text file that contains lines
!# for converting CQL3D source *.f90 files into parallel form,
!# using doparallel.py script.
!# The lines between "!MPII***" 
!# and next "!MPII***" (can use !MPIINSERT_  line for convenience)
!# are substituted into *.f90 files in source directory.
!# Note that it is adapted (YuP[2019-05-31]) to have no "c" in the 1st column,
!# and no continuation sign in the 6th column,
!# so it is suitable for *.f90 modular files.
!##########################################################################
!  DO NOT LEAVE EMPTY LINES IN FRONT OF FIRST OCCURENCE OF !MPI*** LINE.
!
!MPIINSERT_INCLUDE
      include 'mpilib.h'
!MPIINSERT_  ! This is just an indicator of the end of insertion section

!MPIINSERT_IF_RANK_NE_0_RETURN
      if(mpirank.ne.0) return
!MPIINSERT_

!MPIINSERT_IF_RANK_EQ_0
      if(mpirank.eq.0) then
!MPIINSERT_

!MPIINSERT_ENDIF_RANK
      endif  ! for if(mpirank.eq.***)
!MPIINSERT_

!MPIINSERT_START
      call MPI_INIT(mpiierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,mpisize,mpiierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,mpirank,mpiierr)
      if(mpirank.eq.0) PRINT *,'MPISIZE ===',mpisize
      if(mpisize.le.1) stop '===   Run with number of cores >1   ==='
      !PRINT *,'Start mpirank=',mpirank
      if(mpirank.eq.0) then
         mpitime = MPI_WTIME()
      endif
!MPIINSERT_

!MPIINSERT_FINISH
      if(mpirank.eq.0) then
         WRITE(*,*) 'MPI Full time =',MPI_WTIME()-mpitime
      endif
      call MPI_FINALIZE(mpiierr)
      !PRINT *,'close_mpi:  mpirank===',mpirank
!MPIINSERT_

!MPIINSERT_MPIWORKER
      if(soln_method.eq.'direct' .and. setup0%lrzmax.gt.1) then
         ! Parallelization for the impavnc0 solver is limited
         ! for soln_method='direct' (for now)
         if(mpisize.gt.1) then
            mpiworker= MOD(ll-1,mpisize-1)+1  !1...(mpisize-1)
         else
            PRINT*, '------- WARNING: mpisize=1 -------'
            mpiworker=0
         endif
      else
         ! In all other cases, perform calculations
         ! for all flux surfaces on mpirank=0, then broadcast results
         mpiworker=0
      endif
      !if(mpirank.eq.mpiworker) then
      !  write(*,*)'n,n_(ll),ll=',n,n_(ll),ll,' mpiworker=',mpiworker
      !endif
!MPIINSERT_MPIWORKER_KRF
         if(mpisize.gt.1) then
            mpiworker= MOD(krf-1,mpisize-1)+1
         else
            PRINT*, '------- WARNING: mpisize=1 -------'
            mpiworker=0
         endif
!MPIINSERT_MPIWORKER_IRAYKRF
         iraykrf= iray + nrayn*(krf-1)
         if(mpisize.gt.1) then
            mpiworker= MOD(iraykrf-1,mpisize-1)+1
         else
            PRINT*, '------- WARNING: mpisize=1 -------'
            mpiworker=0
         endif
      !if(mpirank.eq.mpiworker) then
      !PRINT *,'n,iray,krf,mpiworker=',n,iray,krf,mpiworker
      !endif
!MPIINSERT_MPIWORKER_LFCT
         if(mpisize.gt.1) then
            mpiworker= MOD(lfct-1,mpisize-1)+1
         else
            PRINT*, '------- WARNING: mpisize=1 -------'
            mpiworker=0
         endif
!MPIINSERT_MPIWORKER_L
         if(mpisize.gt.1) then
            mpiworker= MOD(l-1,mpisize-1)+1 ! poloidal grid
         else
            PRINT*, '------- WARNING: mpisize=1 -------'
            mpiworker=0
         endif
!MPIINSERT_MPIWORKER_I
         if(mpisize.gt.1) then
            mpiworker= MOD(i-1,mpisize-1)+1
         else
            PRINT*, '------- WARNING: mpisize=1 -------'
            mpiworker=0
         endif
!MPIINSERT_MPIWORKER_J
         if(mpisize.gt.1) then
            mpiworker= MOD(j-1,mpisize-1)+1
         else
            PRINT*, '------- WARNING: mpisize=1 -------'
            mpiworker=0
         endif
!MPIINSERT_MPIWORKER_IV
         if(mpisize.gt.1) then
            mpiworker= MOD(iv-1,mpisize-1)+1
         else
            PRINT*, '------- WARNING: mpisize=1 -------'
            mpiworker=0
         endif
!MPIINSERT_MPIWORKER_IJ
         !ij= j + jx*(i-1)
         ij= i + iy*(j-1)
         if(mpisize.gt.1) then
            mpiworker= MOD(ij-1,mpisize-1)+1
         else
            PRINT*, '------- WARNING: mpisize=1 -------'
            mpiworker=0
         endif
!MPIINSERT_BARRIER
      call MPI_BARRIER(MPI_COMM_WORLD,mpiierr)
!MPIINSERT_BCAST_RDC_GRID
      call MPI_BCAST(n_uprp,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(n_upar,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(n_psi,1, MPI_INTEGER,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(vc_cgs,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(upar_min,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(upar_max,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
!MPIINSERT_BCAST_RDC
      call MPI_BCAST(rho_a,n_psi,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(uprp,n_uprp,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(upar,n_upar,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(rdc_cqlb,n_uprp*n_upar*n_psi,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(rdc_cqlc,n_uprp*n_upar*n_psi,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(rdc_cqle,n_uprp*n_upar*n_psi,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(rdc_cqlf,n_uprp*n_upar*n_psi,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
!MPIINSERT_BCAST_NRAYN
      call MPI_BCAST(nrayn,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(nray,nmodsa,MPI_INTEGER,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(nharm,nmodsa,MPI_INTEGER,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(freqcy,nmodsa,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(omega,nmodsa,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
!MPIINSERT_BCAST_NRAYELTS
      call MPI_BCAST(nrayelts,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(nrayelt,nrayn*mrfn,MPI_INTEGER,0,MPI_COMM_WORLD,mpiierr)
!MPIINSERT_BCAST_RAYS_DATA
      call MPI_BCAST(nharm1,nmodsa,MPI_INTEGER,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(nharms,nmodsa,MPI_INTEGER,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(jslofas,nrayn*mrfn,MPI_INTEGER,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(nurefls,nrayn*mrfn,MPI_INTEGER,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(keiks,nrayn*mrfn,MPI_INTEGER,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(jpes,nrayn*mrfn,MPI_INTEGER,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(jpis,nrayn*mrfn,MPI_INTEGER,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(istarts,nrayn*mrfn,MPI_INTEGER,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(iprmt5,nrayn*mrfn,MPI_INTEGER,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(jhlfs,nrayn*mrfn,MPI_INTEGER,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(sxxrt,nrayn*mrfn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(skpsi,nrayn*mrfn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(skth,nrayn*mrfn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(skphi,nrayn*mrfn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(delpwr,nrayelts*nrayn*mrfn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(fluxn,nrayelts*nrayn*mrfn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(seikon,nrayelts*nrayn*mrfn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(spsi,nrayelts*nrayn*mrfn, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(sdpwr,nrayelts*nrayn*mrfn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(sbtot,nrayelts*nrayn*mrfn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(sene,nrayelts*nrayn*mrfn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(salphac,nrayelts*nrayn*mrfn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(salphal,nrayelts*nrayn*mrfn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(ws,nrayelts*nrayn*mrfn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(wr,nrayelts*nrayn*mrfn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(wz,nrayelts*nrayn*mrfn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(wnpar,nrayelts*nrayn*mrfn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(wnper,nrayelts*nrayn*mrfn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(wphi,nrayelts*nrayn*mrfn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(wdnpar,nrayelts*nrayn*mrfn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(cwexde,nrayelts*nrayn*mrfn,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(cweyde,nrayelts*nrayn*mrfn,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(cwezde,nrayelts*nrayn*mrfn,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiierr)
!MPIINSERT_BCAST_DISTRIBUTION
      call MPI_BCAST(f,iyjx2*ngen*lrors,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
!MPIINSERT_BCAST_COLL_COEFFS
      call MPI_BCAST(cal,iyjx*ngen*lrors,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(cbl,iyjx*ngen*lrors,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(ccl,iyjx*ngen*lrors,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(cdl,iyjx*ngen*lrors,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(cel,iyjx*ngen*lrors,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(cfl,iyjx*ngen*lrors,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(eal,iyjx*ngen*2*lrors,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(ebl,iyjx*ngen*2*lrors,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
!MPIINSERT_BCAST_SCAL
      call MPI_BCAST(scal,iyjx*ngen*lrors,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
!MPIINSERT_BCAST_VELSOU
      call MPI_BCAST(velsou,iyjx2*ngen*lrors,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
!MPIINSERT_BCAST_ENTR
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

!MPIINSERT_IF_RANK_EQ_MPIWORKER
      if(mpirank.eq.mpiworker) then

!MPIINSERT_SEND_RECV
      if(soln_method.eq.'direct' .and. setup0%lrzmax.gt.1) then
      ! Parallelization for the impavnc0 solver is limited
      ! for soln_method='direct' (for now)
      if(mpirank.eq.0.or.mpirank.eq.mpiworker) then
         call send_data ! send or recv data on f and coll.coeffs.
      endif
      endif

!MPIINSERT_SEND_RECV_ENTR
      if(mpirank.eq.0.or.mpirank.eq.mpiworker) then
         call send_entr(k,lefct)
         !send/recv entr(k,:,l_),pwrrf,pwrrfs(:,k,l_)
      endif

!MPIINSERT_SEND_URFPWR
         !PRINT*,'SEND_urfpwr: mpirank,krf,iray=',mpirank,krf,iray
         ! Count number of elements for this ray
         irayis=0
         do is=1,nrayelt(iray,krf)! Loop over ray elements
            if(ipwr(is).eq.1) irayis=irayis+1 ! incremented up to nrayis
         enddo
         nrayis=irayis
         mpisz=nrayis ! number of elements
         irayis=0
         do is=1,nrayelt(iray,krf)! Loop over ray elements
           if(ipwr(is).eq.1)  then
             irayis=irayis+1 ! incremented up to nrayis
             urftmp(0*mpisz+irayis)= urfpwr(is,iray,krf)
             urftmp(1*mpisz+irayis)= urfpwrc(is,iray,krf)
             urftmp(2*mpisz+irayis)= urfpwrl(is,iray,krf)
             urftmp(3*mpisz+irayis)= scalurf(is,iray,krf)
             urftmp(4*mpisz+irayis)= salphac(is,iray,krf)
           endif
         enddo
         mpitag= iray + nrayn*(krf-1) ! combined: ray-number + wave-mode
         call MPI_SEND(ipwr,nrayelts, MPI_INTEGER1,0,mpitag,MPI_COMM_WORLD,mpiierr)
         call MPI_SEND(urftmp,5*mpisz,MPI_DOUBLE_PRECISION,0,mpitag,MPI_COMM_WORLD,mpiierr)
         !PRINT*,'SEND_urfpwr: krf,iray,mpirank=',krf,iray,mpirank
      !-----------------------------------------------------------

!MPIINSERT_RECV_URFPWR
      if(mpirank.eq.0) then !-------------------------------------------
         !PRINT*,'recv_urfpwr: mpirank,krf,iray=',mpirank,krf,iray
         call MPI_RECV(ipwr,nrayelts,MPI_INTEGER1,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,mpistatus,mpiierr)
         mpitag=mpistatus(MPI_TAG)
         call MPI_RECV(urftmp,nrayelts*5,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,mpitag,MPI_COMM_WORLD,mpistatus,mpiierr)
         mpiray=MOD(mpitag-1,nrayn)+1  ! determine which ray
         mpikrf=(mpitag-mpiray)/nrayn +1 ! determine which krf wave mode
         ! Get mpisz5 (Number of elements received)
         call MPI_GET_COUNT(mpistatus,MPI_DOUBLE_PRECISION,mpisz5,mpiierr)
         mpisz= mpisz5/5 ! urftmp contains 5 arrays
         irayis=0
         do is=1,nrayelt(mpiray,mpikrf)! Loop over ray elements
           if(ipwr(is).eq.1)  then
             irayis=irayis+1 ! incremented up to nrayis==mpisz
             urfpwr(is,mpiray,mpikrf)=  urftmp(0*mpisz+irayis)
             urfpwrc(is,mpiray,mpikrf)= urftmp(1*mpisz+irayis)
             urfpwrl(is,mpiray,mpikrf)= urftmp(2*mpisz+irayis)
             scalurf(is,mpiray,mpikrf)= urftmp(3*mpisz+irayis)
             salphac(is,mpiray,mpikrf)= urftmp(4*mpisz+irayis)
           endif
         enddo
         !PRINT*,'recv_urfpwr: mpikrf,mpiray,mpisz=', &
         !                    mpikrf,mpiray,mpisz
      endif !-----------------------------------------------------------

!MPIINSERT_BCAST_URFPWR
      call MPI_BCAST(urfpwr, nrayelts*nrayn*mrfn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(urfpwrc,nrayelts*nrayn*mrfn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(urfpwrl,nrayelts*nrayn*mrfn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(scalurf,nrayelts*nrayn*mrfn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(salphac,nrayelts*nrayn*mrfn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)

!MPIINSERT_SEND_URFB0
      if(mpirank.eq.mpiworker) then !-----------------------------------
        mpisz=iyjx*setup0%lrz ! number of elements in urfb(i,j,lr)
        call dcopy(mpisz,urfb(1:iy,1:jx,1:setup0%lrz,krf),1,urfbwk(0*mpisz+1),1) !         1 : mpisz
        call dcopy(mpisz,urfc(1:iy,1:jx,1:setup0%lrz,krf),1,urfbwk(1*mpisz+1),1) ! 1*mpisz+1 : 2*mpisz
        ! urfb and urfc are dimensioned as (1:iy,1:jx,1:setup0%lrz,1:mrfn)
        mpisz3=2*mpisz ! the last elem. in above
        urfbwk(mpisz3+0*setup0%lrz+1 : mpisz3+1*setup0%lrz) = powrfl(1:setup0%lrz,krf) !linear damp.
        urfbwk(mpisz3+1*setup0%lrz+1 : mpisz3+2*setup0%lrz) = powrfc(1:setup0%lrz,krf) !coll.damp.
        mpitag= krf ! wave-mode
        call MPI_SEND(urfbwk,mpisz3+2*setup0%lrz,MPI_DOUBLE_PRECISION,0,mpitag,MPI_COMM_WORLD,mpiierr)
      endif !-----------------------------------------------------------
!MPIINSERT_RECV_URFB0
      if(mpirank.eq.0) then !-------------------------------------------
        mpisz=iyjx*setup0%lrz ! number of elements in urfb(i,j,lr)
        mpisz3=2*mpisz ! storage size for urfb,urfc
        call MPI_RECV(urfbwk, mpisz3+2*setup0%lrz,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,mpistatus,mpiierr)
        mpitag=mpistatus(MPI_TAG)
        mpikrf=mpitag ! determine which krf wave mode
        ij=0
        do ll=1,setup0%lrz
        call tdnflxs(lmdpln(ll))
        do j=1,jx
        do i=1,iy
           ij=ij+1
           urfb(i,j,indxlr_,mpikrf)=urfb(i,j,indxlr_,mpikrf)+urfbwk(0*mpisz+ij)
           urfc(i,j,indxlr_,mpikrf)=urfc(i,j,indxlr_,mpikrf)+urfbwk(1*mpisz+ij)
        enddo
        enddo
        enddo ! ll
      endif !-----------------------------------------------------------

!MPIINSERT_BCAST_URFB0
      call MPI_BCAST(urfb,iyjx*setup0%lrz*mrfn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(urfc,iyjx*setup0%lrz*mrfn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)

!MPIINSERT_SEND_EFLUXWK
      if(mpirank.eq.mpiworker) then !-----------------------------------
        call dcopy(mpisz,efluxwk(1:mpisz,i),1,tem2(1:mpisz),1)
        mpitag= i ! i=1,setup0%lrzmax
        call MPI_SEND(tem2, mpisz,MPI_DOUBLE_PRECISION,0,mpitag,MPI_COMM_WORLD,mpiierr)
      endif !-----------------------------------------------------------
!MPIINSERT_RECV_EFLUXWK
      if(mpirank.eq.0) then !-------------------------------------------
        call MPI_RECV(tem2, mpisz,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,mpistatus,mpiierr)
        mpitag=mpistatus(MPI_TAG)
        mpil_=mpitag ! determine which radial surface sent the data
        call dcopy(mpisz,tem2(1:mpisz),1,efluxwk(1:mpisz,mpil_),1)
      endif !-----------------------------------------------------------

!MPIINSERT_SEND_FUS
      if(mpirank.eq.mpiworker) then !-----------------------------------
        mpisz=4
        call dcopy(mpisz, fuspwrv(1:mpisz,lr_), 1,buff(0*mpisz+1),1)
        call dcopy(mpisz, fuspwrm(1:mpisz,lr_), 1,buff(1*mpisz+1),1)
        call dcopy(mpisz, sigf(1:mpisz,lr_),    1,buff(2*mpisz+1),1)
        call dcopy(mpisz, sigm(1:mpisz,lr_),    1,buff(3*mpisz+1),1)
        mpitag= lr_ ! over flux surfaces =1,setup0%lrz
        call MPI_SEND(buff, mpisz*4,MPI_DOUBLE_PRECISION,0,mpitag,MPI_COMM_WORLD,mpiierr)
      endif !-----------------------------------------------------------
!MPIINSERT_RECV_FUS
      if(mpirank.eq.0) then !-------------------------------------------
        mpisz=4
        call MPI_RECV(buff, mpisz*4,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,mpistatus,mpiierr)
        mpitag=mpistatus(MPI_TAG)
        mpil_=mpitag ! determine which radial surface sent the data
        call dcopy(mpisz, buff(0*mpisz+1),1, fuspwrv(1:mpisz,mpil_),1)
        call dcopy(mpisz, buff(1*mpisz+1),1, fuspwrm(1:mpisz,mpil_),1)
        call dcopy(mpisz, buff(2*mpisz+1),1, sigf(1:mpisz,mpil_),   1)
        call dcopy(mpisz, buff(3*mpisz+1),1, sigm(1:mpisz,mpil_),   1)
      endif !-----------------------------------------------------------

!MPIINSERT_BCAST_EFLUX
      call MPI_BCAST(eflux(1:nena,nn),nena,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
!MPIINSERT_

!MPIINSERT_BCAST_EFLUX_NPA
      call MPI_BCAST(eflux_npa(1:nen_npa,nn),nen_npa,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
!MPIINSERT_

!MPIINSERT_BCAST_FUS
      call MPI_BCAST(fuspwrv,4*lrorsa,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(fuspwrm,4*lrorsa,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(sigf,4*lrorsa,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(sigm,4*lrorsa,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
!MPIINSERT_

!MPIINSERT_STARTTIME
      if(mpirank.eq.0) then
         mpitime = MPI_WTIME()
      endif
!MPIINSERT_

!MPIINSERT_ENDTIME
      if(mpirank.eq.0) then
         WRITE(*,*) 'MPI Full time =',MPI_WTIME()-mpitime
      endif
!MPIINSERT_

!MPIINSERT_SUB_SEND_DATA
      subroutine send_data  !used in tdchief.f90 (only)
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      include 'mpilib.h'

      real(c_double), allocatable :: &
              buff(:) ! (mpifsz+10*mpicsz) Buffer for 11 arrays below

      mpifsz= iyjx2*ngen !For send/recv of f(0:iy+1,0:jx+1,1:ngen,lr_),
                         !and velsou(0:iy+1,0:jx+1,1:ngen,lr_)
      mpicsz= iyjx*ngen  !For send/recv of cal(1:iy,1:jx,1:ngen,lr_)
                         !and other collisional coeffs.,
                         !and scal(1:iyjx*ngen,lr_)
      mpisz= mpifsz+11*mpicsz+mpifsz ! buffer size

      if(mpirank.eq.0) then ! receive data from other ranks
         if (.NOT.ALLOCATED(buff)) allocate(buff(mpisz))
         call MPI_RECV(buff, mpisz, &
              MPI_DOUBLE_PRECISION, &
              MPI_ANY_SOURCE, MPI_ANY_TAG, &
              MPI_COMM_WORLD,mpistatus,mpiierr)
         mpil_=mpistatus(MPI_TAG) ! determine which flux surface
         call dcopy(mpifsz,buff(1:mpifsz),1,           f(0:iy+1,0:jx+1,1:ngen,mpil_),1)
         call dcopy(mpicsz,buff(mpifsz+ 0*mpicsz+1),1, cal(1:iy,1:jx,1:ngen,  mpil_),1)
         call dcopy(mpicsz,buff(mpifsz+ 1*mpicsz+1),1, cbl(1:iy,1:jx,1:ngen,  mpil_),1)
         call dcopy(mpicsz,buff(mpifsz+ 2*mpicsz+1),1, ccl(1:iy,1:jx,1:ngen,  mpil_),1)
         call dcopy(mpicsz,buff(mpifsz+ 3*mpicsz+1),1, cdl(1:iy,1:jx,1:ngen,  mpil_),1)
         call dcopy(mpicsz,buff(mpifsz+ 4*mpicsz+1),1, cel(1:iy,1:jx,1:ngen,  mpil_),1)
         call dcopy(mpicsz,buff(mpifsz+ 5*mpicsz+1),1, cfl(1:iy,1:jx,1:ngen,  mpil_),1)
         call dcopy(mpicsz,buff(mpifsz+ 6*mpicsz+1),1, eal(1:iy,1:jx,1:ngen,1,mpil_),1)
         call dcopy(mpicsz,buff(mpifsz+ 7*mpicsz+1),1, eal(1:iy,1:jx,1:ngen,2,mpil_),1)
         call dcopy(mpicsz,buff(mpifsz+ 8*mpicsz+1),1, ebl(1:iy,1:jx,1:ngen,1,mpil_),1)
         call dcopy(mpicsz,buff(mpifsz+ 9*mpicsz+1),1, ebl(1:iy,1:jx,1:ngen,2,mpil_),1)
         call dcopy(mpicsz,buff(mpifsz+10*mpicsz+1),1, scal(1:iyjx*ngen,      mpil_),1)
         !--- Velocity source for radial transport (Note: size=mpifsz)
         call dcopy(mpifsz,buff(mpifsz+11*mpicsz+1),1,velsou(0:iy+1,0:jx+1,1:ngen,mpil_),1)
         !PRINT*,'recv: mpirank,mpil_=',mpirank,mpil_
      else !-> all other ranks send data to rank 0
         if (.NOT.ALLOCATED(buff)) allocate(buff(mpisz))
         call dcopy(mpifsz, f(0:iy+1,0:jx+1,1:ngen,lr_),1,buff(1),1)
         call dcopy(mpicsz, cal(1:iy,1:jx,1:ngen,  lr_),1,buff(mpifsz+ 0*mpicsz+1),1)
         call dcopy(mpicsz, cbl(1:iy,1:jx,1:ngen,  lr_),1,buff(mpifsz+ 1*mpicsz+1),1)
         call dcopy(mpicsz, ccl(1:iy,1:jx,1:ngen,  lr_),1,buff(mpifsz+ 2*mpicsz+1),1)
         call dcopy(mpicsz, cdl(1:iy,1:jx,1:ngen,  lr_),1,buff(mpifsz+ 3*mpicsz+1),1)
         call dcopy(mpicsz, cel(1:iy,1:jx,1:ngen,  lr_),1,buff(mpifsz+ 4*mpicsz+1),1)
         call dcopy(mpicsz, cfl(1:iy,1:jx,1:ngen,  lr_),1,buff(mpifsz+ 5*mpicsz+1),1)
         call dcopy(mpicsz, eal(1:iy,1:jx,1:ngen,1,lr_),1,buff(mpifsz+ 6*mpicsz+1),1)
         call dcopy(mpicsz, eal(1:iy,1:jx,1:ngen,2,lr_),1,buff(mpifsz+ 7*mpicsz+1),1)
         call dcopy(mpicsz, ebl(1:iy,1:jx,1:ngen,1,lr_),1,buff(mpifsz+ 8*mpicsz+1),1)
         call dcopy(mpicsz, ebl(1:iy,1:jx,1:ngen,2,lr_),1,buff(mpifsz+ 9*mpicsz+1),1)
         call dcopy(mpicsz, scal(1:iyjx*ngen,      lr_),1,buff(mpifsz+10*mpicsz+1),1)
         !--- Velocity source for radial transport (Note: size=mpifsz)
         call dcopy(mpifsz,velsou(0:iy+1,0:jx+1,1:ngen,lr_),1,buff(mpifsz+11*mpicsz+1),1)
         mpitag=lr_ ! tag == flux surface number
         call MPI_SEND(buff, mpisz, &
              MPI_DOUBLE_PRECISION, &
              0, mpitag, &
              MPI_COMM_WORLD,mpiierr)
         !PRINT*,'SEND: mpirank,lr_=',mpirank,lr_
      endif
      return
      end subroutine send_data
!MPIINSERT_

!MPIINSERT_SUB_SEND_ENTR
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
!MPIINSERT_

!MPIINSERT_WTIME
      subroutine mpiwtime(s)  ! Not used but can be, when needed (for debugging/optimization)
      character(*) s
      include 'mpilib.h'
      mpitime1 = MPI_WTIME()
      mpitime = mpitime1
      return
      end subroutine mpiwtime
!MPIINSERT_

