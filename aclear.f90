module aclear_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast

  !---END USE

!
!

contains

      subroutine aclear
      use param_mod
      use comm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!..................................................................
!     This routine clears some arrays, especially radial arrays, before
!     starting the calculus
!..................................................................

!.......................................................................

!     Scalars
      iupdn=0
      timet=0.d0 ! YuP[2017] added
      n=0 ! YuP[2017] added

!.......................................................................
!l    1. 1-D arrays
!.......................................................................

!     lrza
      do 100 l=1,lrza
        currt(l)=0.0
        currtp(l)=0.0
        curtor(l)=0.0
        curpol(l)=0.0
        qsafety(l)=0.0
        curreq(l)=0.0
        psimx(l)=0.0
        r0geom(l)=0.0
        rgeom(l)=0.0
        zgeom(l)=0.0
        rovs(l)=0.0
        rovsn(l)=0.0
        tauee(l)=0.0
        xlbnd(l)=0.0
        xlndn0(l)=0.0
        bdre(l)=0.0
        bdrep(l)=0.0
        sorpwt(l)=0.0
        xlncurt(l)=0.0
        area(l)=0.0
        areamid(l)=0.0
        volmid(l)=0.0
        currtpz(l)=0.0
        currtz(l)=0.0
        currtzi(l)=0.0
        currtpzi(l)=0.0
        bpolsqaz(l)=0.0
        h_r(l)=0.0
        totcurzi(l)=0.0
        totcurz(l)=0.0
        fpsiz2(l)=0.0
        jparb(l)=0.0
        jparbt(l)=0.0
        jparbp(l)=0.0
        prestp(l)=0.0
        prest(l)=0.0
        dvol(l)=0.0
        darea(l)=0.0
        psyncz(l)=0.0
        vol(l)=0.0
        onovrpz(l,1)=0.0
        onovrpz(l,2)=0.0
        areacon(l)=0.0
        bpolsqa(l)=0.0
        epsicon(l)=0.0
        erhocon(l)=0.0
        fpsi(l)=0.0
        flxavgd(l)=0.0
        psiavg(1,l)=0.0
        psiavg(2,l)=0.0
        rpcon(l)=0.0
        rmcon(l)=0.0
        zpcon(l)=0.0
        zmcon(l)=0.0
        es_bmax(l)=0.0
        bpsi_max(l)=0.0
        bpsi_min(l)=0.0
        lbpsi_max(l)=0
        lbpsi_min(l)=0
        z_bmax(l)=0.0
        bpsi_z_bmax(l)=0.0
        lz_bmax(lrza)=0
        volcon(l)=0.0
        powrft(l)=0.0
        currpar(l)=0.
        zreshin(l)=0.
        tauii(l)=0.
        tau_neo(l)=0.
        drr_gs(l)=0.
        rhol(l)=0.
        rhol_pol(l)=0.
        taubi(l)=0.
        tau_neo_b(l)=0.
        drr_gs_b(l)=0.
        rhol_b(l)=0.
        rhol_pol_b(l)=0.
 100  continue

      call bcast(bscurm,zero,size(bscurm))
      call bcast(bscurmi,zero,(lrza+1)*4)
      call bcast(bscurma,zero,4)

      do 101 l=0,lrza
        sorpwti(l)=0.0
        rya(l)=0.0
        rz(l)=0.0
        rrz(l)=0.0
        ccurtor(l)=0.0
        ccurpol(l)=0.0
        do 102 kk=0,nmodsa
          powurfi(l,kk)=0.0
 102    continue
 101  continue
      rya(lrza+1)=0.0

!     lrorsa
      do 110 l=1,lrorsa
        n_(l)=0        ! YuP[2017]added
        time_(l)=0.d0  ! YuP[2017]added
        consn(l)=0.0
        consn0(l)=0.0
        currmtp(l)=0.0
        currmt(l)=0.0
        sgaint1(l)=0.0
        thb(l)=0.0
        tbnd(l)=0.0
        currmtz(l)=0.0
        currmtpz(l)=0.0
        starnue(l)=0.0
        taueeh(l)=0.0
        vfluxz(l)=0.0
        rovsloc(l)=0.0
        irzplt(l)=0
        sptzr(l)=0.0
 110  continue

!     lsa
!%OS  do 120 l=1,lsa
!%OS  120  continue


!.......................................................................
!l    2. 2-D arrays
!.......................................................................

!     (ntotala,lrza) ; (ntotala,lrorsa)
      do 200 k=1,ntotala
        reden(k,0)=0.0
        temp(k,0)=0.0
        do 201 l=1,lrza
          energy(k,l)=0.0
          reden(k,l)=0.0
          temp(k,l)=0.0
          vth(k,l)=0.0
 201    continue
        denpar(k,0)=0.0
        temppar(k,0)=0.0
        do 202 l=0,lsa1
          enrgypa(k,l)=0.0
          vthpar(k,l)=0.0
          if (cqlpmod .eq. "enabled") then
            denpar(k,l)=0.0
            temppar(k,l)=0.0
          endif
 202    continue
 200  continue

!     (ngena,.)
      do 210 k=1,ngena

!     (ngena,lrza)
        do 215 l=1,lrza
          curr(k,l)=0.0
          hnis(k,l)=0.0
          ratio(k,l)=0.0
          tauegy(k,l)=0.0
          eparc(k,l)=0.0
          eperc(k,l)=0.0
          wpar(k,l)=0.0
          wperp(k,l)=0.0
          xlndn00(k,l)=0.0
          xlncur(k,l)=0.0
          xlndn(k,l)=0.0
          xlndnr(k,l)=0.0
          energyr(k,l)=0.0
          currr(k,l)=0.0
          xlndnv(k,l)=0.0
          energyv(k,l)=0.0
          currv_(k,l)=0.0
          currz(k,l)=0.0
          gkpwrz(k,l)=0.0
          rfpwrz(k,l)=0.0
          pegyz(k,l)=0.0
          pplossz(k,l)=0.0
          denra(k,l)=0.0
          curra(k,l)=0.0
          fdenra(k,l)=0.0
          fcurra(k,l)=0.0
 215    continue

!     (ngena,lrorsa)
        do 216 l=1,lrorsa
          currm(k,l)=0.0
          energym(k,l)=0.0
          jchang(k,l)=0
 216    continue
 210  continue

!     (lrza,nmodsa)
      do 220 l=1,nmodsa
        do 221 k=1,lrza
          powrf(k,l)=0.0
          powrfc(k,l)=0.0
          powrfl(k,l)=0.0
 221    continue
 220  continue

!     (nplota)
      do i=1,nplota
         tplot(i)=-one
         tplt3d(i)=-one
      enddo
!     (nsavea)
      do i=1,nsavea
         tsave(i)=-one
      enddo


      return
      end
end module aclear_mod
