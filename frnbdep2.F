!
!
!     ONETWO DIVERGENCE
      subroutine frnbdep2(psi,mi,mj,r,z,potsid,mf,rzpat,nrpat,nzpat,me,
     &  mb,sgxn,vbeam,hxfrac,iopt)
      use param_mod

      implicit none
!-----------------------------------------------------------------------
!     calculates neutral beam deposition on (r,z) grid.  grid size
!     determined by nrpat and nzpat, the number of equally spaced
!     elements in the r and z axes.  outputs the 3rd excited state
!     component, rzhex, to file 'beamdep' for postprocessing.
!-----------------------------------------------------------------------
!     ONETWO DIVERGENCE
#ifdef __MPI
      include 'mpilib.h'
#endif

      integer nbdep ! local 
      parameter (nbdep=29)
      real(c_double) :: psi(ki,kj),r(ki),z(kj),potsid(kf),
     &  rzpat(kix2,kjx2,ke,kb),sgxn(kz,ke,kb),vbeam(ke,kb),
     &  hxfrac(ke,kb)
      real(c_double) :: psipat(kix2,kjx2),rr(kix2),zz(kjx2),
     &  rznub(kix2,kjx2,ke,kb),rzhex(kix2,kjx2,ke,kb),pc(kz),
     &  capsig(kz),rzsig(kikj),frac(kz),rzfrac(kikj),
     &  splnwk(kwork),cspln1(kzm1,3),cspln2(kzm1,3),
     &  bpar1(4),bpar2(4),psikp(kikj)
      integer :: isupp(4,ke,kb),inc(kjx2)
!
      integer i,j,ib,ie,n,nout ! local
      integer :: mi,mj,mf, nrpat,nzpat,me,mb,iopt ! arg. in subr.  
      integer igrid,mfm1,ier,ifail ! local (and arg. in called subroutines)
      real(c_double) :: zero,dr,dz ! local

      zero=0.d0

!     set (r,z) grid modification option
      igrid=iopt-1
!
!     zero out arrays
      do 100 i=1,4
        bpar1(i)=0.
        bpar2(i)=0.
 100  continue
!
!     set up integer array to allow vectorization further on
      do 90 i=1,nzpat
        inc(i)=i-1
 90   continue
!
!     get psi on zone centers
      mfm1=mf-1
      do 110 i=1,mfm1
        pc(i)=0.5*(potsid(i)+potsid(i+1))
 110  continue
!
!     if user defined new (r,z) grid, interpolate psi(i,j) onto it
      if(igrid.eq.1)then
        dr=(r(mi)-r(1))/(nrpat-1)
        rr(1)=r(1)
        do 112 i=2,nrpat
          rr(i)=rr(1)+inc(i)*dr
 112    continue
        dz=(z(mj)-z(1))/(nzpat-1)
        zz(1)=z(1)
        do 114 i=2,nzpat
          zz(i)=zz(1)+inc(i)*dz
 114    continue
        call ibcieu1(psi,ki,r,mi,z,mj,rr,nrpat,zz,nzpat,psipat,kix2,
     &    splnwk,ier)
      endif
!
!     begin loop over beam and beam energy
      do 200 ib=1,mb
        do 201 ie=1,me
!
!     get an estimate of the support of rzpat
          call frsuppor(rzpat,nrpat,nzpat,ie,ib,isupp,ifail)
!
!     get spline fits to macroscopic neutral beam attenuation cross
!     sections as a function of psi
          do 120 i=1,mfm1
            capsig(i)=sgxn(i,ie,ib)
 120      continue
          call icsicu1(pc,capsig,mfm1,bpar1,cspln1,kzm1,ier)
!
!     loop over (r,z) points
          do 130 i=isupp(1,ie,ib),isupp(2,ie,ib)
            do 131 j=isupp(3,ie,ib),isupp(4,ie,ib)
              if(rzpat(i,j,ie,ib).ne.zero)then
                n=n+1
                if(igrid.eq.0)then
                  psikp(n)=psi(i,j)
                else
                  psikp(n)=psipat(i,j)
                endif
              endif
 131        continue
 130      continue
          call icsevu1(pc,capsig,mfm1,cspln1,kzm1,psikp,rzsig,n,ier)
          n=0
          do 140 i=isupp(1,ie,ib),isupp(2,ie,ib)
            do 141 j=isupp(3,ie,ib),isupp(4,ie,ib)
              if(rzpat(i,j,ie,ib).ne.zero)then
                n=n+1
                rznub(i,j,ie,ib)=rzpat(i,j,ie,ib)/(vbeam(ie,ib)
     &            *rzsig(n))
              endif
 141        continue
 140      continue
          do 150 i=isupp(1,ie,ib),isupp(2,ie,ib)
            do 151 j=isupp(3,ie,ib),isupp(4,ie,ib)
              rzhex(i,j,ie,ib)=rznub(i,j,ie,ib)*hxfrac(ie,ib)
 151        continue
 150      continue
 201    continue
 200  continue
!
!     output subset of rzhex to 'beamdep'
      nout=nbdep

#ifdef __MPI
      if(mpirank.ne.0) return
#endif

      open(unit=nout,file='beamdep',status='new')
      write(nout,1000) nrpat,nzpat,me,mb
      if(igrid.eq.0)then
        write(nout,1010) r(1),r(mi),z(1),z(mj)
      else
        write(nout,1010) rr(1),rr(nrpat),zz(1),zz(nzpat)
      endif
      do 300 ib=1,mb
        do 301 ie=1,me
          write(nout,1000) (isupp(i,ie,ib),i=1,4)
          write(nout,1010) hxfrac(ie,ib)
          write(nout,1010) ((rzhex(i,j,ie,ib)
     &      ,i=isupp(1,ie,ib),isupp(2,ie,ib))
     &      ,j=isupp(3,ie,ib),isupp(4,ie,ib))
 301    continue
 300  continue
 1000 format(4(5x,i4))
 1010 format(7(2x,e16.7))
      close(unit=nout)
!
      return
      end subroutine frnbdep2
