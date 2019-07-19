module wptramu_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use r8subs_mod, only : dcopy
  use wpbdry_mod, only : wpbdry
  use wpcheck_mod, only : wpcheck
  use znonsym_mod, only : nonsym

  !---END USE

!
!

contains

      subroutine wptramu
      use param_mod
      use cqlcomm_mod
      use r8subs_mod, only : cvmgt, dcopy

      implicit integer (i-n), real(c_double) (a-h,o-z)
!..............................................................
!     Time advancement for parallel transport. (meshy="fixed_mu" case)
!     Assume non-constant total number of y mesh points vs. l: iy_(l).ne.cst
!     Assume l_lower(i)=1, for all i (checked in wpinitl)
!..............................................................

      parameter(lmatsa=(lsa+1)*nbanda)

      dimension zmat(lmatsa),zxdumy(lsa+1)
      dimension zmat2(lmatsa),zxdumy2(lsa+1)
      dimension ism1(0:lsa+1), isdiag(0:lsa+1), isp1(0:lsa+1)
      dimension imaxper(0:lsa+1), iymxper(0:lsa+1), iyhper(0:lsa+1)

      include 'wpadvnc.h'
      dithta(i,j,l)=0.5-0.5*(1/(i+1))+0.5*(1/(iy_(l)+1-i))
      fpithta(i,j,k,l)=fnhalf(i+1,j,k,l)*(1.-dithta(i,j,l)) + &
        fnhalf(i  ,j,k,l)*dithta(i,j,l)

!.......................................................................
!     0. The parameter nummods determines the numerical model used to solve
!     the equation along s. nummods is divided into classes of multiple
!     of 10: numclas=nummods/10 => [0,9]; [10,19]; [20,29] ...
!     numclas = 0: normal case
!     1: assumes up/down symmetric case and periodic condition
!     but compute equil. parameters on top half cross-section
!     on setup0%ls/2+1 points and then copy the rest in sub. wploweq
!     Otherwise, same as numclas=0 with sbdry="periodic"
!     In each class, nummods is divided into two groups:
!     mod(nummods,10) < 5: 2-D FP velocity equ. solved on nodal s points
!     > 4: velocity eq. solved at s+1/2
!     Then, in these two groups one has the following options,
!     with numindx=mod(mod(nummods,10),5)=0,1,2,3 or 4, we have:
!     If numindx = 0,1: use backward diff. scheme + bound. cond at l=1
!     2: "  back/forward depending on sign(cos(theta))
!     and on numixts
!     3: " forward diff. + bound cond at l=setup0%ls
!     4: " centered diff. + bound. cond. at l=1
!     Note: if sbdry="periodic" there might not be any extra boundary
!     conditions imposed (see wpbdry)
!.......................................................................

!.......................................................................
!l    1. Set the index shift flags and the bandwidths for solving
!     the equation as a function of s only (distance along the magnetic
!     field B).
!     If periodic boundary conditions are used, the unknowns are
!     numbered differently (see wpinitl), which doubles the
!     effective bandwidths.
!     For coss(theta)>0 : df/ds -> (f(s+ipshft2)-f(s+ipshft1))/ds_eff
!     <0 : df/ds -> (f(s+imshft2)-f(s+imshft1))/ds_eff
!.......................................................................

      ileft=1
      iright=1
      if (sbdry .eq. "periodic") then
        ileft=2
        iright=2
      endif
      iband=ileft+iright+1

!     standard case: numindx=0 or 1
      ipshft2=0
      ipshft1=-1
      imshft2=0
      imshft1=-1
!     l loop between 1+i.strt and l_upper(i)+i.end, with .=p or m for i< or >iyh
      ipstrt=+1
      ipend=0
      imstrt=+1
      imend=0
      if (numindx .eq. 2) then
        ipshft2=(numixts+1)/2
        ipshft1=(numixts-1)/2
        imshft2=-(numixts-1)/2
        imshft1=-(numixts+1)/2
        ipstrt=-(numixts-1)/2
        ipend=-(numixts+1)/2
        imstrt=(numixts+1)/2
        imend=(numixts-1)/2
      else if (numindx .eq. 3) then
        ipshft2=+1
        ipshft1=0
        imshft2=+1
        imshft1=0
        ipstrt=0
        ipend=-1
        imstrt=0
        imend=-1
      else if (numindx .eq. 4) then
        ipshft2=+1
        ipshft1=-1
        imshft2=+1
        imshft1=-1
        ipstrt=+1
        ipend=-1
        imstrt=+1
        imend=-1
      endif
      if (sbdry .eq. "periodic") then
        ipstrt=0
        ipend=0
        imstrt=0
        imend=0
      endif
      zip1p2=ipshft1+ipshft2
      zim1p2=imshft1+imshft2
      ztrdrt=0.5/dtreff
      ip1t1m2=ipshft1*(1-ipshft2)
      im1t1m2=imshft1*(1-imshft2)
      if (lmidvel .eq. 0) ip1t1m2=0
      if (lmidvel .eq. 0) im1t1m2=0

!.......................................................................
!l    2. Construct the matrix. The equations are solved at a given species k
!     and theta i, but "simultaneously" for all momentum j. This
!     way, one can have longer inner loops over j than over l, as in
!     general jx > setup0%ls. The index j was chosen instead of i, because
!     y(i,l) mesh does not depend on j. But note that we consider the
!     equ. for i and ii=iymax+1-i "simultaneously".
!.......................................................................

!.......................................................................
!l    2.1 Initialize f_n+1/2, velsou and imaxper, iymxper, iyhper
!.......................................................................

      do 210 l=1,setup0%ls
        imaxper(l)=imax(l,lr_)
        iymxper(l)=iy_(l)
        iyhper(l)=iyh_(l)
 210  continue
      call dcopy(iyjx2*ngen*lrors,f(0:iy+1,0:jx+1,1:ngen,1:lrors),1, &
                             fnhalf(0:iy+1,0:jx+1,1:ngen,1:lrors),1)
      if (sbdry .eq. "periodic") then
        call dcopy(iyjx2*ngen,fnhalf(0:iy+1,0:jx+1,1:ngen,1   ),1, &
                              fnhalf(0:iy+1,0:jx+1,1:ngen,setup0%ls+1),1)
        call dcopy(iyjx2*ngen,fnhalf(0:iy+1,0:jx+1,1:ngen,setup0%ls  ),1, &
                              fnhalf(0:iy+1,0:jx+1,1:ngen,0   ),1)
        call dcopy(iyjx2*ngen,velsou(0:iy+1,0:jx+1,1:ngen,1   ),1, &
                              velsou(0:iy+1,0:jx+1,1:ngen,setup0%ls+1),1)
        call dcopy(iyjx2*ngen,velsou(0:iy+1,0:jx+1,1:ngen,setup0%ls  ),1, &
                              velsou(0:iy+1,0:jx+1,1:ngen,0   ),1)
        imaxper(0)=imaxper(setup0%ls)
        imaxper(setup0%ls+1)=imaxper(1)
        iymxper(0)=iymxper(setup0%ls)
        iymxper(setup0%ls+1)=iymxper(1)
        iyhper(0)=iyhper(setup0%ls)
        iyhper(setup0%ls+1)=iyhper(1)
      else
        call dcopy(iyjx2*ngen,fnhalf(0:iy+1,0:jx+1,1:ngen,1   ),1, &
                              fnhalf(0:iy+1,0:jx+1,1:ngen,0   ),1)
        call dcopy(iyjx2*ngen,fnhalf(0:iy+1,0:jx+1,1:ngen,setup0%ls  ),1, &
                              fnhalf(0:iy+1,0:jx+1,1:ngen,setup0%ls+1),1)
        call dcopy(iyjx2*ngen,velsou(0:iy+1,0:jx+1,1:ngen,1   ),1, &
                              velsou(0:iy+1,0:jx+1,1:ngen,0   ),1)
        call dcopy(iyjx2*ngen,velsou(0:iy+1,0:jx+1,1:ngen,setup0%ls  ),1, &
                              velsou(0:iy+1,0:jx+1,1:ngen,setup0%ls+1),1)
        imaxper(0)=imaxper(1)
        imaxper(setup0%ls+1)=imaxper(setup0%ls)
        iymxper(0)=iymxper(1)
        iymxper(setup0%ls+1)=iymxper(setup0%ls)
        iyhper(0)=iyhper(1)
        iyhper(setup0%ls+1)=iyhper(setup0%ls)
      endif

!.......................................................................
!l    2.1.2 Might need value of fnhalf at l=l+1/2. Node value still in f
!     (Already the case if mod(nummods,10)>4, as 2-D momentum eq.
!     solved at l+1/2)
!.......................................................................

      if (mod(nummods,10).le.4 .and. lmidvel.ne.0) then
        do 2120 k=1,ngen
          do 2121 j=1,jx
            do 2122 l=1,setup0%ls
              do 2123 i=1,min(iyhper(l+1),iyhper(l))
                fnhalf(i,j,k,l)=0.5*(fnhalf(i,j,k,l)+fnhalf(i,j,k,l+1))
                ii=iymxper(l)+1-i
                fnhalf(ii,j,k,l)=0.5*(fnhalf(ii,j,k,l)+ &
                  fnhalf(iymxper(l+1)+1-i,j,k,l+1))
 2123         continue
 2122       continue
 2121     continue
 2120   continue
        if (sbdry .eq. "periodic") then
          do 2124 k=1,ngen
            do 2125 j=1,jx
              do 2126 i=1,iyhper(setup0%ls)
                fnhalf(i,j,k,0)=0.5*(fnhalf(i,j,k,setup0%ls)+fnhalf(i,j,k,1))
                iils=iymxper(setup0%ls)+1-i
                ii1=iymxper(1)+1-i
                fnhalf(iils,j,k,0)=0.5*(fnhalf(iils,j,k,setup0%ls)+ &
                  fnhalf(ii1,j,k,1))
                fnhalf(i,j,k,setup0%ls+1)=fnhalf(i,j,k,0)
                fnhalf(ii1,j,k,setup0%ls+1)=fnhalf(iils,j,k,0)
 2126         continue
 2125       continue
 2124     continue
        else
          call dcopy(iyjx2*ngen,fnhalf(0:iy+1,0:jx+1,1:ngen,1   ),1, &
                                fnhalf(0:iy+1,0:jx+1,1:ngen,0   ),1)
          call dcopy(iyjx2*ngen,fnhalf(0:iy+1,0:jx+1,1:ngen,setup0%ls  ),1, &
                                fnhalf(0:iy+1,0:jx+1,1:ngen,setup0%ls+1),1)
        endif
      endif
      ivelmid=1
      if (mod(nummods,10).le.4 .and. lmidvel.ne.0) ivelmid=0

!.......................................................................
!l    2.1.3 Pre-determine matrix row and band index for each l
!.......................................................................
      do 213 l=1,lrors
        ism1(l)=1
        isdiag(l)=ileft+1
        isp1(l)=iband
        if (sbdry.eq."periodic" .and. (l.le.2.or. l.gt.setup0%ls/2)) then
          if (l .eq. 1) then
            ism1(l)=iband
            isp1(l)=iband-1
          else if (l .eq. 2) then
            ism1(l)=ileft
          else if (l .eq. setup0%ls/2+1) then
            isp1(l)=ileft
          else if (l .eq. setup0%ls/2+2) then
            ism1(l)=isdiag(l)+1
            isp1(l)=1
          else
            ism1(l)=iband
            isp1(l)=1
          endif
        endif
 213  continue

!.......................................................................
!l    2.2 Loop over general species k
!.......................................................................

      do 220 k=1,ngen

!     j=1 => zeroth order eq. to be solved, assuming velsou(j=1)=0
!     => f_n+1(j=1)=f_n+1/2(j=1)
        do 221 ll=1,setup0%ls
          call bcast(fnp1(0:iyp1,1,k,ll),f(1,1,k,ll),iyp1+1)
 221    continue

!.......................................................................
!l    2.3 Loop over theta mesh
!.......................................................................

        do 230 i=1,iyhper(1)

!     special case if only one point: l_lower(i)=l_upper(i)=1
          if (l_upper(1) .eq. 1) then
            do 610 j=2,jx
              fnp1(i,j,k,1)=f(i,j,k,1)
              fnp1(iy_(1)+1-i,j,k,1)=f(iy_(1)+1-i,j,k,1)
 610        continue
            go to 230
          endif

          ii=iymax+1-i
          do 2301 l=1,setup0%ls
            do 2302 icol=1,iband
              do 2303 j=2,jx
                bndmats(j,l,icol,1)=0.0
                bndmats(j,l,icol,2)=0.0
 2303         continue
 2302       continue
            do 2304 j=2,jx
              rhspar(l,j,1)=0.0
              rhspar(setup0%ls-2+l,j,1)=0.0
              rhspar(l,j,2)=0.0
 2304       continue
 2301     continue

!.......................................................................
!l    2.3.1 Construct matrices at each j , l and setup0%ls+2-l (if periodic)
!.......................................................................

          do 231 l=1+ipstrt*imstrt,l_upper(i)-ipend*imend
            ilsrow=lsbtopr(l)

!.......................................................................
!l    2.3.1.1 cos(theta) > 0
!.......................................................................

            if (l.lt.1+ipstrt .or. l.gt.l_upper(i)+ipend .or. &
              i.gt.min(imaxper(l+ipshft1),imaxper(l+ipshft2))) &
              go to 2312

            ztrders=vnorm*coss(i,l)/(ipshft2*dszp5(l)-ipshft1*dszm5(l))
            do 23110 j=2,jx
              bndmats(j,ilsrow,ism1(l),1)=DBLE(ipshft1)* &
                (x(j)*ztrders+lmidvel*zip1p2*ztrdrt)
              bndmats(j,ilsrow,isdiag(l),1)=(2.-lmidvel)*ztrdrt- &
                zip1p2*x(j)*ztrders
              bndmats(j,ilsrow,isp1(l),1)=DBLE(ipshft2)* &
                (x(j)*ztrders+lmidvel*zip1p2*ztrdrt)
              rhspar(ilsrow,j,1)=cvmgt( &
                wprhsmd(ilpm1ef(i,l,ip1t1m2),j,k,lpm1eff(l,ip1t1m2)) &
                ,wprhs(ilpm1ef(i,l,ip1t1m2),j,k,lpm1eff(l,ip1t1m2)) &
                ,ivelmid.eq.1 .or. (ipshft1+ipshft2).eq.0)
23110       continue

!.......................................................................
!l    2.3.1.2 cos(theta) < 0
!.......................................................................

 2312       continue
            if (l.lt.1+imstrt .or. l.gt.l_upper(i)+imend .or. &
              i.gt.min(imaxper(l+imshft1),imaxper(l+imshft2))) &
              go to 2313

            iieff=iymxper(l)+1-i
            iieff12=iymxper(l+im1t1m2)+1-i
            ztrders=vnorm*coss(iieff,l) &
              /(imshft2*dszp5(l)-imshft1*dszm5(l))
            do 23120 j=2,jx
              bndmats(j,ilsrow,ism1(l),2)=DBLE(imshft1)* &
                (x(j)*ztrders+lmidvel*zim1p2*ztrdrt)
              bndmats(j,ilsrow,isdiag(l),2)=(2.-lmidvel)*ztrdrt- &
                zim1p2*x(j)*ztrders
              bndmats(j,ilsrow,isp1(l),2)=DBLE(imshft2)* &
                (x(j)*ztrders+lmidvel*zim1p2*ztrdrt)
              rhspar(ilsrow,j,2)=cvmgt(wprhsmd(iieff12,j,k,l+im1t1m2) &
                ,wprhs(iieff12,j,k,l+im1t1m2) &
                ,ivelmid.eq.1 .or. (imshft1+imshft2).eq.0)
23120       continue

!.......................................................................
!l    2.3.1.3 point l=setup0%ls+2-l, i.e. bottom half cross-section, if not
!     yet treated
!.......................................................................

 2313       if (sbdry.ne."periodic" .or. l_upper(i).eq.setup0%ls .or. l.eq.1) &
              go to 231

            ll=setup0%ls+2-l
            ilsrow=lsbtopr(ll)

!.......................................................................
!l    2.3.1.3.1 cos(theta) > 0, ll=setup0%ls+2-l
!.......................................................................

            if (i.gt.min(imaxper(ll+ipshft1),imaxper(ll+ipshft2))) &
              go to 2314

            ztrders=vnorm*coss(i,ll) &
              /(ipshft2*dszp5(ll)-ipshft1*dszm5(ll))
            do 23131 j=2,jx
              bndmats(j,ilsrow,ism1(ll),1)=DBLE(ipshft1)* &
                (x(j)*ztrders+lmidvel*zip1p2*ztrdrt)
              bndmats(j,ilsrow,isdiag(ll),1)=(2.-lmidvel)*ztrdrt- &
                zip1p2*x(j)*ztrders
              bndmats(j,ilsrow,isp1(ll),1)=DBLE(ipshft2)* &
                (x(j)*ztrders+lmidvel*zip1p2*ztrdrt)
              rhspar(ilsrow,j,1)=cvmgt(wprhsmd(i,j,k,ll+ip1t1m2) &
                ,wprhs  (i,j,k,ll+ip1t1m2) &
                ,ivelmid.eq.1 .or. (ipshft1+ipshft2).eq.0)
23131       continue

!.......................................................................
!l    2.3.1.3.2 cos(theta) < 0, ll=setup0%ls+2-l
!.......................................................................

 2314       continue
            if (i.gt.min(imaxper(ll+imshft1),imaxper(ll+imshft2))) &
              go to 231
            iieff=iymxper(ll)+1-i
            iieff12=iymxper(ll+im1t1m2)+1-i
            ztrders=vnorm*coss(iieff,ll)/(imshft2*dszp5(ll)- &
              imshft1*dszm5(ll))
            do 23132 j=2,jx
              bndmats(j,ilsrow,ism1(ll),2)=DBLE(imshft1)* &
                (x(j)*ztrders+lmidvel*zim1p2*ztrdrt)
              bndmats(j,ilsrow,isdiag(ll),2)=(2.-lmidvel)*ztrdrt- &
                zim1p2*x(j)*ztrders
              bndmats(j,ilsrow,isp1(ll),2)=DBLE(imshft2)* &
                (x(j)*ztrders+lmidvel*zim1p2*ztrdrt)
              rhspar(ilsrow,j,2)=cvmgt(wprhsmd(iieff12,j,k,ll+im1t1m2) &
                ,wprhs(iieff12,j,k,ll+im1t1m2) &
                ,ivelmid.eq.1 .or. (imshft1+imshft2).eq.0)
23132       continue

!.......................................................................
!     end of construction of matrix

 231      continue

!.......................................................................
!l    2.3.2 Apply the boundary conditions
!.......................................................................

          if (sbdry .ne. "periodic" .or. laddbnd.ne.0) &
            call wpbdry(i,k,ileft,iright,2)

!.......................................................................
!l    2.3.3 Solve the linear systems for each theta point
!l    If i>itl, combine equation for i and ii=iy-i+1, assuming
!l    perfect reflection at turning point, i.e. f(iyh)=f(iyh+1)
!.......................................................................

          ilslen=l_upper(i)-l_lower(i)+1
          if (sbdry .eq. "periodic" .and. l_upper(i).ne.setup0%ls) &
            ilslen=2*ilslen-1
          icombin=1
          if (sbdry.eq."periodic" .and. i.gt.itl_(1) .and. numindx.eq.2 &
            .and. laddbnd.eq.0) icombin=2

          do 233 j=2,jx

            do 2330 iic=1,2,icombin

!.......................................................................
!l    2.3.3.1 Arrange matrix according to band matrix solver used
!     Uses the fact that the matrix is tri-diagonal
!.......................................................................

              do 2331 jcolumn=1,iband
                do 2332 l=1,ilslen
                  ijmat=(l-1)*iband+jcolumn
                  zmat(ijmat)=bndmats(j,l,jcolumn,iic)
 2332           continue
 2331         continue

!     When combining eqs. for i and ii, one assumes that iband=5 with
!     only three non-zero element per row, at jcolumn=1, 3 and 5

              if (icombin .eq. 2) then
!     rhs:
                ilento=2*ilslen-2
                rhspar(ilslen-1,j,1)=rhspar(ilslen-1,j,1)+ &
                  rhspar(ilslen-1,j,2)
                rhspar(ilslen,j,1)=rhspar(ilslen,j,1)+rhspar(ilslen,j,2)
                do 2333 l=ilslen+1,ilento-2,2
                  rhspar(l,j,1)=rhspar(ilento-l,j,2)
                  rhspar(l+1,j,1)=rhspar(ilento-l+1,j,2)
 2333           continue
                rhspar(ilento,j,1)=rhspar(1,j,2)
!     matrix:
!     point ilslen-1 (i.e. l_upper(i)):
                ijdiag=(ilslen-2)*iband+3
                zmat(ijdiag)=zmat(ijdiag)+bndmats(j,ilslen-1,3,2)
                zmat(ijdiag+2)=bndmats(j,ilslen-1,1,2)
!     point ilslen (i.e. setup0%ls+2-l_upper(i)):
                ijdiag=(ilslen-1)*iband+3
                zmat(ijdiag)=zmat(ijdiag)+bndmats(j,ilslen,3,2)
                zmat(ijdiag+2)=bndmats(j,ilslen,1,2)
!     points ilslen-2 to 2:
                do 2334 jcolumn=1,iband,2
                  ijeff=(ilslen-2)*iband+(iband-jcolumn+1)
                  do 2335 l=ilslen-3,2,-2
                    ijeff=ijeff+2*iband
                    zmat(ijeff)=bndmats(j,l,jcolumn,2)
                    zmat(ijeff+iband)=bndmats(j,l+1,jcolumn,2)
 2335             continue
 2334           continue
!     new total number of unknown
                ilslen=2*ilslen-2
!     correct point l=3 (i.e. l=setup0%ls):
                ijdiag=(ilslen-2)*iband+3
                zmat(ijdiag+1)=zmat(ijdiag+2)
                zmat(ijdiag+2)=0.0
!     correct point l=2:
                ijdiag=(ilslen-3)*iband+3
                zmat(ijdiag+2)=bndmats(j,2,2,2)
!     point l=1:
                ijdiag=(ilslen-1)*iband+3
                zmat(ijdiag-2)=bndmats(j,1,4,2)
                zmat(ijdiag-1)=bndmats(j,1,5,2)
                zmat(ijdiag  )=bndmats(j,1,3,2)
              endif

!%OS
              call dcopy(ilslen*iband,zmat,1,zmat2,1)
              call dcopy(ilslen,rhspar(1:ilslen,j,iic),1,zxdumy2,1)
!%OS

!.......................................................................
!l    2.3.3.2 Call solver
!.......................................................................

              do 2336 l=ilslen+1,ilslen+iright
                rhspar(l,j,iic)=0.0
 2336         continue
              zerr=1.0e-08
              call nonsym(zmat,zxdumy,rhspar(1:ilslen,j,iic),ilslen,ileft &
                ,iright,zerr,icond)
              if (icond .ne. 0) write(6,'(/," WARNING: bad condition in" &
                ," nonsym: icond = ",i4)') icond

!%OS  check solution
              zdiff=0.0
              do 2337 l=1,ilslen
                zzzl=0.0
                do 2338 jcol=1,iband
                  zzzl=zzzl+zmat2((l-1)*iband+jcol) &
                    *rhspar(l-ileft-1+jcol,j,iic)
 2338           continue
                zdiff=zdiff + abs((zxdumy2(l)-zzzl)/zxdumy2(l))
 2337         continue
              zdiff=zdiff/(1.*ilslen)
!%OS  if (zdiff .gt. 1.0E-08) print *,' zdiff,i,j,n= ',zdiff,i,j,n
!%OS

!     redistribute solution for ii points onto rhspar(.,j,2) and redefine ilslen
              if (icombin .eq. 2) then
                ilslen=ilslen/2+1
                do 2339 l=ilslen-1,ilento-2,2
                  rhspar(ilento-l,j,2)=rhspar(l,j,1)
                  rhspar(ilento-l+1,j,2)=rhspar(l+1,j,1)
 2339           continue
                rhspar(1,j,2)=rhspar(ilento,j,1)
              endif

 2330       continue
 233      continue

!.......................................................................
!l    2.3.4 Copy solution onto fnp1
!.......................................................................

          do 2340 l=1,ilslen
            ileff=lsprtob(l)
            iieff=iymxper(ileff)+1-i
            do 2341 j=2,jx
              fnp1(i,j,k,ileff)=rhspar(l,j,1)
              fnp1(iieff,j,k,ileff)=rhspar(l,j,2)
 2341       continue
 2340     continue

!     end of loop over theta
 230    continue

!     end of loop over general species
 220  continue

      if (sbdry .eq. "periodic") then
        call dcopy(iyjx2*ngen,fnp1(0:iy+1,0:jx+1,1:ngen,1   ),1, &
                              fnp1(0:iy+1,0:jx+1,1:ngen,setup0%ls+1),1)
        call dcopy(iyjx2*ngen,fnp1(0:iy+1,0:jx+1,1:ngen,setup0%ls  ),1, &
                              fnp1(0:iy+1,0:jx+1,1:ngen,0   ),1)
      else
        call dcopy(iyjx2*ngen,fnp1(0:iy+1,0:jx+1,1:ngen,1   ),1, &
                              fnp1(0:iy+1,0:jx+1,1:ngen,0   ),1)
        call dcopy(iyjx2*ngen,fnp1(0:iy+1,0:jx+1,1:ngen,setup0%ls  ),1, &
                              fnp1(0:iy+1,0:jx+1,1:ngen,setup0%ls+1),1)
      endif

!.......................................................................
!l    3. Check the solution
!.......................................................................

      if (scheck .eq. "enabled") call wpcheck

!.......................................................................
!l    4. Redefine fnp1 to be unique at v=0
!.......................................................................

      if (nonadi .ne. 3) then

        do 400 k=1,ngen
          do 410 l=1,setup0%ls
            zs=0.
            zt=0.
            do 420 i=1,iymxper(l)
              zs=zs+cynt2(i,l)
              zt=zt+cynt2(i,l)*fnp1(i,1,k,l)
 420        continue
            do 430 i=1,iymxper(l)
              fnp1(i,1,k,l)=zt/zs
 430        continue
 410      continue

          if (sbdry .eq. "periodic") then
            call dcopy(iyp1+1,fnp1(0:iyp1,1,k,1   ),1, &
                              fnp1(0:iyp1,1,k,setup0%ls+1),1)
            call dcopy(iyp1+1,fnp1(0:iyp1,1,k,setup0%ls  ),1, &
                              fnp1(0:iyp1,1,k,0   ),1)
          else
            call dcopy(iyp1+1,fnp1(0:iyp1,1,k,1   ),1, &
                              fnp1(0:iyp1,1,k,0   ),1)
            call dcopy(iyp1+1,fnp1(0:iyp1,1,k,setup0%ls  ),1, &
                              fnp1(0:iyp1,1,k,setup0%ls+1),1)
          endif

 400    continue

      endif

      return
      end subroutine wptramu

end module wptramu_mod
