module coefwti_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use coefmidt_mod, only : coefmidt

  !---END USE

!
!

contains

      subroutine coefwti(k)
      use param_mod
      use comm_mod
      use advnce_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     This routine computes w = dy(i,l_)*dd(i,j)/df(i,j) and then sets
!     g(w)=di(i,j,k,l_)=1./w - (1./(exp(w)-1.))
!     Note g(minus infinity)=1.; g(0)=.5; g(plus infinity)=0.
!     This is the Chang-Cooper algorithm (2-d form ala Karney)
!     This routine actually has modified the above by subtracting
!     off the RF contribution to the diffusion (see the coding).
!     This procedure, while not yet
!     justified theoretically, helps to keep the distribution from
!     going negative for strongly driven problems.
!..................................................................



!-YuP      call bcast(di(0,1,k,l_),half,(iy+1)*jx) ! could it be error?
      call bcast(di(0:(iy+1)*(jx+2)-1,0,k,l_),half,(iy+1)*(jx+2))
      !Note: di(0:iy,0:jx+1,1:ngen,lrors)

      if (chang .ne. "disabled") then

        if (chang.eq."noneg") then
          if (n.eq.1) jchang(k,l_)=jx
          do 100 j=1,jx
            do 1001 i=1,iy
              if (f(i,j,k,l_).le.0.d0) then
                jval=j
                go to 101
              endif
 1001       continue
 100      continue
          jval=jx
 101      continue
          if (jval.ne.jx) then
            if (jchang(k,l_).gt.jval-5) then
              jchang(k,l_)=jval-5
              if (jchang(k,l_).le.0) jchang(k,l_)=1
            endif
          endif
        elseif (chang.eq."enabled") then
          jchang(k,l_)=jx
        endif

       if (jchang(k,l_).le.0) jchang(k,l_)=1
       !-YuP 101126: added to prevent jchang=0

        op=one+em8
        call coefmidt(dff,3)
        do 10 j=1,jchang(k,l_)
          do 2 i=1,iy
            temc1(i)=dy(i,l_)*dd(i,j)*op*df(i,j)/(op*df(i,j) &
              -dff(i,j))**2
 2        continue

!...............................................................
!     Keep code from blowing up in the 9 loop below.
!..............................................................

          do 8 i=1,iy
            if(abs(temc1(i)).gt.em6) then
              temc3(i)=temc1(i)
            else
              temc3(i)=em6
            endif
 8        continue

          do 4 i=1,iy
            if(temc3(i).lt.sevenhun) then
              temc2(i)=temc3(i)
            else
              temc2(i)=sevenhun
            endif
 4        continue

          do 7 i=1,iy
            if(temc2(i).le.-sevenhun) then
              temc2(i)=-sevenhun
            endif
 7        continue

!..............................................................
!     Evaluate the ratio (modified)
!..............................................................

          do 9 i=1,iy
            di(i,j,k,l_)=1.d0/temc2(i)-(1.d0/(exp(temc2(i))-1.d0))
 9        continue

!.............................................................
!     Correct for errors in 9 above
!.............................................................

          do 21 i=1,iy
            wsub=(three+temc1(i))/(two+temc1(i))/three
            if(temc3(i).eq.em6) then
              di(i,j,k,l_)=wsub
            endif
 21       continue

!..............................................................
!     Limit for large positive or negative temc1 follows
!..............................................................

          do 22 i=1,iy
            if (temc2(i).eq.sevenhun) then
              di(i,j,k,l_)=1.d0/temc1(i)
            elseif (temc2(i).eq.-sevenhun) then
              di(i,j,k,l_)=1.d0+1.d0/temc1(i)
            endif
 22       continue
 10     continue

!...............................................................
!     Now force one-sided differencing if desired.
!...............................................................
       if (jchang(k,l_).le.0) jchang(k,l_)=1
       !-YuP 101126: added to prevent jchang=0

        if (jchang(k,l_).lt.jx) then
          do 60 j=jchang(k,l_),jx
            do 50 i=1,iy
              if(dd(i,j).ge.0.d0) then
                di(i,j,k,l_)=zero
              else
                di(i,j,k,l_)=one
              endif
 50         continue
 60       continue
        endif

      endif

!.......................................................................
!     Ensures correct differentiation at end of intervals
!.......................................................................

      do 3 j=1,jx
        di(0,j,k,l_)=0.d0
        di(iy,j,k,l_)=1.d0
 3    continue


      return
      end
end module coefwti_mod
