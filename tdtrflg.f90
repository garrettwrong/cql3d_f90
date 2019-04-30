module tdtrflg_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use bcast_mod, only : ibcast

  !---END USE

!
!

contains

      subroutine tdtrflg
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..............................................................
!     This routine sets some mesh flag arrays to facilitate
!     coding of the transport routine
!..............................................................


      if (transp.eq."disabled") return

!..............................................................
!     define lpt(i); the transport mesh point l such that
!     the particle characterized by transport pitch angle index i
!     is passing at l-1 and trapped at l. If this never happens
!     set equal to lrors.
!..............................................................

      do 20 i=1,iytr(lrors)/2
         lpt(i)=lrors
         if (idx(i,lrors).lt.itl_(lrors)) go to 20
         do 30 l=lrors-1,2,-1
!         write(*,*)'tdtrflg:l,i,idx(i,l),itl_(l),idx(i,l-1),itl_(l-1)'
!     +           ,l,i,idx(i,l),itl_(l),idx(i,l-1),itl_(l-1)
         if (ipacktp.eq.3) then !Older system for soln_method=direct,tfac=1.
            if (idx(i,l).gt.itl_(l) .and. idx(i,l-1).lt. itl_(l-1)) then
               if (idx(i,l-1).ne.0) lpt(i)=l
!               write(*,*)'lpt(i),i,l',lpt(i),i,l
            endif
         elseif (ipacktp.eq.0) then
            if (idx(i,l).ge.itl_(l) .and. idx(i,l-1).lt. itl_(l-1)) then
               if (idx(i,l-1).ne.0) lpt(i)=l
!               write(*,*)'lpt(i),i,l',lpt(i),i,l
            endif
         endif
 30   continue
 20   continue

!..............................................................
!     Define l_lower(i); the lowest index l for which pitch angle
!     index i is defined on the local transport pitch angle mesh.
!..............................................................

      call ibcast(l_lower,1,iytr(lrors))
      do 50 l=2,lrors
        if (iytr(l-1).lt. iytr(l)) then
          do 80 i=iytr(l-1)/2+1,iytr(lrors)-iytr(l-1)/2
            l_lower(i)=l
            if (lpt(i).eq.l_lower(i)) lpt(i)=lrors
 80       continue
        endif
 50   continue
      do 100 i=1,iytr(lrors)/2
        ii=iytr(lrors)+1-i
        lpt(ii)=lpt(i)
 100  continue

      write(*,*)'tdtrflg:  iytr(1:lrors)=',iytr(1:lrors)
      write(*,*)'tdtrflg:  l_lower(1:iytr)=',l_lower(1:iytr(lrors))
      write(*,*)'tdtrflg:  lpt(1:iytr)=',lpt(1:iytr(lrors))
!.......................................................................
!     initialize coefficients for v to r velocity meshes transformation
!.......................................................................

      call bcast(f_vtor,one,jx*ngen*lrors*18)

      return
      end
end module tdtrflg_mod
