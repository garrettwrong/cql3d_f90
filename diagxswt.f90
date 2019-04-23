module diagxswt_mod

  !---BEGIN USE

  use diagdens_mod, only : diagdens
  use diagentr_mod, only : diagentr
  use r8subs_mod, only : dcopy

  !---END USE

!
!

contains

      subroutine diagxswt(k)
      use param_mod
      use comm_mod
      use r8subs_mod, only : dcopy
      implicit integer (i-n), real*8 (a-h,o-z)
      save

!..................................................................
!     Routine computes density gain due to particle sources and losses
!     during the during theta sweep in exsweep.
!     gains are positive - losses are negative.
!     2- density gained due to setting negative values of
!     distribution function to 0 during theta sweep.
!     3- particle source contribution.
!..................................................................



!..................................................................
!     Add in the source contribution.
!..................................................................

      sgain(3,k)=xlncur(k,lr_)*.5*dtr+sgain(3,k)
      call dcopy(iyjx2,temp2(0:iyjx2-1,0),1,temp1(0:iyjx2-1,0),1)
      s=0.
      if (ineg .eq. "disabled") go to 350

!..................................................................
!     if desired set negative values of distribution function = 0.
!..................................................................

      do 300 j=1,jx
        do 301 i=1,iy
          if(temp2(i,j) .lt. 0.) then
            temp1(i,j)=zero
            temp4(i,j)=-temp2(i,j)
          else ! temp2(i,j) .ge. 0.
            temp1(i,j)=temp2(i,j)
            temp4(i,j)=zero
          endif
 301    continue
 300  continue
      call diagdens(xline,xmidp,eline)
      engain(k)=engain(k)+eline*one_
      sgain(2,k)=xline*one_
 350  continue
      call dcopy(iyjx2,temp1(0:iyjx2-1,0),1,temp4(0:iyjx2-1,0),1)

!..................................................................
!     Compute power from df/dt and from setting neg f to zero.
!..................................................................

      if (n .gt. 0 .and. n/nchec*nchec .eq. n) then
        call diagentr(9,k)
        call diagentr(10,k)
        call dcopy(iyjx2,temp4(0:iyjx2-1,0),1,temp1(0:iyjx2-1,0),1)
      endif
      if (iactst.eq."disabled") go to 500

!..................................................................
!     if debugging, compute density at end of theta split
!..................................................................

      call diagdens(yline,ymidd,eline)
      yline=yline*one_
 500  continue
      return
      end
end module diagxswt_mod
