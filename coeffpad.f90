module coeffpad_mod

!
!

contains

      subroutine coeffpad(k)
      use param_mod
      use comm_mod
      implicit integer (i-n), real*8 (a-h,o-z)
!      real*8,dimension(iy):: prnt1,prnt2,prnt3
!      real*8,dimension(iy):: prnt4,prnt5,prnt6

!..................................................................
!     This routine adds in the collisional contribution to the
!     coefficients employed in time advancement..
!     scatfrac=0. disables pitch angle scattering (along with
!       mx=0, see cqlinput_help. scatfrac=1. by default).
!..................................................................


      do 10 j=1,jx
        do 11 i=1,iy
          da(i,j)=da(i,j)+cal(i,j,k,l_)
          db(i,j)=db(i,j)+cbl(i,j,k,l_)
          dc(i,j)=dc(i,j)+scatfrac*ccl(i,j,k,l_)
          dd(i,j)=dd(i,j)+cdl(i,j,k,l_)
          de(i,j)=de(i,j)+scatfrac*cel(i,j,k,l_)
          df(i,j)=df(i,j)+scatfrac*cfl(i,j,k,l_)
 11     continue
 10   continue
!      do i=1,iy
!         prnt1(i)=da(i,2)
!         prnt2(i)=db(i,2)
!         prnt3(i)=dc(i,2)
!         prnt4(i)=dd(i,2)
!         prnt5(i)=de(i,2)
!         prnt6(i)=df(i,2)
!      enddo


      return
      end
end module coeffpad_mod