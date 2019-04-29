module soup0_mod

  !---BEGIN USE

  !---END USE

!
!

contains

      subroutine soup0
      use param_mod
      use comm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

      if (soucoord.eq."cart") return
      do 10 kk=1,ngen
        do 11 m=1,nso
          do 12 j=1,jx
            tam1(j)=-(x(j)-xem1(kk,m,lr_))**2/xem2(kk,m,lr_)
            sovt(j,kk,m,lr_)=exp(tam1(j))
 12       continue
 11     continue
 10   continue
      return
      end
end module soup0_mod
