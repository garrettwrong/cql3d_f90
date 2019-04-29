module sourc0_mod

  !---BEGIN USE

  !---END USE

!
!

contains

      subroutine sourc0
      use param_mod
      use comm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!..................................................................
!     define source uniquely at v=0
!..................................................................


      do 10 k=1,ngen
        s=0.
        u=0.
        do 11 i=1,iy
          u=u+source(i,1,k,indxlr_)*cynt2(i,l_)*vptb(i,lr_)
          s=s+cynt2(i,l_)*vptb(i,lr_)
 11     continue
        do 12 i=1,iy
          source(i,1,k,indxlr_)=u/s
 12     continue
 10   continue
      return
      end
end module sourc0_mod
