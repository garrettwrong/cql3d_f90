module pltlosc_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use pltdf_mod, only : pltcont

  !---END USE

!
!

contains

      subroutine pltlosc
      use param_mod
      use comm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
!
!     Plot contours of the loss region..
!
!     Modified from Graflib to pgplot calls by Yuri Petrov, 090727,
!     using PGPLOT + GRAFLIBtoPGPLOT.f routines (put in pltmain.f).
!

      if (noplots.eq."enabled1") return

      do 100 k=1,ngen
        suu=0.
        do 92001 i=1,iy
          do 92002 j=1,jx
            if(gone(i,j,k,indxlr_).lt.-.9) then
              temp1(i,j)=vnorm*x(j)*f(i,j,k,l_)/tau(i,lr_)
            elseif (gone(i,j,k,indxlr_).gt.0.) then
              temp1(i,j)=0.
            elseif (gone(i,j,k,indxlr_) .le. 0.) then
              temp1(i,j)=vnorm*x(j)*f(i,j,k,l_) &
                *(-gone(i,j,k,indxlr_))/tau(i,lr_)
            else
              temp1(i,j)=gone(i,j,k,indxlr_)*f(i,j,k,l_)
            endif
            if (temp1(i,j).ne.zero) suu=temp1(i,j)
92002     continue
92001   continue
        if (suu.eq.0.) go to 92003

        write(t_,588) k
 588    format("Loss due to lossmode(k) and torloss(k), k=",i5)
        CALL PGPAGE
        call pltcont(k,1,t_,8) ! itype=8 for pltlosc
        !call GSCPVS(.5,.4)
!$$$        call gxglfr(0)
!$$$        call gscpvs(.15,.85)
!$$$        write(t_,560)
!$$$        call gptx2d(t_)
!$$$ 560    format("Contour values:")
!$$$        write(t_,570) (temp2(jc,1),jc=1,(ncont/2)*2)
!$$$        call gptx2d(t_)
!$$$ 570    format((1x,e16.6,5x,e16.6),"$")
92003   continue
 100  continue
      return
      end
end module pltlosc_mod
