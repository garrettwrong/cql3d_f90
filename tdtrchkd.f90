module tdtrchkd_mod

  !---BEGIN USE

  use bcast_mod, only : bcast

  !---END USE

!
!

contains

      subroutine tdtrchkd(f1,vp,densty)
      use param_mod
      use comm_mod
      implicit integer (i-n), real*8 (a-h,o-z)

!.............................................................
!     This routine computes the density of species f1, given
!     that vpint=vp(lrindx). This will work for arrays defined on
!     radial as well as the velocity mesh.
!.............................................................

      dimension f1(0:iyp1,0:jxp1,ngen,0:*),vp(iy,lrza)
      dimension densty(ngen,lrza)
      if (iactst.eq."disabled") return

      do 10 k=1,ngen
        do 15 l=1,lrors
          call bcast(tam1,zero,jx)
          do 20 i=1,iy_(l)
            do 30 j=1,jx
              tam1(j)=tam1(j)+f1(i,j,k,l)*vp(i,lrindx(l))
 30         continue
 20       continue
          densty(k,l)=0.
          do 40 j=1,jx
            densty(k,l)=densty(k,l)+tam1(j)*cint2(j)
 40       continue
 15     continue
 10   continue
      return
      end
end module tdtrchkd_mod
