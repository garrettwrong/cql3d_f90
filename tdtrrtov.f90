module tdtrrtov_mod

!
!

contains

      subroutine tdtrrtov(f1)
      use param_mod
      use comm_mod
      implicit integer (i-n), real*8 (a-h,o-z)

!..............................................................
!     This routine interpolates onto the full velocity mesh
!     from the transport velocity mesh.
!     f1 ======> f1,    THAT IS, just fix up itl+/-1,itu+/-1
!..............................................................

      dimension f1(0:iyp1,0:jxp1,ngen,0:*),denrad(ngen,lrorsa)
      dimension denvel(ngen,lrorsa)

      call tdtrchkd(f1,vpint_,denrad)

!BH070419:   removing special itl,itu treatment for ipacktp=0
        if (ipacktp.eq.3) then

      do 10 k=1,ngen
        do 11 l=1,lrors
          ilr=lrindx(l)
          itl=itl_(l)
          itu=itu_(l)
          do 12 j=1,jx
            fact1=(vpint_(itl-2,ilr)-vpint(itl-2,ilr))*f1(itl-2,j,k,l) &
              +2.*(vpint_(itl+2,ilr)-vpint(itl+2,ilr))*f1(itl+2,j,k,l) &
              +(vpint_(itu+2,ilr)-vpint(itu+2,ilr))*f1(itu+2,j,k,l)
            fact2=vpint(itl-1,ilr)*f_lm(j,k,l)+ &
              2.*vpint(itl,ilr)+2.*vpint(itl+1,ilr)* &
              f_lp(j,k,l)+vpint(itu+1,ilr)*f_up(j,k,l)
            f1(itl,j,k,l)=fact1/fact2
            f1(itu,j,k,l)=f1(itl,j,k,l)
            f1(itl-1,j,k,l)=f_lm(j,k,l)*f1(itl,j,k,l)
            f1(itu+1,j,k,l)=f_up(j,k,l)*f1(itl,j,k,l)
            f1(itl+1,j,k,l)=f_lp(j,k,l)*f1(itl,j,k,l)
            f1(itu-1,j,k,l)=f1(itl+1,j,k,l)
 12       continue
 11     continue
 10   continue

        elseif (ipacktp.ne.0) then
           write(*,*)'STOP in tdtrrtov:  Check ipacktp'
           stop
        endif

      call tdtrchkd(f1,vpint,denvel)
      return
      end
end module tdtrrtov_mod