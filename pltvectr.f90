module pltvectr_mod

  !---BEGIN USE

  use aminmx_mod, only : aminmx

  !---END USE

!
!

contains

      subroutine pltvectr(xt,yt,xh,yh,rheads,jpxy,ipxy,veclen,noplots)
      use r8subs_mod, only : rbound
      use aminmx_mod, only : aminmx
      implicit integer (i-n), real*8 (a-h,o-z)
      character*8 noplots
      save
!
!     Plots vector field in x(horizontal),y(vertical) space.
!     Tail of vector is at xt(1:jpxy,1:ipxy),yt(1:jpxy,1:ipxy).
!     Vectors have x,y components xy,yh. rheads(1:jpxy) is used
!       for calculation of vector lengths (for each y).
!       Vectors are scaled to length veclen.
!
      dimension tam1(2),tam2(2)
      dimension xt(jpxy,ipxy),yt(jpxy,ipxy)
      dimension xh(jpxy,ipxy),yh(jpxy,ipxy)
      dimension rheads(*)
!
      REAL*4 RPX1,RPY1,RPX2,RPY2
      REAL*4 ANGLE,BARB,SIZE
      !XXX REAL*4 RBOUND
!BH011228 Modifications for plotting with PGPLOT,  011228.

      ep90=1.d+90
      rscale=-ep90
      do 100 i=1,ipxy
         do 90 j=1,jpxy
            rheads(j)=sqrt(xh(j,i)**2+yh(j,i)**2)
 90      continue
         call aminmx(rheads,1,jpxy,1,dummy,rscalexx,kmin,kmax)
         rscale=max(rscale,rscalexx)
 100  continue
      if(rscale.eq.0.) rscale=1.
!
!     SET ARROW HEAD STYLE:
      ANGLE=90.
      BARB=0.7
      CALL PGSAH(1,ANGLE,BARB)
!     SAVE PGPLOT attributes, and reset character/arrowhead size:
      CALL PGSAVE
      SIZE=10./jpxy
      CALL PGSCH(SIZE)
!      write(*,*)''
!      write(*,*)'pltvectr: j,i,x1,y1,x2,y2='
      do 200 i=1,ipxy
         do 201 j=1,jpxy
!            if(xh(j,i).eq.0. .and. yh(j,i).eq.0.) goto 201
            tam1(1)=xt(j,i)
            tam2(1)=yt(j,i)
            tam1(2)=xt(j,i)+veclen*xh(j,i)/rscale
            tam2(2)=yt(j,i)+veclen*yh(j,i)/rscale
            RPX1=RBOUND(tam1(1))
            RPY1=max(RBOUND(tam2(1)),0.)
            RPX2=RBOUND(tam1(2))
            RPY2=max(RBOUND(tam2(2)),0.)
            if(abs(RPX2).lt. 1.e-20) RPX2=0.  !Helps near vpar=0.
!          write(*,5) j,i,RPX1,RPY1,RPX2,RPY2
! 5        format(2i5,4(1pe12.5))
!          if ((RPX1.ne.0.) .and. (RPY1.ne.0.) .and.
!     +        (RPX2.ne.0.) .and. (RPY2.ne.0.)) then
          if ((RPY1.ne.0.) .and. &
              (RPY2.ne.0.)) then
              CALL PGARRO(RPX1,RPY1,RPX2,RPY2)
           endif
 201    continue
 200  continue
!     RESTORE PGPLOT ATTRIBUTES
      CALL PGUNSA
!
      return
      end





end module pltvectr_mod
