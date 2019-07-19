! Copyright 2019 Garrett Wright, Princeton Plasma Physics Laboratory,
!    contracted by the U.S. Department of Energy (putnumberhere).
!
! This file is part of cql3d_f90. See LICENSE.
!
! cql3d_f90 is free software: you can redistribute it and/or modify it
! under the terms of the GNU Affero General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! cql3d_f90 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with cql3d_f90.  If not, see <https://www.gnu.org/licenses/>.

module pltvectr_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use aminmx_mod, only : aminmx

  !---END USE

!
!

contains

      subroutine pltvectr(xt,yt,xh,yh,rheads,jpxy,ipxy,veclen,noplots)
      use r8subs_mod, only : rbound
      use aminmx_mod, only : aminmx
      implicit integer (i-n), real(c_double) (a-h,o-z)
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
      real(c_float) RPX1,RPY1,RPX2,RPY2
      real(c_float) ANGLE,BARB,SIZE
      !XXX real(c_float) RBOUND
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
#ifndef NOPGPLOT
      CALL PGSAH(1,ANGLE,BARB)
#endif
!     SAVE PGPLOT attributes, and reset character/arrowhead size:
#ifndef NOPGPLOT
      CALL PGSAVE
#endif
      SIZE=10./jpxy
#ifndef NOPGPLOT
      CALL PGSCH(SIZE)
#endif
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
#ifndef NOPGPLOT
              CALL PGARRO(RPX1,RPY1,RPX2,RPY2)
#endif
           endif
 201    continue
 200  continue
!     RESTORE PGPLOT ATTRIBUTES
#ifndef NOPGPLOT
      CALL PGUNSA
#endif
!
      return
      end subroutine pltvectr





end module pltvectr_mod
