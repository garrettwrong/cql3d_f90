!
!
module pltdf_mod

  use iso_c_binding, only : c_double
  integer, parameter, public :: nconta=100
  ! these are shared with other routines (outside the contains)
  real(c_double), public :: cont(nconta),tempcntr(nconta)
  ! PGLOCAL1 from pltcont.f
  real(c_double), pointer, public :: wx(:), wy(:)
  integer, public :: IIX, IIY, JXQ

  save

contains 

  subroutine pltdf
    use param_mod
    use cqcomm_mod
    use r8subs_mod, only : dcopy
    implicit integer (i-n), real*8 (a-h,o-z)

    !..................................................................
    !     if (pltd.eq."enabled" or pltd.eq."color") then
    !     subroutine pltdf contour plots the distribution function, f
    !     if (pltd.eq."df"  or "df_color") then
    !     subroutine pltdf also plots the difference between f
    !     at the current time and f at the previous time
    !     step.
    !..................................................................


    REAL RILIN

    if (noplots.eq."enabled1") return

    if (pltd.eq."disabled") return


    do 10 k=1,ngen

       ! This part is plotted for any pltd.ne.'disabled'

       call dcopy(iyjx2,f(0:iyjx2-1,0,k,l_),1,temp1(0:iyjx2-1,0),1)
       write(t_,550) k
550    format(1x,"Species ",i2," Distribution Function Contour Plot")
       !PGPLT        CALL PGPAGE
       itype=1 ! means: plots are made for distr.func f
       call pltcont(k,2,t_,itype) ! for f()
       write(t_,560)
560    format("Contour values:")
       RILIN=10.
       !PGPLT        CALL PGMTXT('B',RILIN,-.2,0.,t_)


       do 11 jcs=1,ncont,4
          write(t_,570) (tempcntr(jc),jc=jcs,min(jcs+3,ncont))
          if ((ncont/4)*4.ne.ncont .and. ncont-jcs.le.2) then
             icend=4 * 16 + 1
             t_(icend:icend)="$"
          endif
          RILIN=RILIN+1.
          !PGPLT          CALL PGMTXT('B',RILIN,-.2,0.,t_)
11     end do


       if (n.eq.0) goto 10
       
       ! Additionally, plot f(n+1)-f(n)
       if (pltd.eq."df" .or. pltd.eq."df_color")then
          
          do 20 i=1,iy
             do 21 j=1,jx
                temp1(i,j)=f(i,j,k,l_)-f_(i,j,k,l_)
21           end do
20        end do
          write(t_,530) k,n
530       format(1x, "Contours of df/dt for species",i2,1x,"during timestep",i5)
          !PGPLT        CALL PGPAGE
          itype=2 ! means: plots are made for df
          !PGPLT        call pltcont(k,1,t_,itype) ! for df
          RILIN=10.
          write(t_,560)
          !PGPLT        CALL PGMTXT('B',RILIN,-.2,0.,t_)
          
          do 12 jcs=1,ncont,4
             write(t_,570) (tempcntr(jc),jc=jcs,min(jcs+3,ncont))
             if ((ncont/4)*4.ne.ncont .and. ncont-jcs.le.2) then
                icend=4 * 16 + 1
                t_(icend:icend)="$"
             endif
             RILIN=RILIN+1.
             !PGPLT          CALL PGMTXT('B',RILIN,-.2,0.,t_)
12        end do

       endif !  pltd.eq."df" .or. pltd.eq."df_color"
       
       
10  end do

570 format(4(1pe16.6))
    
    return

!bug end do
  end subroutine pltdf

end module pltdf_mod
