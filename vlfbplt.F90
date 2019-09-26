! Copyright 2019 Garrett Wright, Princeton Plasma Physics Laboratory,
!    contracted by the U.S. Department of Energy (DE-AC02-09CH11466).
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

module vlfbplt_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use pltdf_mod, only : pltcont
  use urfwrong_mod, only : urfwrong

  !---END USE


contains

  subroutine vlfbplt
    use cqlconf_mod, only : setup0
      use param_mod
      use cqlcomm_mod
      use r8subs_mod, only : luf
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     plots rf cqlb coefficient as a contour plot.
!..................................................................

#ifdef __MPI
      include 'cql3d_mpilib.h'
#endif

      character*8 pltvlfb
      character*8 pltovlp
      real(c_float) RTAB1(iy),RTAB2(iy) ! local, for PGPLOT

      save pltvlfb,pltovlp
      data pltvlfb /'enabled'/
      data pltovlp /'enabled'/


#ifdef __MPI
      if(mpirank.ne.0) return
#endif
 ! make plots on mpirank.eq.0 only

      if (setup0%noplots.eq."enabled1") return
      if (pltvlfb.ne."enabled") return

      if (pltovlp.eq."enabled".and. mrfn.gt.1) then

      if (setup0%cqlpmod.ne."enabled") then
         call bcast(temp1(0:iy+1,0:jx+1),zero,iyjx2)

         do 560 k=1,mrfn
            do 561 j=1,jx
               do 562 i=1,iy
                  temp1(i,j)=cqlb(i,j,indxlr_,k)
 562           continue
 561        continue
#ifndef NOPGPLOT
            CALL PGPAGE
#endif
            itype=5 ! means: plots are made for vlfb
            call pltcont(1,1,'Contours of CqlB vs. v_parallel,v_perp', &
                           itype)
            write(t_,552) lr_
 552        format("Flux surface number",i3,"; all modes")
#ifndef NOPGPLOT
            CALL PGMTXT('B',10.,0.,0.,t_)
#endif
 560     continue

      elseif (setup0%cqlpmod.eq."enabled") then

!       Save l_ (just in case), and pass l_=l to pltcont.
!       Restore at end of loop.
        l_save=l_
        do l=1,lz
           l_=l

!cc           call bcast(temp1(1,0),zero,iy*(jx+1)) ! YuP-101215: error?
           call bcast(temp1(0:iy+1,0:jx+1),zero,iyjx2)  !temp1(0:iyp1,0:jxp1)

        do 570 k=1,mrfn
          do 571 j=1,jx
            do 572 i=1,iy
              temp1(i,j)=wcqlb(i,j,k,l)
 572        continue
 571      continue
#ifndef NOPGPLOT
        CALL PGPAGE
#endif
          itype=5 ! means: plots are made for urfb
          call pltcont(1,1,'Contours of CqlB vs. v_parallel,v_perp', &
                       itype)
!$$$          call gstxno(80.)
          write(t_,552) lr_
#ifndef NOPGPLOT
          CALL PGMTXT('B',10.,0.,0.,t_)
#endif
 570    continue

        enddo
        l_=l_save

      endif ! setup0%cqlpmod selection

      endif !  mrfn.gt.1


      ! Next part - for both mrfn=1 and mrfn>1
      if (setup0%cqlpmod.ne."enabled") then


      do 680 k=1,mrfn

!..................................................................
!     Compute vpar21/vte and vpar11/vte normalized velocities
!     following vlf.f (or urfpack).
!..................................................................
        rr=solrz(1,lr_)
        if (vlfnpvar.eq."constant") then
          vlfnpar=vlfnp(k)
          vlfdnpar=vlfdnp(k)
          vlfddnpa=vlfddnp(k)
          vlfnper=vlfnperp(k)
        elseif (vlfnpvar.eq."1/R") then
          vlfnpar=vlfnp(k)*rpcon(lr_)/rr
          vlfdnpar=vlfdnp(k)*rpcon(lr_)/rr
          vlfddnpa=vlfddnp(k)*rpcon(lr_)/rr
          vlfnper=vlfnperp(k)*rpcon(lr_)/rr
        else
          stop 'error in vlf input'
        endif
        ks=1
        bcnst=abs(bnumb(ks))*charge/(fmass(ks)*clight)
        wx=0.5
        wce=bcnst*bmidplne(lr_)

        signn=1.
        if (vlfnp(k).lt.0.) signn=-1.
        omomn=1.
        signom=1.

        if (nharm(k).eq.0 .or. kiong(1).eq.1) then
          if (kiong(1).eq.1) then
            omn=nharm(k)*bcnst*bmidplne(lr_)/omega(k)
            omomn=1.0-omn
!990131            signom=sign(1.0,omomn)
            signom=sign(one,omomn)
          endif

          vparu=clight*signom*omomn/ &
            (signn*vlfnpar-wx*(vlfdnpar+vlfddnpa))
          vparl=clight*signom*omomn/ &
            (signn*vlfnpar+wx*(vlfdnpar+vlfddnpa))

          vpar21dv=vparl/vth(1,lr_)
          vpar11dv=vparu/vth(1,lr_)

        elseif (nharm(k).gt.0) then

          rnpar1=signn*vlfnpar-wx*(vlfdnpar+vlfddnpa)
          rnpar2=signn*vlfnpar+wx*(vlfdnpar+vlfddnpa)
          r1mn12=1.-rnpar1*rnpar1
          r1mn22=1.-rnpar2*rnpar2
          omn=nharm(k)*bcnst*bmidplne(lr_)/omega(k)
          omn2=omn*omn
          rad1=omn2-r1mn12
          rad2=omn2-r1mn22

!..................................................................
!     Ray outside of pinch point. No reson. (omega(k) will be .gt. omega_ce).
!..................................................................

          if (rad1.le.0.0.and.rad2.le.0.0)  then

!..................................................................
!     set flag to skip this ray element
!..................................................................

            upar21=0.0
            upar11=0.0
            go to 20
          endif

!..................................................................
!     The following case should not be possible:
!..................................................................

          if (rad1.gt.0.0.and.rad2.lt.0.0)  call urfwrong(6)

          if (rad2.gt.0.0.and.rad1.le.0.0)  then
            upar01=0.0
            uu01=0.0
          else
            upar01=clight*omn*rnpar1/(r1mn12*vnorm)
            uu01=clight*sqrt(rad1)/(r1mn12*vnorm)
          endif
          upar02=clight*omn*rnpar2/(r1mn22*vnorm)
          uu02=clight*sqrt(rad2)/(r1mn22*vnorm)

!..................................................................
!     upar11 and upar12 are parallel momentum-per-mass limits of the
!     resonance ellipse associated with rnpar1..., upar12 and upar22,
!     associated with rnpar2.
!..................................................................

          upar11=upar01-uu01
          upar12=upar01+uu01
          upar21=upar02-uu02
          upar22=upar02+uu02

 20       continue

          j1=luf(abs(upar21),x,jx)
          j2=luf(abs(upar11),x,jx)
          vpar21dv=upar21*vnorm/gamma(j1)/vth(1,lr_)
          vpar11dv=upar11*vnorm/gamma(j2)/vth(1,lr_)

        endif

        do  j=1,jx
           do  i=1,iy
              temp1(i,j)=cqlb(i,j,indxlr_,k)
!              if (mod(i,1).eq.0 .and. mod(j,1).eq.0) then
!                 write(*,*) 'i,j,cqlb(i,j,indxlr_,k)',i,j,cqlb(i,j,k)
!              endif
           enddo
        enddo

#ifndef NOPGPLOT
        CALL PGPAGE
#endif
        itype=5 ! means: plots are made for vlfb
        call pltcont(1,1,'Contours of CqlB vs. v_parallel,v_perp',itype)
!$$$        call gstxno(80.)
        write(t_,660) lr_,k
#ifndef NOPGPLOT
        CALL PGMTXT('B',10.,0.,0.,t_)
#endif
        write(t_,661) vpar21dv,vpar11dv
#ifndef NOPGPLOT
        CALL PGMTXT('B',11.,0.,0.,t_)
#endif

 680  continue

!..................................................................
!
!     setup0%cqlpmod.eq."enabled" case
!
!..................................................................

      elseif (setup0%cqlpmod.eq."enabled") then

!     Save l_ (just in case), and pass l_=l to pltcont.
!     Restore at end of loop.
      l_save=l_
      do l=1,lz

      do 780 k=1,mrfn

!..................................................................
!     Compute vpar21/vte and vpar11/vte normalized velocities
!     following vlf.f (or urfpack).
!..................................................................
        rr=solrz(1,lr_)
        if (vlfnpvar.eq."constant") then
          vlfnpar=vlfnp(k)
          vlfdnpar=vlfdnp(k)
          vlfddnpa=vlfddnp(k)
          vlfnper=vlfnperp(k)
        elseif (vlfnpvar.eq."1/R") then
          vlfnpar=vlfnp(k)*rpcon(lr_)/rr
          vlfdnpar=vlfdnp(k)*rpcon(lr_)/rr
          vlfddnpa=vlfddnp(k)*rpcon(lr_)/rr
          vlfnper=vlfnperp(k)*rpcon(lr_)/rr
        else
          stop 'error in vlf input'
        endif
        ks=1
        bcnst=abs(bnumb(ks))*charge/(fmass(ks)*clight)
        wx=0.5
        wce=bcnst*bmidplne(lr_)

        signn=1.
        if (vlfnp(k).lt.0.) signn=-1.
        omomn=1.
        signom=1.

        if (nharm(k).eq.0 .or. kiong(1).eq.1) then
          if (kiong(1).eq.1) then
            omn=nharm(k)*bcnst*bmidplne(lr_)/omega(k)
            omomn=1.0-omn
!990131            signom=sign(1.0,omomn)
            signom=sign(one,omomn)
          endif

          vparu=clight*signom*omomn/ &
            (signn*vlfnpar-wx*(vlfdnpar+vlfddnpa))
          vparl=clight*signom*omomn/ &
            (signn*vlfnpar+wx*(vlfdnpar+vlfddnpa))

          vpar21dv=vparl/vth(1,lr_)
          vpar11dv=vparu/vth(1,lr_)

        elseif (nharm(k).gt.0) then

          rnpar1=signn*vlfnpar-wx*(vlfdnpar+vlfddnpa)
          rnpar2=signn*vlfnpar+wx*(vlfdnpar+vlfddnpa)
          r1mn12=1.-rnpar1*rnpar1
          r1mn22=1.-rnpar2*rnpar2
          omn=nharm(k)*bcnst*bmidplne(lr_)/omega(k)
          omn2=omn*omn
          rad1=omn2-r1mn12
          rad2=omn2-r1mn22

!..................................................................
!     Ray outside of pinch point. No reson. (omega(k) will be .gt. omega_ce).
!..................................................................

          if (rad1.le.0.0.and.rad2.le.0.0)  then

!..................................................................
!     set flag to skip this ray element
!..................................................................

            upar21=0.0
            upar11=0.0
            go to 30
          endif

!..................................................................
!     The following case should not be possible:
!..................................................................

          if (rad1.gt.0.0.and.rad2.lt.0.0)  call urfwrong(6)

          if (rad2.gt.0.0.and.rad1.le.0.0)  then
            upar01=0.0
            uu01=0.0
          else
            upar01=clight*omn*rnpar1/(r1mn12*vnorm)
            uu01=clight*sqrt(rad1)/(r1mn12*vnorm)
          endif
          upar02=clight*omn*rnpar2/(r1mn22*vnorm)
          uu02=clight*sqrt(rad2)/(r1mn22*vnorm)

!..................................................................
!     upar11 and upar12 are parallel momentum-per-mass limits of the
!     resonance ellipse associated with rnpar1..., upar21 and upar22,
!     associated with rnpar2.
!..................................................................

          upar11=upar01-uu01
          upar12=upar01+uu01
          upar21=upar02-uu02
          upar22=upar02+uu02

 30       continue

          j1=luf(abs(upar21),x,jx)
          j2=luf(abs(upar11),x,jx)
          vpar21dv=upar21*vnorm/gamma(j1)/vth(1,lr_)
          vpar11dv=upar11*vnorm/gamma(j2)/vth(1,lr_)

        endif


        do 761 j=1,jx
          do 762 i=1,iy
            temp1(i,j)=wcqlb(i,j,k,l)
 762      continue
 761    continue
#ifndef NOPGPLOT
        CALL PGPAGE
#endif
        itype=5 ! means: plots are made for vlfb
        call pltcont(1,1,'Contours of CqlB vs. v_parallel,v_perp',itype)
!$$$        call gstxno(80.)

        write(t_,660) lr_,k
#ifndef NOPGPLOT
        CALL PGMTXT('B',10.,0.,0.,t_)
#endif
        write(t_,661) vpar21dv,vpar11dv
#ifndef NOPGPLOT
        CALL PGMTXT('B',11.,0.,0.,t_)
#endif

 780  continue

      enddo
      l_=l_save

      endif


 660  format("Flux surface number",1x,i3,"   mode=",i1)
 661  format("vpar21/vth=",1pe15.7,"   vpar11/vth=",1pe15.7)


      return
      end subroutine vlfbplt


end module vlfbplt_mod
