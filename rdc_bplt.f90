module rdc_bplt_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast
  use pltdf_mod, only : pltcont

  !---END USE


contains

      subroutine rdc_bplt(krf)
      use param_mod
      use comm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!..................................................................
!     plots rf rdcb coefficient as a contour plot on cql3d grid.
!..................................................................

!MPIINSERT_INCLUDE

      character*8 pltvlhb
      character*8 pltovlp

      data pltvlhb /'enabled'/
      data pltovlp /'enabled'/


!MPIINSERT_IF_RANK_NE_0_RETURN
 ! make plots on mpirank.eq.0 only

      if (noplots.eq."enabled1") return
      if (pltvlhb.ne."enabled") return
!$$$      if (pltovlp.eq."enabled".and. mrfn.gt.1) then
!cc      call bcast(temp1(1,0),zero,iy*(jx+1)) ! YuP-101215: error?
      call bcast(temp1(0:iyjx2-1,0),zero,iyjx2)  !temp1(0:iyp1,0:jxp1)

      write(*,*)'rdc_bplt(krf): mrfn =',mrfn,' krf=',krf

!$$   do 560 k=1,mrfn
      do 561 j=1,jx
         do 562 i=1,iy
            temp1(i,j)=rdcb(i,j,lr_,krf)
 562     continue
 561  continue
      CALL PGPAGE
      itype=7 ! means: plots are made for rdcb
      call pltcont(nrdcspecies(krf),1, &
           'Contours of RdcB vs. v_parallel,v_perp',7)
      write(t_,552) lr_
 552  format(" Flux surface number",i3,";   all modes, krf=",i2)
      CALL PGMTXT('B',10.,0.,0.,t_)

!$$$  560    continue

!$$$  endif


      return
      end


end module rdc_bplt_mod
