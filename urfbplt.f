      subroutine urfbplt
      use param_mod
      use comm_mod
      implicit integer (i-n), real*8 (a-h,o-z)
      save

c..................................................................
c     plots rf urfb coefficient as a contour plot.
c..................................................................

CMPIINSERT_INCLUDE

      character*8 pltovlp

      data pltovlp /'enabled'/


CMPIINSERT_IF_RANK_NE_0_RETURN
 ! make plots on mpirank.eq.0 only

      if (noplots.eq."enabled1") return
cBH081105      iiplt3d=0
      iiplt3d=-1
      do i=1,nplota
         if (n.eq.nplt3d(i)) iiplt3d=n
      enddo
         
c      if (n/nplt3d*nplt3d.ne.n .and. n.ne.1) return
cBH081105      if (iiplt3d.eq.0 .and. n.ne.1) return
      if (iiplt3d.eq.-1) return
      if (mplot(l_).eq."disabled") return
      if (urfmod.eq."disabled") return
      if (plturfb.eq."disabled") return
            
      if (pltovlp.eq."enabled".and. mrfn.gt.1) then
      
c       This do 560 loop was overlapping the modes on one plot.
        !YuP[2016] instead of overlapping contour plots,
        !sum them and then plot contours for the total 
        !UrfB(all modes) at a given surface (and for each k-species).
        do k=1,ngen ! YuP[10-2016] scan general species: can be more than one
           !initialize for each species: one plot for each k
           call bcast(temp1(0,0),zero,iyjx2)  !temp1(0:iyp1,0:jxp1)
           do 560 krf=1,mrfn
             if (nrfspecies(krfn(krf)) .eq.k) then 
               !sum-up modes for a given species only
               do j=1,jx
               do i=1,iy
                 !temp1(i,j)=urfb(i,j,indxlr_,krf) ! YuP: original
                 temp1(i,j)= temp1(i,j)+urfb(i,j,indxlr_,krf) ! YuP: new version: sum-up
               enddo
               enddo
               !call pltcont(1,1,'Contours of UrfB vs. v_parallel,v_perp')  !YuP:original
             endif
 560       continue ! krf mode (usually = harmonic number)
           if( MAXVAL(temp1)-MINVAL(temp1) .gt. 0.d0 ) then
             CALL PGPAGE ! new page for each k 
	       itype=4 ! means: plots are made for urfb
             call pltcont(k,1,'Contours of UrfB vs. v_parallel,v_perp',
     +         itype) !YuP:summed-up
             write(t_,552) 
             CALL PGMTXT('B',10.,0.,0.,t_)
             write(t_,553) lr_
             CALL PGMTXT('B',11.,0.,0.,t_)
             write(t_,692) MAXVAL(temp1) !YuP[10-2016] max value for this krf
             CALL PGMTXT('B',12.,0.,0.,t_)
             write(t_,693) k 
             CALL PGMTXT('B',13.,0.,0.,t_)
           endif
        enddo ! k species
        
      endif
      
c     This do 680 loop plots the individual mode contributions:
      do k=1,ngen ! YuP[10-2016] scan general species
      do 680 krf=1,mrfn
         do j=1,jx
            do i=1,iy
               temp1(i,j)=urfb(i,j,indxlr_,krf)
            enddo
         enddo
         if (nrfspecies(krfn(krf)) .eq. k) then
         if( MAXVAL(temp1)-MINVAL(temp1) .gt. 0.d0 ) then
          CALL PGPAGE ! opens new page for each krf-mode
          itype=4 ! means: plots are made for urfb
          call pltcont(k,1,'Contours of UrfB vs. v_parallel,v_perp',
     +      itype)
          write(t_,690) 
          CALL PGMTXT('B',10.,0.,0.,t_)
          ! write flux surface number and mode number;
          ! also harmonic number and species number (added YuP[10-2016])
          write(t_,691) lr_ ,krf,nharm(krf),k 
          CALL PGMTXT('B',11.,0.,0.,t_)
          write(t_,692) MAXVAL(temp1) !YuP[10-2016] max value for this krf
          CALL PGMTXT('B',12.,0.,0.,t_)
         endif
         endif
 680  continue ! krf
      enddo ! k species

 552  format("Contours of the rf (v,v) diffusion coefficient, urfb")
 553  format(" Flux surface number",1x,i3,"; all modes")
 690  format("Contours of the rf (v,v) diffusion coefficient, urfb")
 691  format("Flux surf.N",i3,";  mode,nharm=",2i5,";  Species k=",i1)
 692  format("Max value for this surface/mode:",e13.3)
 693  format("Species k=",i1)
      return
      end
