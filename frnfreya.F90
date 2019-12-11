!
!
      subroutine frnfreya(frmod_,fr_gyro_,beamplse_,beampon_,beampoff_, &
           hibrz_,mfm1_,noplots)
      use cqlconf_mod, only : setup0
      use param_mod
      use frplteq_mod, only :frplteq
      implicit none

!.................................................................
      include 'frcomm.h'
!     ONETWO DIVERGENCE: SEE COMMENTS AT BEGINNING OF FREYA
      character*8 frmod_,fr_gyro_,beamplse_,noplots,codeid
      real(c_double) :: beampon_,beampoff_
      real(c_double), intent(out) :: hibrz_(kz,ke,kb)
      integer, intent(out) :: mfm1_
      
      integer ib,ie,i ! local
      integer ipts ! arg. in freya(), to be found
      integer :: mi,mj ! arg. in frstup()      
      real(c_double) :: rin,rmax,zmin,zmax,zax,zshift !arg. in frstup()
      real(c_double) :: curdep ! local

      character*8 ifirst
      save ifirst
      data ifirst/"first"/

!     To pass these freya namelist in frcomm.h to comm.h
      frmod_=frmod
      fr_gyro_=fr_gyro
      beamplse_=beamplse
      beampon_=beampon
      beampoff_=beampoff


!..................................................................
!     This routine controls all of the NFREYA routines.
!     Called from subroutine tdintl and tdchief
!..................................................................

!..................................................................
!     Return if input variable frmod is "disabled"
!..................................................................

      if (frmod.ne."enabled") return
      !elong=0 ! YuP: not used here?

!..................................................................
!     First initialize some (iteration or time dependent) data.
!     [Possible dependence is on equilibrium, average temp,
!      density, etc.   Check in frstup.]
!..................................................................

      call frstup(mf,mfm1,mi,mj,nion,potsid,codeid,rin,rmax, &
        zax,zmin,zmax,zni,zne,zte,zti,zzi,p,psivol,xxx,yyy, &
        nprim,nimp,zeffctv,zshift)

!.......................................................................
!     zshift if frstup above is passed back from comm.h data to
!     adjust beam pivot height according to any shift in the equilibrim.
!.......................................................................

      if (ifirst.eq."first") then
         do ib=1,nbeams
            zpivot(ib)=zpivot(ib)-zshift
         enddo
         ifirst="notfirst"
      endif

!..................................................................
!     Call NFREYA (cray32 from ONETWO)
!..................................................................

      if(setup0%verbose>0) write(*,*)'frnfreya:mi,mj,codeid',mi,mj,codeid
      if(setup0%verbose>0) write(*,*)'frnfreya:rin,rmax',rin,rmax
      if(setup0%verbose>0) write(*,*)'frnfreya:zax,zmin,zmax',zax,zmin,zmax
      call freya(ipts,mi,mj,codeid,rin,rmax,zax,zmin,zmax)
      if(setup0%verbose>0) write(*,*) 'Done calling freya...'
      if (ipts.eq.0)write(*,*)'frnfreya: WARNING, ipts=0, NB missed plasma?'

!..................................................................
!     Compute the total number of particles/sec actually deposited in
!     the plasma. This requires subtracting off from the neutral current
!     the fraction that are lost at the aperture and the fraction
!     lost at the walls. ORBIT EFFECTS ARE NOT CURRENTLY CONSIDERED.
!     [2014-5: cql3d-fow now accounts for gc orbits and gyro-radius
!      offsets.]
!     Note that in the event that Freya input variable bptor is
!     used, we are specifying the total power injected into the tokamak
!     disregarding any losses at the apertures. The expression for
!     the deposited current below takes this into account. While fap may
!     not be zero, in routine FREYA the total number of neutralized
!     particles has been normalized up to account for losses at
!     the apertures. In the event that bptor is not used, bneut is
!     not renormalized, and the expression below makes sense as it
!     stands.
!..................................................................

      curdep=0.
      do 10 ib=1,nbeams
        do 20 ie=1,3
          curdep=curdep+bneut(ie,ib)*(1.-fap(ie,ib))*(1.-fwall(ie,ib))
 20     continue
 10   continue

!..................................................................
!     Now determine the source arrays used in CQL3D
!..................................................................

      call freyasou(xpts,ypts,zpts,rpts,vx,vy,vz,ipts,curdep, &
        bmsprd,multiply,multiplyn)


!..................................................................
!     Pass additional data to cql3d through subroutine arguments
!..................................................................

      mfm1_=mfm1
      do i=1,kz
         do ie=1,ke
            do ib=1,kb
               hibrz_(i,ie,ib)=hibrz(i,ie,ib)
            enddo
         enddo
      enddo

!..................................................................
!     Plot the FREYA birth points.
!..................................................................

        call frplteq(xpts,ypts,zpts,rpts,vx,vy,vz,ipts,curdep, &
          nfrplt,frplt)

      return
      end subroutine frnfreya
