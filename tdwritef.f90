module tdwritef_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

  subroutine tdwritef
    use cqlconf_mod, only : setup0
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
!.......................................................................
!     Save current parameters ,distribution function and spatial source
!     to enable complete restart of run
!.......................................................................

      include 'frname_decl.h'
      include 'frname.h'
!MPIINSERT_INCLUDE
!.......................................................................

!MPIINSERT_IF_RANK_NE_0_RETURN

      iunwrif=19
      open(unit=iunwrif,delim='apostrophe',file='distrfunc')

!.......................................................................
!l    1. Write current time-step and main diagnostics
!     inamlin=line number of first namelist (i.e. number of lines to
!     skip when reading this file). [innamlin NOT USED, BH070508].]
!.......................................................................

      write(iunwrif,9100) lrors+5,n,dtr,lrors,setup0%lrz,setup0%cqlpmod
      if (setup0%cqlpmod .ne. "enabled") then

        write(iunwrif,9101)
        do k=1,ngen
        do 101 l=1,lrors
!BH070408 write(iunwrif,9102) l,rovera(l),iy_(l),reden(kelecg,setup0%lrindx(l))
          write(iunwrif,9102) l,rovera(l),iy_(l),reden(k,setup0%lrindx(l)) &
            ,energy(k,setup0%lrindx(l)),totcurz(setup0%lrindx(l)),rovs(setup0%lrindx(l))
 101    continue
        enddo
!MPIINSERT_IF_RANK_EQ_0
        WRITE(*,*) &
        'tdwritef[nlwritf.ne."ncdfdist"]:Writing data into distrfunc.nc'
        WRITE(*,*)'tdwritef_43: For checkup SUM(reden),SUM(energy)=', &
         SUM(reden),SUM(energy)
        WRITE(*,*)'tdwritef_44: For checkup SUM(totcurz),SUM(rovs)=', &
         SUM(totcurz),SUM(rovs)
!MPIINSERT_ENDIF_RANK

      else

        write(iunwrif,9103)
        do k=1,ngen
        do 102 l=1,lrors
!BH070408 write(iunwrif,9102) l,sz(l),iy_(l),denpar(kelecg,setup0%lsindx(l)),
          write(iunwrif,9102) l,sz(l),iy_(l),denpar(k,setup0%lsindx(l)), &
            enrgypa(k,setup0%lsindx(l)),currmtpz(l),rovsloc(l)
 102    continue
        enddo

      endif

!.......................................................................
!l    2. Write namelists and distribution function
!.......................................................................

      write(iunwrif,'(" ")')
      ! nml name now private, can make a function in the module for this (iunwrif,setup0)
      write(iunwrif,setup)
      ! nml name now private, can make a function in the module for this write(iunwrif,trsetup)
      !nml name now private, can make a function in the module for this write(iunwrif,sousetup)
      !nml name now private, can make a function in the module for this write(iunwrif,eqsetup)
      !nml name now private, can make a function in the module for this write(iunwrif,rfsetup)
      write(iunwrif,frsetup)
!.......................................................................
!BH070408:  Have added write of frsetup here.  Probably should
!           also write the seed number for the random number generator,
!           so can achief same NBI results with restart.
!.......................................................................



!.......................................................................
!l    2.2 f(i,j,k,l)
!.......................................................................

!     Don't bother writing f, if nlwritf="ncdfdist", indicating
!     will use the netcdf setup0%mnemonic.nc file as source for f for
!     restart (nlrestrt="ncdfdist").

      if (setup0%nlwritf.ne."ncdfdist") then

!MPIINSERT_IF_RANK_EQ_0
         WRITE(*,*) &
           'tdwritef[nlwritf.ne."ncdfdist"]:Writing f into distrfunc.nc'
         WRITE(*,*)'tdwritef_83: For checkup SUM(f)=', SUM(f)
!MPIINSERT_ENDIF_RANK
!BH050328:  Problem with reading the given format with index
!BH050328:  numbers .lt.1.e-99
        do 220 k=1,ngen
        do 221 il=1,lrors
          do 222 j=1,jx
             do i=1,iy_(il)
                f(i,j,k,il)=max(f(i,j,k,il),1.d-99)
             enddo
             write(iunwrif,9220) (f(i,j,k,il),i=1,iy_(il))
 222      continue
 221    continue
 220    continue
!MPIINSERT_IF_RANK_EQ_0
         WRITE(*,*)'tdwritef_98: For checkup SUM(f)=', SUM(f)
!MPIINSERT_ENDIF_RANK

      endif ! nlwritf.ne."ncdfdist"

      if (setup0%cqlpmod.eq."enabled") then
!.......................................................................
!l    2.3 Spatial source term: spasou(i,j,k,l)
!.......................................................................
         do 230 k=1,ngen
            do 231 il=1,lrors
               do 232 j=1,jx
                  write(iunwrif,9220) (spasou(i,j,k,il),i=1,iy_(il))
 232           continue
 231        continue
 230     continue
!.......................................................................
!l    2.4 Velocity source term: velsou(i,j,k,l)
!.......................................................................
         do 240 k=1,ngen
            do 241 il=1,lrors
               do 242 j=1,jx
                  write(iunwrif,9220) (velsou(i,j,k,il),i=1,iy_(il))
 242           continue
 241        continue
 240     continue
      endif ! setup0%cqlpmod=enabled

!.......................................................................
 9100 format(i3," more lines to skip before reading namelist and f",//, &
        "time-step n=",i4," dtr=",1pe12.5," lrors=",i3," setup0%lrz=",i3, &
        " setup0%cqlpmod=",a)
 9101 format(/,"  l",3x,"rovera",4x," iy ",3x,"density",5x,"energy",4x, &
        "tot.curr.",2x,"res/spitzr")
 9102 format(i3,1pe12.4,i4,4e12.4)
 9103 format(/,"  l",6x,"s",6x," iy ",3x,"density",5x,"energy",4x, &
        "tot.curr.",2x,"res/spitzr")
!BH080201 9220 format(1p10e13.6)
 9220 format(1p10e14.6)
!.......................................................................
      close(unit=iunwrif)
      !pause
      return
      end subroutine tdwritef

end module tdwritef_mod
