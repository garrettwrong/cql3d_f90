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

module achief1_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use achiefn_mod, only : achiefn
  use ainalloc_mod, only : ainalloc
  use aindflt1_mod, only : aindflt1
  use aindfpa_mod, only : ain_transcribe
  use aingeom_mod, only : aingeom
  use ainitial_mod, only : ainitial
  use ainpla_mod, only : ainpla
  use ainsetva_mod, only : ainsetva
  use ainspec_mod, only : ainspec
  use ainvnorm_mod, only : ainvnorm
  use cfpgamma_mod, only : cfpgamma
  use coefmidt_mod, only : coefmidt
  use coefmidv_mod, only : coefmidv
  use coefstup_mod, only : coefstup
  use coefwti_mod, only : coefwti
  use coefwtj_mod, only : coefwtj
  use diagimpd_mod, only : diagimpd
  use diagscal_mod, only : diagscal
  use eqinitl_mod, only : eqinitl
  use micxinit_mod, only : micxinit
  use micxiniz_mod, only : micxiniz
  use ntloop_mod, only : ntloop
  use pltmain_mod, only : pltmain
  use profiles_mod, only : profiles
  use r8subs_mod, only : dcopy
  use tdnflxs_mod, only : tdnflxs
  use wploweq_mod, only : wploweq
  use netcdfrw2_mod, only : netcdfrw2

  !---END USE

  !
  !

contains

  subroutine achief1(nml_file)
    use param_mod
    use cqlconf_mod, only : setup0
    use cqlconf_mod, only : get_setup_from_nml
    use cqlconf_mod, only : get_eqsetup_from_nml
    use cqlconf_mod, only : get_rfsetup_from_nml
    use cqlconf_mod, only : get_trsetup_from_nml
    use cqlconf_mod, only : get_sousetup_from_nml
    use cqlconf_mod, only : print_all_conf_nml
    use cqlcomm_mod
    use pltmain_mod, only : pltmain
    use r8subs_mod, only : dcopy
    use ainvnorm_mod, only : ainvnorm
    use ainspec_mod, only : ainspec
    use ainsetva_mod, only : ainsetva
    use ainpla_mod, only : ainpla
    use aingeom_mod, only : aingeom
    use ainitial_mod, only : ainitial
    use aindflt1_mod, only : aindflt1
    use ainalloc_mod, only : ainalloc
    use aindfpa_mod , only : ain_transcribe
    use achiefn_mod, only : achiefn

    implicit integer (i-n), real(c_double) (a-h,o-z)
    character(len=*), intent(in), optional :: nml_file
    save

    !..................................................................
    !     This routine directs the calculation for setup0%lrzmax=1
    !..................................................................

    !..................................................................
    !     Set defaults not in the cql_conf types
    !..................................................................

    call aindflt1

    !.....................................................................
    !     Read in driver input namelist setup
    !.....................................................................
    if(present(nml_file)) then
       call get_setup_from_nml(nml_file, close_nml_file=.FALSE., debug_print=.FALSE.)
       call get_trsetup_from_nml(nml_file, close_nml_file=.FALSE., debug_print=.FALSE.)
       call get_sousetup_from_nml(nml_file, close_nml_file=.FALSE., debug_print=.FALSE.)
       call get_eqsetup_from_nml(nml_file, close_nml_file=.FALSE., debug_print=.FALSE.)
       call get_rfsetup_from_nml(nml_file, close_nml_file=.TRUE., debug_print=.FALSE.)
    end if
    !..................................................................
    !     Call routine which finds electron and ion species indices.
    !..................................................................

    call ainspec

    !.......................................................................
    !     set variables dependent on input variables
    !.......................................................................

    call ainsetva

    !..................................................................
    !     Allocate arrays , if required
    !..................................................................

    call ainalloc

    !.......................................................................
    !     print namelists
    !.......................................................................

    if (setup0%nmlstout.eq."enabled") then
       if(setup0%verbose>0) write(6,*)'  In achief1: '
       call print_all_conf_nml
    elseif (setup0%nmlstout.eq."trnscrib") then
       if(setup0%verbose>0) write(6,*)'  In achief1: '
       call ain_transcribe(nml_file)
    else
       if(setup0%verbose>0) write(6,*)
       if(setup0%verbose>0) write(6,*) 'setup0%mnemonic = ',setup0%mnemonic
       if(setup0%verbose>0) write(6,*)
    endif

    !.....................................................................
    !     Determine mesh normalization constant vnorm.
    !.....................................................................

    call ainvnorm

    !.....................................................................
    !     Call the initialization routines for the appended modules..
    !.....................................................................

    call eqinitl
    call frinitl

    !NMLXXX, frsetup nml not converted yet
    if(present(nml_file)) then
       open(unit=2,file=nml_file,delim='apostrophe',status="old")
       call frset(setup0%lrz,setup0%noplots,setup0%nmlstout)   ! Uses unit 2
       close(2)
    end if

    !..................................................................
    !     Call an initialization routine which determines flux surface
    !     geometry and magnetic field structure.
    !..................................................................

    call aingeom

    !.......................................................................
    !     Initialize mesh along magnetic field line
    !.......................................................................

    if (setup0%cqlpmod.eq."enabled" .and. numclas.eq.1 .and. setup0%ls.eq.setup0%lsmax)then
       lz=lz/2+1
       setup0%lsmax=setup0%lsmax/2+1
       setup0%ls=setup0%ls/2+1
    endif

    call micxiniz

    if (setup0%cqlpmod.eq."enabled" .and. numclas.eq.1 .and. setup0%ls.eq.setup0%lsmax)then
       lz=2*(lz-1)
       setup0%lsmax=2*(setup0%lsmax-1)
       setup0%ls=2*(setup0%ls-1)
       call wploweq
    endif

    !.......................................................................
    !     Initialize some plasma parameters
    !.......................................................................

    call ainpla

    !.......................................................................
    !     call a routine to determine meshes y, x and related quantities
    !.......................................................................

    call micxinit

    ieq_tot=inewjx_(1) ! inewjx_() is defined in micxinit
    ieq_(1)=1 ! Eqn no. at beginning of each flux surface
    ieq_(lrors+1)=ieq_tot ! lrors+1 should be 2 here

    !............................................................
    !     call main initialization routine.
    !............................................................

    call ainitial

    if (nstop.eq.0) then
       call pltmain
#ifndef NOPGPLOT
       if(setup0%verbose>0) write(*,*) 'In ACHIEF1, before call pgend'
#endif
#ifndef NOPGPLOT
       call pgend
#endif
       stop 'achief1: nstop=0'
    endif

    !..................................................................
    !     Initialize main netCDF write, if netcdfnm.ne."disabled"
    !..................................................................

    if (netcdfnm.ne."disabled") then
       call netcdfrw2(0)
    endif

    !.......................................................................
    !     Solve equations on the flux surface
    !.......................................................................

    call tdnflxs(1)
    !     Copy current distribution f into f_
    call dcopy(iyjx2*ngen*lrors, &
          f(0:iy+1,0:jx+1,1:ngen,1:lrors),1, &
         f_(0:iy+1,0:jx+1,1:ngen,1:lrors),1)
    !     bring background profiles up to time step n
    if(nefiter.eq.1) call profiles
    ! Reset time step if (n+1).eq.nondtr1(i). .AND. setup0%lrzmax=1
    do i=1,ndtr1a
       if ((n+1).eq.nondtr1(i)) then
          dtr=dtr1(i)
          dtreff=dtr
          dttr=dtr*nrstrt
       endif
    enddo
    !-------------------------------------------!
    call achiefn(0) ! get solution for new f. !
    !-------------------------------------------!
    ! Start time advancement:
    if(nefiter.eq.1) then
       n=n+1
       n_(1)=n ! new time-step for this flux surface
       ! for 2-d (v_par,v_perp) calculation ntloop controls
       ! end of run or restart.
       ! Also updates time.
       call ntloop
    endif

    call tdnflxs(1)
    call cfpgamma ! Re-calc. Coul.Log for the new distr.func.
    do k=1,ngen  ! Compute density gains and losses, and powers.
       ! For lbdry0='disabled',  Redefine f at v=0 so it is unique:
       ! (For lbdry0='enabled', coeff matrix is set up
       !   to automatically maintain unicity.)
       if (lbdry0.ne."enabled") then !-YuP: moved here from impavnc0
          call dcopy(iyjx2,f(0:iy+1,0:jx+1,k,l_),1, &
                        fxsp(0:iy+1,0:jx+1,k,l_),1)
          s=0.
          t=0.
          do 2100 i=1,iy
             s=s+vptb(i,lr_)*cynt2(i,l_)
             t=t+vptb(i,lr_)*cynt2(i,l_)*f(i,1,k,l_)
2100      end do
          do 2200 i=1,iy
             f(i,1,k,l_)=t/s
2200      end do
       endif
       call diagscal(k) !-> renorm f() if lbdry(k)="scale"
       call coefstup(k) ! To define da...df coeffs, gon(i,j), etc
       call coefmidv(da,1)
       call coefmidv(db,2)
       call coefmidv(dc,3)
       call coefmidt(dd,1)
       call coefmidt(de,2)
       call coefmidt(df,3)
       call coefwtj(k)
       call coefwti(k)
       call diagimpd(k)
    enddo ! k
    call achiefn(1)  !Compute plasma energy, density and energy transfer


    return
  end subroutine achief1

end module achief1_mod
