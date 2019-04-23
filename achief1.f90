module achief1_mod

  !---BEGIN USE

  use achiefn_mod, only : achiefn
  use ainalloc_mod, only : ainalloc
  use aindflt1_mod, only : aindflt1
  use aindflt_mod, only : aindflt
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
  use eqindflt_mod, only : eqindflt
  use eqinitl_mod, only : eqinitl
  use micxinit_mod, only : micxinit
  use micxiniz_mod, only : micxiniz
  use ntloop_mod, only : ntloop
  use pltmain_mod, only : pltmain
  use profiles_mod, only : profiles
  use r8subs_mod, only : dcopy
  use tdnflxs_mod, only : tdnflxs
  use wploweq_mod, only : wploweq

  !---END USE

!
!

contains

subroutine achief1
use param_mod
use comm_mod
use pltmain_mod, only : pltmain
use r8subs_mod, only : dcopy
use ainvnorm_mod, only : ainvnorm
use ainspec_mod, only : ainspec
use ainsetva_mod, only : ainsetva
use ainpla_mod, only : ainpla
use aingeom_mod, only : aingeom
use ainitial_mod, only : ainitial
use aindflt_mod, only : aindflt
use aindflt1_mod, only : aindflt1
use ainalloc_mod, only : ainalloc
use aindfpa_mod , only : ain_transcribe
use abchief_mod, only : abchief
use achiefn_mod, only : achiefn
 
implicit integer (i-n), real*8 (a-h,o-z)
save

!..................................................................
!     This routine directs the calculation for lrzmax=1
!..................................................................


include 'name.h90'
!.......................................................................

!..................................................................
!     Set defaults - for main code + "eq" module.
!..................................................................
call aindflt
call eqindflt
call aindflt1

!.....................................................................
!     Read in driver input namelist setup
!.....................................................................
open(unit=2,file="cqlinput",status="old")
read(2,setup)
read(2,trsetup)
read(2,sousetup)
read(2,eqsetup)
read(2,rfsetup)
!lose(2)

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

if (nmlstout.eq."enabled") then
write(6,*)'  In achief1: '
write(6,setup0)
write(6,setup)
write(6,trsetup)
write(6,sousetup)
write(6,eqsetup)
write(6,rfsetup)
elseif (nmlstout.eq."trnscrib") then
write(6,*)'  In achief1: '
call ain_transcribe("cqlinput")
else
write(6,*)
write(6,*) 'mnemonic = ',mnemonic
write(6,*)
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

open(unit=2,file="cqlinput",delim='apostrophe',status="old")
call frset(lrz,noplots,nmlstout)   ! Uses unit 2
!lose(2)

!..................................................................
!     Call an initialization routine which determines flux surface
!     geometry and magnetic field structure.
!..................................................................

call aingeom

!.......................................................................
!     Initialize mesh along magnetic field line
!.......................................................................

if (cqlpmod.eq."enabled" .and. numclas.eq.1 .and. ls.eq.lsmax)then
lz=lz/2+1
lsmax=lsmax/2+1
ls=ls/2+1
endif

call micxiniz

if (cqlpmod.eq."enabled" .and. numclas.eq.1 .and. ls.eq.lsmax)then
lz=2*(lz-1)
lsmax=2*(lsmax-1)
ls=2*(ls-1)
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
write(*,*) 'In ACHIEF1, before call pgend'
call pgend
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
     f(0:iyjx2*ngen*lrors-1,0,1,1),1, &
     f_(0:iyjx2*ngen*lrors-1,0,1,1),1)
!     bring background profiles up to time step n
if(nefiter.eq.1) call profiles
! Reset time step if (n+1).eq.nondtr1(i). .AND. LRZMAX=1
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
call dcopy(iyjx2,f(0:iyjx2-1,0,k,l_),1, &
     fxsp(0:iyjx2-1,0,k,l_),1)
s=0.
t=0.
do 2100 i=1,iy
s=s+vptb(i,lr_)*cynt2(i,l_)
t=t+vptb(i,lr_)*cynt2(i,l_)*f(i,1,k,l_)
2100        continue
do 2200 i=1,iy
f(i,1,k,l_)=t/s
2200        continue
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
end
end module achief1_mod
