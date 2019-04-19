module achief1_mod

!
!

contains

subroutine achief1
use param_mod
use comm_mod
use pltmain_mod, only : pltmain
use r8subs_mod, only : dcopy
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
!all aindflt
!all eqindflt
!all aindflt1

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

!all ainspec

!.......................................................................
!     set variables dependent on input variables
!.......................................................................

!all ainsetva

!..................................................................
!     Allocate arrays , if required
!..................................................................

!all ainalloc

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
!all ain_transcribe("cqlinput")
else
write(6,*)
write(6,*) 'mnemonic = ',mnemonic
write(6,*)
endif

!.....................................................................
!     Determine mesh normalization constant vnorm.
!.....................................................................

!all ainvnorm

!.....................................................................
!     Call the initialization routines for the appended modules..
!.....................................................................

!all eqinitl
!all frinitl

open(unit=2,file="cqlinput",delim='apostrophe',status="old")
!all frset(lrz,noplots,nmlstout)   ! Uses unit 2
!lose(2)

!..................................................................
!     Call an initialization routine which determines flux surface
!     geometry and magnetic field structure.
!..................................................................

!all aingeom

!.......................................................................
!     Initialize mesh along magnetic field line
!.......................................................................

if (cqlpmod.eq."enabled" .and. numclas.eq.1 .and. ls.eq.lsmax)then
lz=lz/2+1
lsmax=lsmax/2+1
ls=ls/2+1
endif

!all micxiniz

if (cqlpmod.eq."enabled" .and. numclas.eq.1 .and. ls.eq.lsmax)then
lz=2*(lz-1)
lsmax=2*(lsmax-1)
ls=2*(ls-1)
!all wploweq
endif

!.......................................................................
!     Initialize some plasma parameters
!.......................................................................

!all ainpla

!.......................................................................
!     call a routine to determine meshes y, x and related quantities
!.......................................................................

!all micxinit

ieq_tot=inewjx_(1) ! inewjx_() is defined in micxinit
ieq_(1)=1 ! Eqn no. at beginning of each flux surface
ieq_(lrors+1)=ieq_tot ! lrors+1 should be 2 here

!............................................................
!     call main initialization routine.
!............................................................

!all ainitial

if (nstop.eq.0) then
!all pltmain
write(*,*) 'In ACHIEF1, before call pgend'
!all pgend
stop 'achief1: nstop=0'
endif

!..................................................................
!     Initialize main netCDF write, if netcdfnm.ne."disabled"
!..................................................................

if (netcdfnm.ne."disabled") then
!all netcdfrw2(0)
endif

!.......................................................................
!     Solve equations on the flux surface
!.......................................................................

!all tdnflxs(1)
!     Copy current distribution f into f_
!all dcopy(iyjx2*ngen*lrors,f(0:iyjx2*ngen*lrors-1,0,1,1),1,
!+     f_(0:iyjx2*ngen*lrors-1,0,1,1),1)
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
!all achiefn(0) ! get solution for new f. !
!-------------------------------------------!
! Start time advancement:
if(nefiter.eq.1) then
n=n+1
n_(1)=n ! new time-step for this flux surface
! for 2-d (v_par,v_perp) calculation ntloop controls
! end of run or restart.
! Also updates time.
!all ntloop
endif

!all tdnflxs(1)
!all cfpgamma ! Re-calc. Coul.Log for the new distr.func.
do k=1,ngen  ! Compute density gains and losses, and powers.
! For lbdry0='disabled',  Redefine f at v=0 so it is unique:
! (For lbdry0='enabled', coeff matrix is set up
!   to automatically maintain unicity.)
if (lbdry0.ne."enabled") then !-YuP: moved here from impavnc0
!all dcopy(iyjx2,f(0:iyjx2-1,0,k,l_),1,
!+             fxsp(0:iyjx2-1,0,k,l_),1)
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
!all diagscal(k) !-> renorm f() if lbdry(k)="scale"
!all coefstup(k) ! To define da...df coeffs, gon(i,j), etc
!all coefmidv(da,1)
!all coefmidv(db,2)
!all coefmidv(dc,3)
!all coefmidt(dd,1)
!all coefmidt(de,2)
!all coefmidt(df,3)
!all coefwtj(k)
!all coefwti(k)
!all diagimpd(k)
enddo ! k
!all achiefn(1)  !Compute plasma energy, density and energy transfer


return
end
end module achief1_mod
