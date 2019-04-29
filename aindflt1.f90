module aindflt1_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine aindflt1
      use param_mod
      use comm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     Set defaults for some basic variables depending on the
!     namelist values (all sets except frsetup) or parameter values.
!..................................................................


      do ll=0,lrorsa
        indxlr(ll)=ll
        indxls(ll)=ll
!BH080305        mplot(ll)="disabled"
      enddo
      do ll=1,lrorsa
        mplot(ll)="disabled"
      enddo

      mbet=mbeta
      ntotal=ntotala
      mxp1=mx+1

!     analegco="enabled" => ngauss should be .le. 0
      analegco="enabled"
      elpar0=0.
      lmdpln_=lmidpln
      ipxy=min(51,iy)
      jpxy=min(101,jx+1)
      if (mod(jpxy,2).eq.0) jpxy=jpxy-1
      irstart=0
      l_=1
      lr_=l_
      ls_=1
      indxlr_=1
      indxls_=1
      impadi=0
      xmax=1.
!.......................................................................
!     lrza arrays
!.......................................................................
      do ll=1,lrza
        lorbit(ll)=lfielda
        currxj0(ll)=0.0
      enddo
      currxj0(0)=0.0

!.......................................................................
!     lrorsa arrays
!.......................................................................
      do 105 ll=1,lrorsa
        itl_(ll)=iy/2
        itu_(ll)=iy-itl_(ll)+1
        n_(ll)=0
        nch(ll)=1
        iy_(ll)=iy
        iyh_(ll)=iy/2
        iyjx_(ll)=iy*jx
        time_(ll)=0.d0
 105  continue
      timet=0.d0
      n=0

      cursign=+1.0

      frmodp="disabled"
      beamplsep="disabled"
      fr_gyrop="disabled"
      beamponp=zero
      beampoffp=zero

      nnz=nnza
      nnr=nnra

      return
      end

end module aindflt1_mod
