c
c
      subroutine aindflt1
      implicit integer (i-n), real*8 (a-h,o-z)

c..................................................................
c     Set defaults for some basic variables depending on the
c     namelist values (all sets except frsetup) or parameter values.
c..................................................................

      include 'param.h'
      include 'comm.h'

      do ll=0,lrorsa
        indxlr(ll)=ll
        indxls(ll)=ll
cBH080305        mplot(ll)="disabled"
      enddo
      do ll=1,lrorsa
        mplot(ll)="disabled"
      enddo
      
      mbet=mbeta
      ntotal=ntotala
      mxp1=mx+1

c     analegco="enabled" => ngauss should be .le. 0
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
c.......................................................................
c     lrza arrays
c.......................................................................
      do ll=1,lrza
        lorbit(ll)=lfielda
        currxj0(ll)=0.0
      enddo
      currxj0(0)=0.0

c.......................................................................
c     lrorsa arrays
c.......................................................................
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

