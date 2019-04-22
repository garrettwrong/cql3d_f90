module zblock_mod


contains

      subroutine block(limits,abd,lda,lm,m,k)
      implicit integer (i-n), real*8 (a-h,o-z)
      dimension abd(lda,1), tmpsww1(65),tmpsww2(65)
      integer limits(3)

!-sww By including FORTRAN code for sectioning, two 64-element blocks
!-sww of local storage can be used to retain the needed multipliers
!-sww for processing an entire block (up to 64 rows)  using "SAXPY"
!-sww like operations.  Variables used for the sectioning process are :
!-sww IISAVE - the number of rows so far completed
!-sww LIMIT - the number of rows to process this pass
!-sww LMXX  - the number of rows left to processes after this pass

      mmsave=limits(3)
      iisave=0
      lmxx=lm

 683  limit=lmxx
      mm=mmsave
!-sww The section size, LIMIT, has an upper bound of 64 due to
!-sww the Cray hardware vector lengths.
      if(limit.gt.64) limit=64
      lmxx=lmxx-limit

!-sww Loop 687 copies the two (sectioned) multiplier columns (K,K+1)
!-sww into variables which are in local storage (see REGFILE statement)
!DIR$ SHORTLOOP
      do 687 ii=1,limit
        tmpsww1(ii) = abd(m+1+ii+iisave,k)
        tmpsww2(ii) = abd(m+ii+iisave,k+1)
 687  continue

!-sww The twelve "top-of-the-column" values are loaded outside of the
!-sww inner loop.  T1 to T6 are elements which are in the six columns
!-sww being processed AND they are in the Kth row.  T1b to T6b are
!-sww the elements immediately below T1 to T6.
!-sww In processing the residue columns, KP2 was updated to point
!-sww to the first column to be processed by the main loop.
      do 685 j=limits(1),limits(2),6
        t1=abd(mm-2,j)
        t2=abd(mm-3,j+1)
        t3=abd(mm-4,j+2)
        t4=abd(mm-5,j+3)
        t5=abd(mm-6,j+4)
        t6=abd(mm-7,j+5)
        t1b=abd(mm-1,j)
        t2b=abd(mm-2,j+1)
        t3b=abd(mm-3,j+2)
        t4b=abd(mm-4,j+3)
        t5b=abd(mm-5,j+4)
        t6b=abd(mm-6,j+5)
        mmiisave=mm+iisave
!-sww The directive "shortloop" lets the compiler know that the
!-sww number of iterations in loop 684 is 64 or less so that the
!-sww compiler will not generate its own sectioning code.
!DIR$ IVDEP
!DIR$ SHORTLOOP
        do 684 ii=1,limit
!-sww t0=tmpsww1(ii)
!-sww t0b=tmpsww2(ii)
          y1=abd(mmiisave+ii-1,j)
          y2=tmpsww1(ii)*t1
          y3=y1+y2
          y4=t1b*tmpsww2(ii)
          y5=abd(mmiisave+ii-2,j+1)
          y1=y3+y4
          y2=tmpsww1(ii)*t2
          y6=abd(mmiisave+ii-3,j+2)
          abd(mmiisave+ii-1,j) = y1
          y3=y2+y5
          y4=tmpsww2(ii)*t2b
          y1=y3+y4
          y2=tmpsww1(ii)*t3
          y56=abd(mmiisave+ii-6,j+5)
          abd(mmiisave+ii-2,j+1)=y1
          y3=y6+y2
          y4=tmpsww2(ii)*t3b

          y1=y3+y4
          y2=tmpsww1(ii)*t4
          y5=abd(mmiisave+ii-4,j+3)
          abd(mmiisave+ii-3,j+2) = y1
          y3=y2+y5
          y4=tmpsww2(ii)*t4b
          y1=y3+y4
          y2=tmpsww1(ii)*t5
          y6=abd(mmiisave+ii-5,j+4)
          abd(mmiisave+ii-4,j+3)=y1
          y3=y6+y2
          y4=tmpsww2(ii)*t5b

          y1=y3+y4
          y2=tmpsww1(ii)*t6
          abd(mmiisave+ii-5,j+4) = y1
          y3=y2+y56
          y4=tmpsww2(ii)*t6b
          y1=y3+y4
          abd(mmiisave+ii-6,j+5)=y1

 684    continue
!-sww Update banded matrix row index for having processed 6 columns
        mm=mm-6
 685  continue

!-sww IISAVE (number of rows processed) is updated
!-sww If more rows remain (lower sections); goto 683
      iisave=iisave+limit
      if(lmxx.ne.0) goto 683
      return
      end

end module zblock_mod
