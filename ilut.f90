module ilut_mod

  !---BEGIN USE

  use iso_c_binding, only : c_double
  use r8subs_mod, only : daxpy

  !---END USE

!----------------------------------------------------------------------c
!                          S P A R S K I T                             c
!----------------------------------------------------------------------c
!                   ITERATIVE SOLVERS MODULE                           c
!----------------------------------------------------------------------c
! This Version Dated: August 13, 1996. Warning: meaning of some        c
! ============ arguments have changed w.r.t. earlier versions. Some    c
!              Calling sequences may also have changed                 c
!----------------------------------------------------------------------c
! Contents:                                                            c
!-------------------------preconditioners------------------------------c
!                                                                      c
! ILUT    : Incomplete LU factorization with dual truncation strategy  c
! ILUTP   : ILUT with column  pivoting                                 c
! ILUD    : ILU with single dropping + diagonal compensation (~MILUT)  c
! ILUDP   : ILUD with column pivoting                                  c
! ILUK    : level-k ILU                                                c
! ILU0    : simple ILU(0) preconditioning                              c
! MILU0   : MILU(0) preconditioning                                    c
!                                                                      c
!----------sample-accelerator-and-LU-solvers---------------------------c
!                                                                      c
! PGMRES  : preconditioned GMRES solver                                c
! LUSOL   : forward followed by backward triangular solve (Precond.)   c
! LUTSOL  : solving v = (LU)^{-T} u (used for preconditioning)         c
!                                                                      c
!-------------------------utility-routine------------------------------c
!                                                                      c
! QSPLIT  : quick split routine used by ilut to sort out the k largest c
!           elements in absolute value                                 c
!                                                                      c
!----------------------------------------------------------------------c
!                                                                      c
! Note: all preconditioners are preprocessors to pgmres.               c
! usage: call preconditioner then call pgmres                          c
!                                                                      c
!----------------------------------------------------------------------c

contains

      subroutine ilut(n,a,ja,ia,lfil,droptol,alu,jlu,ju,iwk,w,jw,ierr)
!-----------------------------------------------------------------------
      implicit none
      integer n
!BH070419
      real(c_double) a(*),alu(*),w(n+1),droptol
!BH070419
      integer ja(*),ia(n+1),jlu(*),ju(n),jw(2*n),lfil,iwk,ierr
!yup      real(c_double) a(*),alu(*),w(*),droptol
!yup      integer ja(*),ia(*),jlu(*),ju(*),jw(*),lfil,iwk,ierr
!----------------------------------------------------------------------*
!                      *** ILUT preconditioner ***                     *
!      incomplete LU factorization with dual truncation mechanism      *
!----------------------------------------------------------------------*
!     Author: Yousef Saad *May, 5, 1990, Latest revision, August 1996  *
!----------------------------------------------------------------------*
! PARAMETERS
!-----------
!
! on entry:
!==========
! n       = integer. The row dimension of the matrix A. The matrix
!
! a,ja,ia = matrix stored in Compressed Sparse Row format.
!
! lfil    = integer. The fill-in parameter. Each row of L and each row
!           of U will have a maximum of lfil elements (excluding the
!           diagonal element). lfil must be .ge. 0.
!           ** WARNING: THE MEANING OF LFIL HAS CHANGED WITH RESPECT TO
!           EARLIER VERSIONS.
!
! droptol = real(c_double). Sets the threshold for dropping small terms in the
!           factorization. See below for details on dropping strategy.
!
!
! iwk     = integer. The lengths of arrays alu and jlu. If the arrays
!           are not big enough to store the ILU factorizations, ilut
!           will stop with an error message.
!
! On return:
!===========
!
! alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
!           the L and U factors together. The diagonal (stored in
!           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
!           contains the i-th row of L (excluding the diagonal entry=1)
!           followed by the i-th row of U.
!
! ju      = integer array of length n containing the pointers to
!           the beginning of each row of U in the matrix alu,jlu.
!
! ierr    = integer. Error message with the following meaning.
!           ierr  = 0    --> successful return.
!           ierr .gt. 0  --> zero pivot encountered at step number ierr.
!           ierr  = -1   --> Error. input matrix may be wrong.
!                            (The elimination process has generated a
!                            row in L or U whose length is .gt.  n.)
!           ierr  = -2   --> The matrix L overflows the array al.
!           ierr  = -3   --> The matrix U overflows the array alu.
!           ierr  = -4   --> Illegal value for lfil.
!           ierr  = -5   --> zero row encountered.
!
! work arrays:
!=============
! jw      = integer work array of length 2*n.
! w       = real work array of length n+1.
!
!----------------------------------------------------------------------
! w, ju (1:n) store the working array [1:ii-1 = L-part, ii:n = u]
! jw(n+1:2n)  stores nonzero indicators
!
! Notes:
! ------
! The diagonal elements of the input matrix must be  nonzero (at least
! 'structurally').
!
!----------------------------------------------------------------------*
!---- Dual drop strategy works as follows.                             *
!                                                                      *
!     1) Theresholding in L and U as set by droptol. Any element whose *
!        magnitude is less than some tolerance (relative to the abs    *
!        value of diagonal element in u) is dropped.                   *
!                                                                      *
!     2) Keeping only the largest lfil elements in the i-th row of L   *
!        and the largest lfil elements in the i-th row of U (excluding *
!        diagonal elements).                                           *
!                                                                      *
! Flexibility: one  can use  droptol=0  to get  a strategy  based on   *
! keeping  the largest  elements in  each row  of L  and U.   Taking   *
! droptol .ne.  0 but lfil=n will give  the usual threshold strategy   *
! (however, fill-in is then mpredictible).                             *
!----------------------------------------------------------------------*
!     locals
      integer ju0,k,j1,j2,j,ii,i,lenl,lenu,jj,jrow,jpos,len
      real(c_double) tnorm, t, abs, s, fact
      if (lfil .lt. 0) goto 998
!-----------------------------------------------------------------------
!     initialize ju0 (points to next element to be added to alu,jlu)
!     and pointer array.
!-----------------------------------------------------------------------
      ju0 = n+2
      jlu(1) = ju0
!
!     initialize nonzero indicator array.
!
      do 1 j=1,n
         jw(n+j)  = 0
 1    continue
!-----------------------------------------------------------------------
!     beginning of main loop.
!-----------------------------------------------------------------------
      do 500 ii = 1, n
         j1 = ia(ii)
         j2 = ia(ii+1) - 1
         tnorm = 0.0d0
         do 501 k=j1,j2
            tnorm = tnorm+abs(a(k))
 501     continue
         if (tnorm .eq. 0.0) then
            write(*,*)'ilut (in SPARSKIT):ii,ia(ii),ia(ii+1),a(j1:j2)=', &
                                          ii,ia(ii),ia(ii+1),a(j1:j2)
            goto 999
         endif
         tnorm = tnorm/DBLE(j2-j1+1)
!
!     unpack L-part and U-part of row of A in arrays w
!
         lenu = 1
         lenl = 0
         jw(ii) = ii
         w(ii) = 0.0
         jw(n+ii) = ii
!
         do 170  j = j1, j2
            k = ja(j)
            t = a(j)
            if (k .lt. ii) then
               lenl = lenl+1
               jw(lenl) = k
               w(lenl) = t
               jw(n+k) = lenl
            else if (k .eq. ii) then
               w(ii) = t
            else
               lenu = lenu+1
               jpos = ii+lenu-1
               jw(jpos) = k
               w(jpos) = t
               jw(n+k) = jpos
            endif
 170     continue
         jj = 0
         len = 0
!
!     eliminate previous rows
!
 150     jj = jj+1
         if (jj .gt. lenl) goto 160
!-----------------------------------------------------------------------
!     in order to do the elimination in the correct order we must select
!     the smallest column index among jw(k), k=jj+1, ..., lenl.
!-----------------------------------------------------------------------
         jrow = jw(jj)
         k = jj
!
!     determine smallest column index
!
         do 151 j=jj+1,lenl
            if (jw(j) .lt. jrow) then
               jrow = jw(j)
               k = j
            endif
 151     continue
!
         if (k .ne. jj) then
!     exchange in jw
            j = jw(jj)
            jw(jj) = jw(k)
            jw(k) = j
!     exchange in jr
            jw(n+jrow) = jj
            jw(n+j) = k
!     exchange in w
            s = w(jj)
            w(jj) = w(k)
            w(k) = s
         endif
!
!     zero out element in row by setting jw(n+jrow) to zero.
!
         jw(n+jrow) = 0
!
!     get the multiplier for row to be eliminated (jrow).
!
         fact = w(jj)*alu(jrow)
         if (abs(fact) .le. droptol) goto 150
!
!     combine current row and row jrow
!
         do 203 k = ju(jrow), jlu(jrow+1)-1
            s = fact*alu(k)
            j = jlu(k)
            jpos = jw(n+j)
            if (j .ge. ii) then
!
!     dealing with upper part.
!
               if (jpos .eq. 0) then
!
!     this is a fill-in element
!
                  lenu = lenu+1
                  if (lenu .gt. n) goto 995
                  i = ii+lenu-1
                  jw(i) = j
                  jw(n+j) = i
                  w(i) = - s
               else
!
!     this is not a fill-in element
!
                  w(jpos) = w(jpos) - s

               endif
            else
!
!     dealing  with lower part.
!
               if (jpos .eq. 0) then
!
!     this is a fill-in element
!
                  lenl = lenl+1
                  if (lenl .gt. n) goto 995
                  jw(lenl) = j
                  jw(n+j) = lenl
                  w(lenl) = - s
               else
!
!     this is not a fill-in element
!
                  w(jpos) = w(jpos) - s
               endif
            endif
 203     continue
!
!     store this pivot element -- (from left to right -- no danger of
!     overlap with the working elements in L (pivots).
!
         len = len+1
         w(len) = fact
         jw(len)  = jrow
         goto 150
 160     continue
!
!     reset double-pointer to zero (U-part)
!
         do 308 k=1, lenu
            jw(n+jw(ii+k-1)) = 0
 308     continue
!
!     update L-matrix
!
         lenl = len
         len = min0(lenl,lfil)
!
!     sort by quick-split
!
         call qsplit (w,jw,lenl,len)
!
!     store L-part
!
         do 204 k=1, len
            if (ju0 .gt. iwk) goto 996
            alu(ju0) =  w(k)
            jlu(ju0) =  jw(k)
            ju0 = ju0+1
 204     continue
!
!     save pointer to beginning of row ii of U
!
         ju(ii) = ju0
!
!     update U-matrix -- first apply dropping strategy
!
         len = 0
         do k=1, lenu-1
            if (abs(w(ii+k)) .gt. droptol*tnorm) then
               len = len+1
               w(ii+len) = w(ii+k)
               jw(ii+len) = jw(ii+k)
            endif
         enddo
         lenu = len+1
         len = min0(lenu,lfil)
!
         call qsplit (w(ii+1), jw(ii+1), lenu-1,len)
!
!     copy
!
         t = abs(w(ii))
         if (len + ju0 .gt. iwk) goto 997
         do 302 k=ii+1,ii+len-1
            jlu(ju0) = jw(k)
            alu(ju0) = w(k)
            t = t + abs(w(k) )
            ju0 = ju0+1
 302     continue
!
!     store inverse of diagonal element of u
!
         if (w(ii) .eq. 0.0) w(ii) = (0.0001 + droptol)*tnorm
!
         alu(ii) = 1.0d0/ w(ii)
!
!     update pointer to beginning of next row of U.
!
         jlu(ii+1) = ju0
!-----------------------------------------------------------------------
!     end main loop
!-----------------------------------------------------------------------
 500  continue
      ierr = 0
      return
!
!     incomprehensible error. Matrix must be wrong.
!
 995  ierr = -1
      return
!
!     insufficient storage in L.
!
 996  ierr = -2
      return
!
!     insufficient storage in U.
!
 997  ierr = -3
      return
!
!     illegal lfil entered.
!
 998  ierr = -4
      return
!
!     zero row encountered
!
 999  ierr = -5
      return
!----------------end-of-ilut--------------------------------------------
!-----------------------------------------------------------------------
      end
!----------------------------------------------------------------------

!-----------------------------------------------------------------------
       subroutine pgmres(n, im, rhs, sol, vv, eps, maxits, iout, &
                          aa, ja, ia, alu, jlu, ju, ierr)
!-----------------------------------------------------------------------
       use r8subs_mod, only : ddot, daxpy
       implicit real(c_double) (a-h,o-z)
!BH070419
       integer n, im, maxits, iout, ierr, ja(*), ia(n+1), jlu(*), ju(n)
!BH070419
       real(c_double) vv(n,*), rhs(n), sol(n), aa(*), alu(*), eps
!yup       integer n, im, maxits, iout, ierr, ja(*), ia(*), jlu(*), ju(*)
!yup       real(c_double) vv(n,*), rhs(*), sol(*), aa(*), alu(*), eps
!----------------------------------------------------------------------*
!                                                                      *
!                 *** ILUT - Preconditioned GMRES ***                  *
!                                                                      *
!----------------------------------------------------------------------*
! This is a simple version of the ILUT preconditioned GMRES algorithm. *
! The ILUT preconditioner uses a dual strategy for dropping elements   *
! instead  of the usual level of-fill-in approach. See details in ILUT *
! subroutine documentation. PGMRES uses the L and U matrices generated *
! from the subroutine ILUT to precondition the GMRES algorithm.        *
! The preconditioning is applied to the right. The stopping criterion  *
! utilized is based simply on reducing the residual norm by epsilon.   *
! This preconditioning is more reliable than ilu0 but requires more    *
! storage. It seems to be much less prone to difficulties related to   *
! strong nonsymmetries in the matrix. We recommend using a nonzero tol *
! (tol=.005 or .001 usually give good results) in ILUT. Use a large    *
! lfil whenever possible (e.g. lfil = 5 to 10). The higher lfil the    *
! more reliable the code is. Efficiency may also be much improved.     *
! Note that lfil=n and tol=0.0 in ILUT  will yield the same factors as *
! Gaussian elimination without pivoting.                               *
!                                                                      *
! ILU(0) and MILU(0) are also provided for comparison purposes         *
! USAGE: first call ILUT or ILU0 or MILU0 to set up preconditioner and *
! then call pgmres.                                                    *
!----------------------------------------------------------------------*
! Coded by Y. Saad - This version dated May, 7, 1990.                  *
!----------------------------------------------------------------------*
! parameters                                                           *
!-----------                                                           *
! on entry:                                                            *
!==========                                                            *
!                                                                      *
! n     == integer. The dimension of the matrix.                       *
! im    == size of krylov subspace:  should not exceed 50 in this      *
!          version (can be reset by changing parameter command for     *
!          kmax below)                                                 *
! rhs   == real vector of length n containing the right hand side.     *
!          Destroyed on return.                                        *
! sol   == real vector of length n containing an initial guess to the  *
!          solution on input. approximate solution on output           *
! eps   == tolerance for stopping criterion. process is stopped        *
!          as soon as ( ||.|| is the euclidean norm):                  *
!          || current residual||/||initial residual|| <= eps           *
! maxits== maximum number of iterations allowed                        *
! iout  == output unit number number for printing intermediate results *
!          if (iout .le. 0) nothing is printed out.                    *
!                                                                      *
! aa, ja,                                                              *
! ia    == the input matrix in compressed sparse row format:           *
!          aa(1:nnz)  = nonzero elements of A stored row-wise in order *
!          ja(1:nnz) = corresponding column indices.                   *
!          ia(1:n+1) = pointer to beginning of each row in aa and ja.  *
!          here nnz = number of nonzero elements in A = ia(n+1)-ia(1)  *
!                                                                      *
! alu,jlu== A matrix stored in Modified Sparse Row format containing   *
!           the L and U factors, as computed by subroutine ilut.       *
!                                                                      *
! ju     == integer array of length n containing the pointers to       *
!           the beginning of each row of U in alu, jlu as computed     *
!           by subroutine ILUT.                                        *
!                                                                      *
! on return:                                                           *
!==========                                                            *
! sol   == contains an approximate solution (upon successful return).  *
! ierr  == integer. Error message with the following meaning.          *
!          ierr = 0 --> successful return.                             *
!          ierr = 1 --> convergence not achieved in itmax iterations.  *
!          ierr =-1 --> the initial guess seems to be the exact        *
!                       solution (initial residual computed was zero)  *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! work arrays:                                                         *
!=============                                                         *
! vv    == work array of length  n x (im+1) (used to store the Arnoli  *
!          basis)                                                      *
!----------------------------------------------------------------------*
! subroutines called :                                                 *
! amux   : SPARSKIT routine to do the matrix by vector multiplication  *
!          delivers y=Ax, given x  -- see SPARSKIT/BLASSM/amux         *
! lusol : combined forward and backward solves (Preconditioning ope.) *
! BLAS1  routines.                                                     *
!----------------------------------------------------------------------*
       parameter (kmax=50)
       real(c_double) hh(kmax+1,kmax), c(kmax), s(kmax), rs(kmax+1),t
!-------------------------------------------------------------
! arnoldi size should not exceed kmax=50 in this version..
! to reset modify paramter kmax accordingly.
!-------------------------------------------------------------
       data epsmac/1.d-16/
       n1 = n + 1
       its = 0

!TEMPORARY PRINTOUT
!$$$       iunit=41
!$$$       open(iunit)
!$$$       nnz=ia(n+1)-ia(n)
!$$$       nnz1=ju(n)-ju(1)
!$$$       write(iunit,*)'pgmres: n,nnz,nnz1,im,maxits,eps',
!$$$     +      n,nnz,nnz1,im,maxits,eps
!$$$       write(iunit,*)'rhs(1:n)'
!$$$       write(iunit,101) (rhs(i),i=1,n)
!$$$       write(iunit,*)'sol(1:n)'
!$$$       write(iunit,101) (sol(i),i=1,n)
!$$$       write(iunit,*)'aa(1:nnz)'
!$$$       write(iunit,101) (aa(i),i=1,nnz)
!$$$       write(iunit,*)'ja(1:nnz)'
!$$$       write(iunit,102) (ja(i),i=1,nnz)
!$$$       write(iunit,*)'ia(1:nnz)'
!$$$       write(iunit,102) (ia(i),i=1,n)
!$$$       write(iunit,*)'alu(1:nnz1)'
!$$$       write(iunit,101) (alu(i),i=1,nnz1)
!$$$       write(iunit,*)'jlu(1:nnz1)'
!$$$       write(iunit,102) (jlu(i),i=1,nnz1)
!$$$       write(iunit,*)'ju(1:n)'
!$$$       write(iunit,102) (ju(i),i=1,n)
!$$$       close(iunit)
!$$$ 101   format(10(1pe9.2))
!$$$ 102   format(10i9)


!-------------------------------------------------------------
! outer loop starts here..
!-------------- compute initial residual vector --------------
       call amux (n, sol, vv, aa, ja, ia)
       do 21 j=1,n
          vv(j,1) = rhs(j) - vv(j,1)
 21    continue
!-------------------------------------------------------------
 20    ro = dnrm2(n, vv, 1)
       if (iout .gt. 0 .and. its .eq. 0) &
            write(iout, 199) its, ro
       if (ro .eq. 0.0d0) goto 999
       t = 1.0d0/ ro
       do 210 j=1, n
          vv(j,1) = vv(j,1)*t
 210   continue
       if (its .eq. 0) eps1=eps*ro
!     ** initialize 1-st term  of rhs of hessenberg system..
       rs(1) = ro
       i = 0
 4     i=i+1
       its = its + 1
       i1 = i + 1
       call lusol (n, vv(1,i), rhs, alu, jlu, ju)
       call amux (n, rhs, vv(1,i1), aa, ja, ia)
!-----------------------------------------
!     modified gram - schmidt...
!-----------------------------------------
       do 55 j=1, i
          t = ddot(n, vv(1,j),1,vv(1,i1),1)
          hh(j,i) = t
          call daxpy(n, -t, vv(1,j), 1, vv(1,i1), 1)
 55    continue
       t = dnrm2(n, vv(1,i1), 1)
       hh(i1,i) = t
       if ( t .eq. 0.0d0) goto 58
       t = 1.0d0/t
       do 57  k=1,n
          vv(k,i1) = vv(k,i1)*t
 57    continue
!
!     done with modified gram schimd and arnoldi step..
!     now  update factorization of hh
!
 58    if (i .eq. 1) goto 121
!--------perfrom previous transformations  on i-th column of h
       do 66 k=2,i
          k1 = k-1
          t = hh(k1,i)
          hh(k1,i) = c(k1)*t + s(k1)*hh(k,i)
          hh(k,i) = -s(k1)*t + c(k1)*hh(k,i)
 66    continue
 121   gam = sqrt(hh(i,i)**2 + hh(i1,i)**2)
!
!     if gamma is zero then any small value will do...
!     will affect only residual estimate
!
       if (gam .eq. 0.0d0) gam = epsmac
!
!     get  next plane rotation
!
       c(i) = hh(i,i)/gam
       s(i) = hh(i1,i)/gam
       rs(i1) = -s(i)*rs(i)
       rs(i) =  c(i)*rs(i)
!
!     detrermine residual norm and test for convergence-
!
       hh(i,i) = c(i)*hh(i,i) + s(i)*hh(i1,i)
       ro = abs(rs(i1))
 131   format(1h ,2e14.4)
       if (iout .gt. 0) &
            write(iout, 199) its, ro
       if (i .lt. im .and. (ro .gt. eps1))  goto 4
!
!     now compute solution. first solve upper triangular system.
!
       rs(i) = rs(i)/hh(i,i)
       do 30 ii=2,i
          k=i-ii+1
          k1 = k+1
          t=rs(k)
          do 40 j=k1,i
             t = t-hh(k,j)*rs(j)
 40       continue
          rs(k) = t/hh(k,k)
 30    continue
!
!     form linear combination of v(*,i)'s to get solution
!
       t = rs(1)
       do 15 k=1, n
          rhs(k) = vv(k,1)*t
 15    continue
       do 16 j=2, i
          t = rs(j)
          do 161 k=1, n
             rhs(k) = rhs(k)+t*vv(k,j)
 161      continue
 16    continue
!
!     call preconditioner.
!
       call lusol (n, rhs, rhs, alu, jlu, ju)
       do 17 k=1, n
          sol(k) = sol(k) + rhs(k)
 17    continue
!
!     restart outer loop  when necessary
!
       if (ro .le. eps1) goto 990
       if (its .ge. maxits) goto 991
!
!     else compute residual vector and continue..
!
       do 24 j=1,i
          jj = i1-j+1
          rs(jj-1) = -s(jj-1)*rs(jj)
          rs(jj) = c(jj-1)*rs(jj)
 24    continue
       do 25  j=1,i1
          t = rs(j)
          if (j .eq. 1)  t = t-1.0d0
          call daxpy (n, t, vv(1,j), 1,  vv, 1)
 25    continue
 199   format('   its =', i4, ' res. norm =', d20.6)
!     restart outer loop.
       goto 20
 990   ierr = 0
       return
 991   ierr = 1
       return
 999   continue
       ierr = -1
       return
!-----------------end of pgmres ---------------------------------------
!-----------------------------------------------------------------------
       end
!-----------------------------------------------------------------------
        subroutine lusol(n, y, x, alu, jlu, ju)
        real(c_double) x(n), y(n), alu(*)
        integer n, jlu(*), ju(*)

!-----------------------------------------------------------------------
!
! This routine solves the system (LU) x = y,
! given an LU decomposition of a matrix stored in (alu, jlu, ju)
! modified sparse row format
!
!-----------------------------------------------------------------------
! on entry:
! n   = dimension of system
! y   = the right-hand-side vector
! alu, jlu, ju
!     = the LU matrix as provided from the ILU routines.
!
! on return
! x   = solution of LU x = y.
!-----------------------------------------------------------------------
!
! Note: routine is in place: call lusol (n, x, x, alu, jlu, ju)
!       will solve the system with rhs x and overwrite the result on x .
!
!-----------------------------------------------------------------------
! local variables
!
        integer i,k
!
! forward solve
!
        do 40 i = 1, n
           x(i) = y(i)
           do 41 k=jlu(i),ju(i)-1
              x(i) = x(i) - alu(k)* x(jlu(k))
 41        continue
 40     continue
!
!     backward solve.
!
        do 90 i = n, 1, -1
           do 91 k=ju(i),jlu(i+1)-1
              x(i) = x(i) - alu(k)*x(jlu(k))
 91        continue
           x(i) = alu(i)*x(i)
 90     continue
!
        return
!----------------end of lusol ------------------------------------------
!-----------------------------------------------------------------------
        end
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        subroutine qsplit(a,ind,n,ncut)
        real(c_double) a(n)
        integer ind(n), n, ncut
!-----------------------------------------------------------------------
!     does a quick-sort split of a real array.
!     on input a(1:n). is a real array
!     on output a(1:n) is permuted such that its elements satisfy:
!
!     abs(a(i)) .ge. abs(a(ncut)) for i .lt. ncut and
!     abs(a(i)) .le. abs(a(ncut)) for i .gt. ncut
!
!     ind(1:n) is an integer array which permuted in the same way as a(*).
!-----------------------------------------------------------------------
        real(c_double) tmp, abskey
        integer itmp, first, last
!-----
        first = 1
        last = n
        if (ncut .lt. first .or. ncut .gt. last) return
!
!     outer loop -- while mid .ne. ncut do
!
 1      mid = first
        abskey = abs(a(mid))
        do 2 j=first+1, last
           if (abs(a(j)) .gt. abskey) then
              mid = mid+1
!     interchange
              tmp = a(mid)
              itmp = ind(mid)
              a(mid) = a(j)
              ind(mid) = ind(j)
              a(j)  = tmp
              ind(j) = itmp
           endif
 2      continue
!
!     interchange
!
        tmp = a(mid)
        a(mid) = a(first)
        a(first)  = tmp
!
        itmp = ind(mid)
        ind(mid) = ind(first)
        ind(first) = itmp
!
!     test for while loop
!
        if (mid .eq. ncut) return
        if (mid .gt. ncut) then
           last = mid-1
        else
           first = mid+1
        endif
        goto 1
!----------------end-of-qsplit------------------------------------------
!-----------------------------------------------------------------------
        end


!----------------------------------------------------------------------c
!                          S P A R S K I T                             c
!----------------------------------------------------------------------c
!          BASIC MATRIX-VECTOR OPERATIONS - MATVEC MODULE              c
!         Matrix-vector Mulitiplications and Triang. Solves            c
!----------------------------------------------------------------------c
! contents: (as of Nov 18, 1991)                                       c
!----------                                                            c
! 1) Matrix-vector products:                                           c
!---------------------------                                           c
! amux  : A times a vector. Compressed Sparse Row (CSR) format.        c
! amuxms: A times a vector. Modified Compress Sparse Row format.       c
! atmux : Transp(A) times a vector. CSR format.                        c
! atmuxr: Transp(A) times a vector. CSR format. A rectangular.         c
! amuxe : A times a vector. Ellpack/Itpack (ELL) format.               c
! amuxd : A times a vector. Diagonal (DIA) format.                     c
! amuxj : A times a vector. Jagged Diagonal (JAD) format.              c
! vbrmv : Sparse matrix-full vector product, in VBR format             c
!                                                                      c
! 2) Triangular system solutions:                                      c
!-------------------------------                                       c
! lsol  : Unit Lower Triang. solve. Compressed Sparse Row (CSR) format.c
! ldsol : Lower Triang. solve.  Modified Sparse Row (MSR) format.      c
! lsolc : Unit Lower Triang. solve. Comp. Sparse Column (CSC) format.  c
! ldsolc: Lower Triang. solve. Modified Sparse Column (MSC) format.    c
! ldsoll: Lower Triang. solve with level scheduling. MSR format.       c
! usol  : Unit Upper Triang. solve. Compressed Sparse Row (CSR) format.c
! udsol : Upper Triang. solve.  Modified Sparse Row (MSR) format.      c
! usolc : Unit Upper Triang. solve. Comp. Sparse Column (CSC) format.  c
! udsolc: Upper Triang. solve.  Modified Sparse Column (MSC) format.   c
!----------------------------------------------------------------------c
! 1)     M A T R I X    B Y    V E C T O R     P R O D U C T S         c
!----------------------------------------------------------------------c
      subroutine amux (n, x, y, a,ja,ia)
      real(c_double)  x(*), y(*), a(*)
      integer n, ja(*), ia(*)
!-----------------------------------------------------------------------
!         A times a vector
!-----------------------------------------------------------------------
! multiplies a matrix by a vector using the dot product form
! Matrix A is stored in compressed sparse row storage.
!
! on entry:
!----------
! n     = row dimension of A
! x     = real array of length equal to the column dimension of
!         the A matrix.
! a, ja,
!    ia = input matrix in compressed sparse row format.
!
! on return:
!-----------
! y     = real array of length n, containing the product y=Ax
!
!-----------------------------------------------------------------------
! local variables
!
      real(c_double) t
      integer i, k
!-----------------------------------------------------------------------
      do 100 i = 1,n
!
!     compute the inner product of row i with vector x
!
         t = 0.0d0
         do 99 k=ia(i), ia(i+1)-1
            t = t + a(k)*x(ja(k))
 99      continue
!
!     store result in y(i)
!
         y(i) = t
 100  continue
!
      return
!---------end-of-amux---------------------------------------------------
!-----------------------------------------------------------------------
      end
!-----------------------------------------------------------------------



!----------------------------------------------------------------------c
!                          S P A R S K I T                             c
!----------------------------------------------------------------------c
!        BASIC LINEAR ALGEBRA FOR SPARSE MATRICES. BLASSM MODULE       c
!----------------------------------------------------------------------c
! amub   :   computes     C = A*B                                      c
! aplb   :   computes     C = A+B                                      c
! aplb1  :   computes     C = A+B  [Sorted version: A, B, C sorted]    c
! aplsb  :   computes     C = A + s B                                  c
! aplsb1 :   computes     C = A+sB  [Sorted version: A, B, C sorted]   c
! apmbt  :   Computes     C = A +/- transp(B)                          c
! aplsbt :   Computes     C = A + s * transp(B)                        c
! diamua :   Computes     C = Diag * A                                 c
! amudia :   Computes     C = A* Diag                                  c
! aplsca :   Computes     A:= A + s I    (s = scalar)                  c
! apldia :   Computes     C = A + Diag.                                c
!----------------------------------------------------------------------c
! Note: this module still incomplete.                                  c
!----------------------------------------------------------------------c

!-----------------------------------------------------------------------
      subroutine aplb (nrow,ncol,job,a,ja,ia,b,jb,ib, &
           c,jc,ic,nzmax,iw,ierr)
      real(c_double) a(*), b(*), c(*)
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(nrow+1),ic(nrow+1), &
           iw(ncol)
!-----------------------------------------------------------------------
! performs the matrix sum  C = A+B.
!-----------------------------------------------------------------------
! on entry:
! ---------
! nrow  = integer. The row dimension of A and B
! ncol  = integer. The column dimension of A and B.
! job   = integer. Job indicator. When job = 0, only the structure
!                  (i.e. the arrays jc, ic) is computed and the
!                  real values are ignored.
!
! a,
! ja,
! ia   = Matrix A in compressed sparse row format.
!
! b,
! jb,
! ib    =  Matrix B in compressed sparse row format.
!
! nzmax = integer. The  length of the arrays c and jc.
!         amub will stop if the result matrix C  has a number
!         of elements that exceeds exceeds nzmax. See ierr.
!
! on return:
!----------
! c,
! jc,
! ic    = resulting matrix C in compressed sparse row sparse format.
!
! ierr  = integer. serving as error message.
!         ierr = 0 means normal return,
!         ierr .gt. 0 means that amub stopped while computing the
!         i-th row  of C with i=ierr, because the number
!         of elements in C exceeds nzmax.
!
! work arrays:
!------------
! iw    = integer work array of length equal to the number of
!         columns in A.
!
!-----------------------------------------------------------------------
      logical values
      values = (job .ne. 0)
      ierr = 0
      len = 0
      ic(1) = 1
      do 1 j=1, ncol
         iw(j) = 0
 1    continue
!
      do 500 ii=1, nrow
!     row i
         do 200 ka=ia(ii), ia(ii+1)-1
            len = len+1
            jcol    = ja(ka)
            if (len .gt. nzmax) goto 999
            jc(len) = jcol
            if (values) c(len)  = a(ka)
            iw(jcol)= len
 200     continue
!
         do 300 kb=ib(ii),ib(ii+1)-1
            jcol = jb(kb)
            jpos = iw(jcol)
            if (jpos .eq. 0) then
               len = len+1
               if (len .gt. nzmax) goto 999
               jc(len) = jcol
               if (values) c(len)  = b(kb)
               iw(jcol)= len
            else
               if (values) c(jpos) = c(jpos) + b(kb)
            endif
 300     continue
         do 301 k=ic(ii), len
            iw(jc(k)) = 0
 301     continue
         ic(ii+1) = len+1
 500  continue
      return
 999  ierr = ii
      return
!------------end of aplb -----------------------------------------------
!-----------------------------------------------------------------------
      end
!-----------------------------------------------------------------------




!----------------------------------------------------------------------c
!                          S P A R S K I T                             c
!----------------------------------------------------------------------c
!                    FORMAT CONVERSION MODULE                          c
!----------------------------------------------------------------------c

!-----------------------------------------------------------------------
      subroutine bndcsr (n,abd,nabd,lowd,ml,mu,a,ja,ia,len,ierr)
      real(c_double) a(*),abd(nabd,*), t
      integer ia(n+1),ja(*)
!-----------------------------------------------------------------------
! Banded (Linpack ) format   to    Compressed Sparse Row  format.
!-----------------------------------------------------------------------
! on entry:
!----------
! n     = integer,the actual row dimension of the matrix.
!
! nabd  = first dimension of array abd.
!
! abd   = real array containing the values of the matrix stored in
!         banded form. The j-th column of abd contains the elements
!         of the j-th column of  the original matrix,comprised in the
!         band ( i in (j-ml,j+mu) ) with the lowest diagonal located
!         in row lowd (see below).
!
! lowd  = integer. this should be set to the row number in abd where
!         the lowest diagonal (leftmost) of A is located.
!         lowd should be s.t.  ( 1  .le.  lowd  .le. nabd).
!         The subroutines dgbco, ... of linpack use lowd=2*ml+mu+1.
!
! ml    = integer. equal to the bandwidth of the strict lower part of A
! mu    = integer. equal to the bandwidth of the strict upper part of A
!         thus the total bandwidth of A is ml+mu+1.
!         if ml+mu+1 is found to be larger than nabd then an error
!         message is set. see ierr.
!
! len   = integer. length of arrays a and ja. bndcsr will stop if the
!         length of the arrays a and ja is insufficient to store the
!         matrix. see ierr.
!
! on return:
!-----------
! a,
! ja,
! ia    = input matrix stored in compressed sparse row format.
!
! lowd  = if on entry lowd was zero then lowd is reset to the default
!         value ml+mu+l.
!
! ierr  = integer. used for error message output.
!         ierr .eq. 0 :means normal return
!         ierr .eq. -1 : means invalid value for lowd.
!         ierr .gt. 0 : means that there was not enough storage in a and ja
!         for storing the ourput matrix. The process ran out of space
!         (as indicated by len) while trying to fill row number ierr.
!         This should give an idea of much more storage might be required.
!         Moreover, the first irow-1 rows are correctly filled.
!
! notes:  the values in abd found to be equal to zero
! -----   (actual test: if (abd(...) .eq. 0.0d0) are removed.
!         The resulting may not be identical to a csr matrix
!         originally transformed to a bnd format.
!
!-----------------------------------------------------------------------
      ierr = 0
!-----------
      if (lowd .gt. nabd .or. lowd .le. 0) then
         ierr = -1
         return
      endif
!-----------
      ko = 1
      ia(1) = 1
      do 30 irow=1,n
!-----------------------------------------------------------------------
         i = lowd
          do  20 j=irow-ml,irow+mu
             if (j .le. 0 ) goto 19
             if (j .gt. n) goto 21
             t = abd(i,j)
             if (t .eq. 0.0d0) goto 19
             if (ko .gt. len) then
               ierr = irow
               return
            endif
            a(ko) = t
            ja(ko) = j
            ko = ko+1
 19         i = i-1
 20      continue
!     end for row irow
 21      ia(irow+1) = ko
 30   continue
      return
!------------- end of bndcsr -------------------------------------------
!-----------------------------------------------------------------------
      end
!-----------------------------------------------------------------------

      real(c_double) function dnrm2 ( n, dx, incx)
      integer          next
!BH100826  double precision   dx(1), cutlo, cuthi, hitest, sum, xmax,zero,one
      real(c_double)   dx(*), cutlo, cuthi, hitest, sum, xmax,zero,one
      data   zero, one /0.0d0, 1.0d0/
!
!     euclidean norm of the n-vector stored in dx() with storage
!     increment incx .
!     if    n .le. 0 return with result = 0.
!     if n .ge. 1 then incx must be .ge. 1
!
!           c.l.lawson, 1978 jan 08
!
!     four phase method     using two built-in constants that are
!     hopefully applicable to all machines.
!         cutlo = maximum of  dsqrt(u/eps)  over all known machines.
!         cuthi = minimum of  dsqrt(v)      over all known machines.
!     where
!         eps = smallest no. such that eps + 1. .gt. 1.
!         u   = smallest positive no.   (underflow limit)
!         v   = largest  no.            (overflow  limit)
!
!     brief outline of algorithm..
!
!     phase 1    scans zero components.
!     move to phase 2 when a component is nonzero and .le. cutlo
!     move to phase 3 when a component is .gt. cutlo
!     move to phase 4 when a component is .ge. cuthi/m
!     where m = n for x() real and m = 2*n for complex.
!
!     values for cutlo and cuthi..
!     from the environmental parameters listed in the imsl converter
!     document the limiting values are as follows..
!     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are
!                   univac and dec at 2**(-103)
!                   thus cutlo = 2**(-51) = 4.44089e-16
!     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec.
!                   thus cuthi = 2**(63.5) = 1.30438e19
!     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec.
!                   thus cutlo = 2**(-33.5) = 8.23181d-11
!     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19
!     data cutlo, cuthi / 8.232d-11,  1.304d19 /
!     data cutlo, cuthi / 4.441e-16,  1.304e19 /
      data cutlo, cuthi / 8.232d-11,  1.304d19 /
!
      cutlo = dsqrt(EPSILON(one))

      if(n .gt. 0) go to 10
         dnrm2  = zero
         go to 300
!
   10 assign 30 to next
      sum = zero
      nn = n * incx
!                                                 begin main loop
      i = 1
   20 continue   ! next i
      go to next,(30, 50, 70, 110)
   30 if( dabs(dx(i)) .gt. cutlo) go to 85
      assign 50 to next
      xmax = zero
!
!                        phase 1.  sum is zero
!
   50 if( dx(i) .eq. zero) go to 200
      if( dabs(dx(i)) .gt. cutlo) go to 85
!
!                                prepare for phase 2.
      assign 70 to next
      go to 105
!
!                                prepare for phase 4.
!
  100 i = j
      assign 110 to next
      sum = (sum / dx(i)) / dx(i)
  105 xmax = dabs(dx(i))
      go to 115
!
!                   phase 2.  sum is small.
!                             scale to avoid destructive underflow.
!
   70 if( dabs(dx(i)) .gt. cutlo ) go to 75
!
!                     common code for phases 2 and 4.
!                     in phase 4 sum is large.  scale to avoid overflow.
!
  110 if( dabs(dx(i)) .le. xmax ) go to 115
         sum = one + sum * (xmax / dx(i))**2
         xmax = dabs(dx(i))
         go to 200
!
  115 continue
      sum = sum + (dx(i)/xmax)**2
      go to 200
!
!
!                  prepare for phase 3.
!
   75 sum = (sum * xmax) * xmax
!
!
!     for real or d.p. set hitest = cuthi/n
!     for complex      set hitest = cuthi/(2*n)
!
   85 hitest = cuthi/float( n )
!
!                   phase 3.  sum is mid-range.  no scaling.
!
      do 95 j =i,nn,incx
      if(dabs(dx(j)) .ge. hitest) go to 100
   95    sum = sum + dx(j)**2
      dnrm2 = dsqrt( sum )
      go to 300
!
  200 continue
      i = i + incx
      if ( i .le. nn ) go to 20
!
!              end of main loop.
!
!              compute square root and adjust for scaling.
!
      dnrm2 = xmax * dsqrt(sum)
  300 continue
      return
      end
!----------------------------------------------------------------------



end module ilut_mod
