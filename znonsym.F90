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

module znonsym_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double
  use cqlconf_mod, only : setup0

  !---END USE

!
!
!

contains

      subroutine nonsym(A,X,C,N,MLEFT,MRIGHT,EPS,NCOND)
      implicit integer (i-n), real(c_double) (a-h,o-z)
!     *************************************************
!
!     5.3  SOLVES A REAL VALUED NONSYMMETRIC LINEAR SYSTEM  A . X  =  C
!
!     VERSION 1               APRIL 1988       KA        LAUSANNE
!     VERSION 2               MAI 1992         AJ        LAUSANNE
!
!-----------------------------------------------------------------------
      DIMENSION   A((MLEFT+MRIGHT+1)*N),    X(N),    C(N)
!MPLX COMPLEX     A,       X,       C,        ZPIVOT,   ZSUM,     ZTOP
      real(c_double)        A,       X,       C,        ZPIVOT,   ZSUM,     ZTOP
      DATA        IMESS  / 0 /
!-----------------------------------------------------------------------
!
!     A IS A BANDMATRIX OF RANK N WITH MLEFT/MRIGHT OFF-DIAGONAL
!     ELEMENTS TO THE LEFT/RIGHT.
!     :A: AND :C: ARE DESTROYED BY :NONCYM:
!     AT THE RETURN :C: CONTAINS THE SOLUTION VECTOR, :X: THE SQUARE
!     ROOT OF THE DIAGONAL ELEMENTS
!
!     THE ELEMENTS ARE HORIZONTALLY NUMBERED ROW AFTER ROW
!
!     FOR EASY NUMBERING, ZERO ELEMENTS HAVE BEEN INTRODUCED IN THE
!     UPPER LEFT-HAND CORNER AND IN THE LOWER RIGHT-HAND CORNER.
!     THESE ELEMENTS
!     M U S T  B E  S E T  T O  Z E R O
!     BEFORE CALLING :NONCYM:
!
!     EXAMPLE FOR NUMBERING   MLEFT=2,  MRIGHT=3
!     -------
!
!     A(1)   A(2) I A(3)   A(4)   A(5)   A(6)
!     A(7) I A(8)   A(9)   A(10)  A(11)  A(12)
!     I                 .
!     I                        .
!
!     A PIVOT IS CONSIDERED TO BE BAD IF IT IS BY A FACTOR :EPS: SMALLER
!     THAN ITS OFF-DIAGONAL ELEMENTS. IF A BAD PIVOT IS ENCOUNTERED
!     A MESSAGE IS ISSUED (AT MOST TEN TIMES IN A RUN).
!     THE FLAG :NCOND: IS SET TO -1, OTHERWISE 0.
!
!-----------------------------------------------------------------------
!L    0.        INITIALIZATION
!
      MOFFDI=MLEFT+MRIGHT
      MBAND=MOFFDI+1
      RATIO=0.
      NCOND=0
!     PRELIMINARY CHECK OF PIVOTS
      IPIV1=MLEFT+1
      IPIV2=IPIV1+(N-1)*MBAND
      DO 10 J=IPIV1,IPIV2,MBAND
!MPLX IF(CABS(A(J)) .EQ. 0.) GO TO 510
        IF( ABS(A(J)) .EQ. 0.) GO TO 510
 10   CONTINUE
!
!-----------------------------------------------------------------------
!L    1.        PRECONDITIONNING: GET ONES ON THE DIAGONAL
!
      DO 130 JP=1,N
        JPABS=(MLEFT+1)+(JP-1)*MBAND
!MPLX X(JP)=SQRT(CABS(A(JPABS)))
        X(JP)=SQRT(ABS(A(JPABS)))
!
!     DIVIDE THE EQUATIONS BY THE SQUARE ROOT OF THE PIVOT
        DO 110 J=-MLEFT,MRIGHT
          A(JPABS+J)=A(JPABS+J)/X(JP)
 110    CONTINUE
        C(JP)=C(JP)/X(JP)
!
!     CHANGE VARIABLES, DIVIDE COLUMN BY SQUARE ROOT OF THE PIVOT
        DO 120 JA=JPABS-MRIGHT*MOFFDI,JPABS+MLEFT*MOFFDI,MOFFDI
          IF (JA.GE.IPIV1 .AND. JA.LE.IPIV2) A(JA)=A(JA)/X(JP)
 120    CONTINUE
 130  CONTINUE
!
!-----------------------------------------------------------------------
!L    2.        CONSTRUCT UPPER TRIANGULAR MATRIX
!
!     INITIALIZATION OF VERTICAL COUNTER
      IDOWN=MLEFT
!     PIVOT INDEX FOR FIRST INCOMPLETE COLUMN
      INCMPL=(N-MLEFT)*MBAND+MLEFT+1
!     FIRST PIVOT AND LAST BUT ONE
      IPIV1=MLEFT+1
      IPIV2=IPIV1+(N-2)*MBAND
!
      DO 220 JPIVOT=IPIV1,IPIV2,MBAND
!     COLUMN LENGTH
        IF(JPIVOT.GE.INCMPL) IDOWN=IDOWN-1
        ZPIVOT=A(JPIVOT)
!     FIRST AND LAST HORIZONTAL ELEMENT INDEX
        IHOR1=JPIVOT+1
        IHOR2=JPIVOT+MRIGHT
!     CHECK THE PIVOT
        ZMAX=0.
        DO 205 JHORIZ=IHOR1,IHOR2
!MPLX ZMAX=AMAX1(ZMAX,CABS(A(JHORIZ)))
!990131          ZMAX=AMAX1(ZMAX, ABS(A(JHORIZ)))
          ZMAX=MAX(ZMAX, ABS(A(JHORIZ)))
 205    CONTINUE
!MPLX IF(CABS(ZPIVOT).EQ.0.) GO TO 510
        IF( ABS(ZPIVOT).EQ.0.) GO TO 510
!MPLX RATIO=AMAX1(RATIO,ZMAX/CABS(ZPIVOT))
!990131        RATIO=AMAX1(RATIO,ZMAX/ ABS(ZPIVOT))
        RATIO=MAX(RATIO,ZMAX/ ABS(ZPIVOT))
!
        DO 210 JHORIZ=IHOR1,IHOR2
          ZTOP=-A(JHORIZ)/ZPIVOT
          A(JHORIZ)=ZTOP
!     INITIALIZATION OF ELEMENT INDEX
          IELEM=JHORIZ
!     INITIALIZATION OF RECTANGULAR RULE CORNER ELEMENT
          ICORN=JPIVOT
!     LOOP DOWN THE COLUMN
          DO 211 JDOWN=1,IDOWN
            IELEM=IELEM+MOFFDI
            ICORN=ICORN+MOFFDI
            A(IELEM)=A(IELEM)+ZTOP*A(ICORN)
 211      CONTINUE
 210    CONTINUE
!     TREAT THE CONSTANTS
!     INDEX OF THE FIRST ONE
        ICONST=JPIVOT/MBAND+1
        ZTOP=-C(ICONST)/ZPIVOT
        C(ICONST)=ZTOP
!     INITIALIZATION OF RECTANGULAR RULE CORNER ELEMENT
        ICORN=JPIVOT
!     LOOP DOWN THE CONSTANTS
        DO 221 JDOWN=1,IDOWN
          ICONST=ICONST+1
          ICORN=ICORN+MOFFDI
          C(ICONST)=C(ICONST)+ZTOP*A(ICORN)
 221    CONTINUE
 220  CONTINUE
!
!-----------------------------------------------------------------------
!L    3.        BACKSUBSTITUTION
!
 300  CONTINUE
!     INITIALIZATION OF SOLUTION INDEX
      ISOLUT=N
!     LAST PIVOT
      JPIVOT=IPIV2+MBAND
!     CHECK LAST PIVOT
!MPLX IF(CABS(A(JPIVOT)).EQ.0.) RATIO=1.E100
!%OS  IF( ABS(A(JPIVOT)).EQ.0.) RATIO=1.E100
      IF( ABS(A(JPIVOT)).EQ.0.) go to 510
!     LAST UNKNOWN
      C(ISOLUT)=C(ISOLUT)/A(JPIVOT)
!     INITIALIZATION OF HORIZONTAL RANGE
      IHOR1=JPIVOT+1
      IHOR2=JPIVOT+MRIGHT
!
      INM1=N-1
      DO 320 J=1,INM1
        ISOLUT=ISOLUT-1
        IHOR1=IHOR1-MBAND
        IHOR2=IHOR2-MBAND
        ZSUM=-C(ISOLUT)
        IKNOWN=ISOLUT
        DO 310 JHORIZ=IHOR1,IHOR2
          IKNOWN=IKNOWN+1
          IF(IKNOWN.GT.N) GO TO 310
          ZSUM=ZSUM+A(JHORIZ)*C(IKNOWN)
 310    CONTINUE
        C(ISOLUT)=ZSUM
 320  CONTINUE
!
!-----------------------------------------------------------------------
!L    4.        PRECONDITIONNING: BACK TO ORIGINAL VARIABLES
!
      DO 410 J=1,N
        C(J)=C(J)/X(J)
 410  CONTINUE
!
!-----------------------------------------------------------------------
!L    5.        HAVE WE ENCOUNTERED A BAD PIVOT ?
!
!%OS  IF(1./RATIO .GT. EPS) RETURN
      IF(RATIO .LT. 1./EPS) RETURN
      NCOND=-1
!     COUNT THE NUMBER OF MESSAGES ISSUED
      IMESS=IMESS+1
      IF(IMESS.LE.10) PRINT 9500, RATIO
      RETURN
!
!     ZERO PIVOT ENCOUNTERED
 510  CONTINUE
      if(setup0%verbose>0) PRINT 9510
      STOP 'ZERO PIVOT IN NONCYM'
!
!-----------------------------------------------------------------------
!L    9.        FORMATS
!

 9500 FORMAT(/////20X,10(1H*),21H  MESSAGE FROM NONCYM,2X,10(1H*)/// &
        20X,'BAD PIVOT ENCOUNTERED, RATIO: ',1PE12.3/// &
        20X,43(1H*)///)
 9510 FORMAT(/////20X,10(1H*),'   ZERO PIVOT IN NONCYM'///)
!
      END subroutine nonsym

end module znonsym_mod
