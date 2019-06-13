CDOC BEGIN_SUBROUTINE RVZERO
CDOC Zero a double precision array
CDOC
CDOC This routine sets to zero the N components of the double precision
CDOC array argument V.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION V      <  Zeroed double precision array.
CDOC INTEGER          N      >  Dimension of V.
CDOC END_PARAMETERS
CDOC
      SUBROUTINE RVZERO
     1(   V          ,N          )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V(N)
      DATA R0/0.0D0/
C***********************************************************************
C INITIALISES TO ZERO A DOUBLE PRECISION ARRAY OF DIMENSION N
C***********************************************************************
      DO 10 I=1,N
        V(I)=R0
   10 CONTINUE
      RETURN
      END
CDOC END_SUBROUTINE RVZERO
