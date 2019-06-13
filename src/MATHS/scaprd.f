CDOC BEGIN_DOUBLE_PRECISION_FUNCTION SCAPRD
CDOC Scalar product of double precision vectors
CDOC
CDOC This function returns the scalar product between its two double
CDOC precision vector arguments U and V.
CDOC
CDOC BEGIN_PARAMETERS 
CDOC DOUBLE_PRECISION U      >  Array of components of a double
CDOC C                          precision vector.
CDOC DOUBLE_PRECISION V      >  Array of components of a double
CDOC C                          precision vector.
CDOC INTEGER          N      >  Dimension of U and V.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, May 1996: Initial coding
CHST
      DOUBLE PRECISION FUNCTION SCAPRD(U  ,V  ,N  ) 
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION U(N), V(N)
      DATA  R0  / 0.0D0 /
C***********************************************************************
C SCALAR PRODUCT OF DOUBLE PRECISION VECTORS U AND V OF DIMENSION N
C***********************************************************************
      SCAPRD=R0
      DO 10 I=1,N
        SCAPRD=SCAPRD+U(I)*V(I)
  10  CONTINUE
      RETURN
      END
CDOC END_DOUBLE_PRECISION_FUNCTION SCAPRD
