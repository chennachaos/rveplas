CDOC BEGIN_SUBROUTINE RVSCAL
CDOC Multiplies a double precision vector by a scalar.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION V      <> Double precision vector.
CDOC INTEGER          N      >  Dimension of V.
CDOC DOUBLE_PRECISION SCAL   >  Scalar by which V will be multiplied.
CDOC END_PARAMETERS
CDOC
      SUBROUTINE RVSCAL
     1(   V          ,N          ,SCAL       )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V(N)
C***********************************************************************
C MULTIPLIES THE DOUBLE PRECISION VECTOR 'V', OF DIMENSION 'N',
C BY THE SCALAR 'SCAL'
C***********************************************************************
      DO 10 I=1,N
        V(I)=SCAL*V(I)
   10 CONTINUE
      RETURN
      END
CDOC END_SUBROUTINE RVSCAL
