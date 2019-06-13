CDOC BEGIN_SUBROUTINE SYMTA
CDOC Converts the matrix representation (2-index array) of a symmetric 
CDOC second order tensor to its array representation.
CDOC Implemented for both 2-D and 3-D tensors.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION AMATX  >  Matrix representation of a symmetric 
CDOC                            second order tensor.
CDOC DOUBLE_PRECISION A      <  Array representation of the given 
CDOC C                          symmetric second order tensor.
CDOC INTEGER          NDIM   >  Number of spatial dimensions of tensor.
CDOC LOGICAL          OUTOFP >  .TRUE. if the out-of-plane component 
CDOC C                          (3,3) of a 2-D tensor is required.
CDOC END_PARAMETERS
CHST 
CHST D. de Bortoli  , March  2015: Initial coding
CHST
      SUBROUTINE SYMTA
     1(   AMATX       ,A       ,NDIM       ,OUTOFP      )
      IMPLICIT NONE
C Arguments
      DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN)  :: AMATX
      DOUBLE PRECISION, DIMENSION(6)  , INTENT(OUT) :: A
      INTEGER                         , INTENT(IN)  :: NDIM
      LOGICAL                         , INTENT(IN)  :: OUTOFP
C
      DOUBLE PRECISION, PARAMETER :: R0=0.0D0
C***********************************************************************
C RE-ARRANGES THE MATRIX REPRESENTATION OF A SYMMETRIC SECOND ORDER
C TENSOR TO ITS ARRAY REPRESENTATION (COMPONENT ORDERING (11,22,12,33) 
C IN 2-D, OR (11,22,33,12,23,13) IN 3-D).
C
C REFERENCE: Section D.1
C***********************************************************************
C
      IF(NDIM==2)THEN
        A(1)=AMATX(1,1)
        A(2)=AMATX(2,2)
        A(3)=AMATX(1,2)
C out-of-plane component
        IF(OUTOFP)THEN
          A(4)=AMATX(3,3)
        ELSE
          A(4)=R0
        ENDIF
        A(5:6)=R0
      ELSEIF(NDIM==3)THEN
        A(1)=AMATX(1,1)
        A(2)=AMATX(2,2)
        A(3)=AMATX(3,3)
        A(4)=AMATX(1,2)
        A(5)=AMATX(2,3)
        A(6)=AMATX(1,3)
      ELSE
        CALL ERRPRT('EI0068')
      ENDIF
C
      END
CDOC END_SUBROUTINE SYMTA
      