CDOC BEGIN_SUBROUTINE ATASYM
CDOC Converts the array representation of an asymmetric second order
CDOC tensor to its matrix representation (2-index array).
CDOC Implemented for both 2-D and 3-D tensors.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION A      >  Array representation of a asymmetric 
CDOC C                          second order tensor.
CDOC DOUBLE_PRECISION AMATX  <  Matrix representation of the given 
CDOC C                          asymmetric second order tensor.
CDOC INTEGER          NDIM   >  Number of spatial dimensions of tensor.
CDOC LOGICAL          OUTOFP >  .TRUE. if the out-of-plane component 
CDOC C                          (3,3) of a 2-D tensor is required.
CDOC END_PARAMETERS
CHST 
CHST D. de Bortoli  , March  2015: Initial coding
CHST
      SUBROUTINE ATASYM
     1(   A       ,AMATX       ,NDIM       ,OUTOFP      )
      IMPLICIT NONE
C Arguments
      DOUBLE PRECISION, DIMENSION(9)  , INTENT(IN)  :: A
      DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: AMATX
      INTEGER                         , INTENT(IN)  :: NDIM
      LOGICAL                         , INTENT(IN)  :: OUTOFP
C
      DOUBLE PRECISION, PARAMETER :: R0=0.0D0
C***********************************************************************
C RE-ARRANGES THE ARRAY REPRESENTATION OF AN ASYMMETRIC SECOND ORDER
C TENSOR, FROM COMPONENT ORDERING (11,21,12,22,33) IN 2-D, OR
C (11,21,31,12,22,32,13,23,33) IN 3-D, INTO ITS MATRIX REPRESENTATION.
C
C REFERENCE: Section D.2.1
C***********************************************************************
C
      IF(NDIM==2)THEN
        AMATX(1,1)=A(1)
        AMATX(2,1)=A(2)
        AMATX(1,2)=A(3)
        AMATX(2,2)=A(4)
        AMATX(1,3)=R0
        AMATX(2,3)=R0
        AMATX(3,3)=R0
        AMATX(3,1)=R0
        AMATX(3,2)=R0
C out-of-plane component
        IF(OUTOFP)THEN
          AMATX(3,3)=A(5)
        ELSE
          AMATX(3,3)=R0
        ENDIF
      ELSEIF(NDIM==3)THEN
        AMATX(1,1)=A(1)
        AMATX(2,1)=A(2)
        AMATX(3,1)=A(3)
        AMATX(1,2)=A(4)
        AMATX(2,2)=A(5)
        AMATX(3,2)=A(6)
        AMATX(1,3)=A(7)
        AMATX(2,3)=A(8)
        AMATX(3,3)=A(9)
      ELSE
        CALL ERRPRT('EI0071')
      ENDIF
C
      END
CDOC END_SUBROUTINE ATASYM
      