CDOC BEGIN_DOUBLE_PRECISION_FUNCTION DETM23
CDOC Calculates the determinant of a 2x2 or 3x3 matrix.
CDOC
CDOC BEGIN_PARAMETERS
CDOC INTEGER          DIMA   >  Size of square matrix AMATX.
CDOC DOUBLE_PRECISION AMATX  >  Matrix whose determinant is to be 
CDOC C                          calculated.
CDOC INTEGER          NDIM   >  Part of AMATX used for the determinant
CDOC C                          calculation: if NDIM=2, the 2x2 part of
CDOC C                          AMATX is used in the determinant
CDOC C                          calculation; if NDIM=3, the 3x3 part is
CDOC C                          is used instead.
CDOC LOGICAL          OUTOFP >  .TRUE. if the out-of-plane component 
CDOC C                          (3,3) of a 2x2 matrix is used in the
CDOC C                          determinant calculation.
CDOC END_PARAMETERS
CHST 
CHST D. de Bortoli  , March  2015: Initial coding
CHST
      DOUBLE PRECISION FUNCTION DETM23
     1(   DIMA      ,AMATX       ,NDIM       ,OUTOFP      )
      IMPLICIT NONE
C Arguments
      INTEGER                               , INTENT(IN)  :: DIMA
      DOUBLE PRECISION, DIMENSION(DIMA,DIMA), INTENT(IN)  :: AMATX
      INTEGER                               , INTENT(IN)  :: NDIM
      LOGICAL                               , INTENT(IN)  :: OUTOFP
C
C***********************************************************************
C CALCULATES THE DETERMINANT OF A 2X2 OR 3X3 MATRIX.
C***********************************************************************
C
      IF(NDIM > DIMA)THEN
        CALL ERRPRT('EI0074')
      ENDIF
C 2x2 determinant required
      IF(NDIM==2)THEN
        DETM23=AMATX(1,1)*AMATX(2,2)
     1        -AMATX(1,2)*AMATX(2,1)
C out-of-plane component of 2x2 determinant
        IF(OUTOFP)THEN
          DETM23=DETM23*AMATX(3,3)
        ENDIF
C 3x3 determinant required
      ELSEIF(NDIM==3)THEN
        DETM23=AMATX(1,1)*AMATX(2,2)*AMATX(3,3)
     1        +AMATX(1,2)*AMATX(2,3)*AMATX(3,1)
     2        +AMATX(1,3)*AMATX(2,1)*AMATX(3,2)
     3        -AMATX(1,2)*AMATX(2,1)*AMATX(3,3)
     4        -AMATX(1,1)*AMATX(2,3)*AMATX(3,2)
     5        -AMATX(1,3)*AMATX(2,2)*AMATX(3,1)
      ELSE
        CALL ERRPRT('EI0075')
      ENDIF
C
      END FUNCTION DETM23
CDOC END_DOUBLE_PRECISION_FUNCTION DETM23