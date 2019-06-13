CDOC BEGIN_SUBROUTINE RTV
CDOC Matrix-vectot product s.Rt V
CDOC
CDOC This routine performs the matrix-vector product s Rt V, where s is
CDOC a scalar, R a real rectangular matrix and v a real vector.
CDOC Rt denotes the transpose of R.
CDOC
CDOC BEGIN_PARAMETERS
CDOC INTEGER          MODE   >  If set to 1, the argument V
CDOC C                          returns the resulting vector s Rt V.
CDOC C                          Otherwise, s Rt v is added to the input
CDOC C                          value of P.
CDOC INTEGER          MROWR  >  Dimensioning parameter: maximum number
CDOC C                          of rows of R.
CDOC INTEGER          NCOLR  >  Number of columns of R.
CDOC INTEGER          NROWR  >  Number of rows of R.
CDOC DOUBLE_PRECISION P      <> Vector where results are stored.
CDOC DOUBLE_PRECISION R      >  Rectangular real matrix.
CDOC DOUBLE_PRECISION V      >  Real vector.
CDOC DOUBLE_PRECISION SCAL   >  Real scalar.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, August 1996: Initial coding
CHST
      SUBROUTINE RTV
     1(   MODE       ,MROWR      ,NCOLR      ,NROWR      ,P          ,
     2    R          ,V          ,SCAL       )  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION
     1    P(NCOLR)           ,R(MROWR,NCOLR)     ,V(NROWR)
      DATA  R0   /0.0D0/
C***********************************************************************
C PERFORMS THE PRODUCT
C                                  T
C                     P := SCAL * R  V          (IF MODE=1)
C OR
C                                      T
C                     P := P + SCAL * R  V      (OTHERWISE)
C
C WHERE 'R' IS A REAL RECTANGULAR MATRIX, 'V' A REAL VECTOR AND
C 'SCAL' A SCALAR.
C***********************************************************************
      IF(MODE.EQ.1)CALL RVZERO(P,NCOLR)
      DO 30 I=1,NCOLR
        DO 20 J=1,NROWR
          IF(R(J,I).NE.R0)THEN
            P(I)=P(I)+SCAL*R(J,I)*V(J)
          ENDIF
   20   CONTINUE
   30 CONTINUE
C
      RETURN
      END
CDOC END_SUBROUTINE RTV
