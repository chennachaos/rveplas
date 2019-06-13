CDOC BEGIN_SUBROUTINE RTSR
CDOC Matrix product s.Rt S R
CDOC
CDOC This routine performs the matrix product s Rt S R, where s is a
CDOC scalar, R a rectangular real matrix and S a square real matrix.
CDOC Rt denotes the transpose of R.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION AUXM   <  Auxiliary matrix used to store partial
CDOC C                          results of the calculation.
CDOC INTEGER          MODE   >  If set to 1, the argument Q
CDOC C                          returns the resulting matrix Rt S R.
CDOC C                          Otherwise, Rt S R is added to the input
CDOC C                          value of Q.
CDOC INTEGER          MROWQ  >  Dimensioning parameter: maximum
CDOC C                          dimension of the square matrix
CDOC C                          Q.
CDOC INTEGER          MROWR  >  Dimensioning parameter: maximum number
CDOC C                          of rows of R (same as the
CDOC C                          maximum dimension of square matrix S).
CDOC INTEGER          NCOLR  >  Number of columns of R.
CDOC INTEGER          NROWR  >  Number of rows of R.
CDOC DOUBLE_PRECISION Q      <> Matrix where results are stored.
CDOC DOUBLE_PRECISION R      >  Rectangular real matrix.
CDOC DOUBLE_PRECISION S      >  Square real matrix.
CDOC DOUBLE_PRECISION SCAL   >  Real scalar.
CDOC LOGICAL          UNSYM  >  Unsymmetry flag.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, August 1996: Initial coding
CHST
      SUBROUTINE RTSR
     1(   AUXM       ,MODE       ,MROWQ      ,MROWR      ,NCOLR      ,
     2    NROWR      ,Q          ,R          ,S          ,SCAL       ,
     3    UNSYM      )  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL UNSYM
      DIMENSION
     1    AUXM(NCOLR,NROWR)  ,Q(MROWQ,MROWQ)     ,R(MROWR,NCOLR)     ,
     2    S(MROWR,MROWR)
      DATA  R0   /0.0D0/
C***********************************************************************
C PERFORMS THE MATRIX PRODUCTS
C                                 T
C                    Q := SCAL * R  S R        (IF MODE=1)
C OR
C                                     T
C                    Q := Q + SCAL * R  S R    (OTHERWISE)
C
C WHERE 'R' IS A REAL RECTANGULAR MATRIX, 'S' A REAL SQUARE MATRIX
C AND 'SCAL' A SCALAR.
C***********************************************************************
      CALL RVZERO(AUXM,NCOLR*NROWR)
      DO 30 I=1,NCOLR
        DO 20 K=1,NROWR
          IF(R(K,I).NE.R0)THEN
            DO 10 J=1,NROWR
              AUXM(I,J)=AUXM(I,J)+SCAL*R(K,I)*S(K,J)
   10       CONTINUE
          ENDIF
   20   CONTINUE
   30 CONTINUE
C
      IF(MODE.EQ.1)THEN
        DO 50 I=1,NCOLR
          DO 40 J=1,NCOLR
            Q(I,J)=R0
   40     CONTINUE
   50   CONTINUE
      ENDIF
C
      IF(UNSYM)THEN
C Construct the whole matrix Q at once
        DO 80 J=1,NCOLR
          DO 70 K=1,NROWR
            IF(R(K,J).NE.R0)THEN
              DO 60 I=1,NCOLR
                Q(I,J)=Q(I,J)+AUXM(I,K)*R(K,J)
   60         CONTINUE
            ENDIF
   70     CONTINUE
   80   CONTINUE
      ELSE
C Construct the lower triangle of Q first
        DO 110 J=1,NCOLR
          DO 100 K=1,NROWR
            IF(R(K,J).NE.R0)THEN
              DO 90 I=J,NCOLR
                Q(I,J)=Q(I,J)+AUXM(I,K)*R(K,J)
   90         CONTINUE
            ENDIF
  100     CONTINUE
  110   CONTINUE
C and then assemble the upper triangle
        DO 130 I=1,NCOLR
          DO 120 J=I+1,NCOLR
            Q(I,J)=Q(J,I)
  120     CONTINUE
  130   CONTINUE
      ENDIF
C
      RETURN
      END
CDOC END_SUBROUTINE RTSR
