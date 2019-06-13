CDOC BEGIN_SUBROUTINE RTSX
CDOC Matrix product s.Rt S X
CDOC
CDOC This routine performs the matrix product s Rt S X, where s is a
CDOC scalar, R and X rectangular real matrices of identical dimensions
CDOC and S a square real matrix.  Rt denotes the transpose of R.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION AUXM   <  Auxiliary matrix used to store partial
CDOC C                          results of the calculation.
CDOC INTEGER          MODE   >  If set to 1, the argument Q
CDOC C                          returns the resulting matrix Rt S R.
CDOC C                          Otherwise, Rt S R is added to the input
CDOC C                          value of Q.
CDOC INTEGER          MROWQ  >  Dimensioning parameter: maximum
CDOC C                          dimension of the square matrix Q.
CDOC INTEGER          MROWR  >  Dimensioning parameter: maximum number
CDOC C                          of rows of R (same as the
CDOC C                          maximum dimension of square matrix S).
CDOC INTEGER          NCOLR  >  Number of columns of R.
CDOC INTEGER          NROWR  >  Number of rows of R.
CDOC DOUBLE_PRECISION Q      <> Matrix where results are stored.
CDOC DOUBLE_PRECISION R      >  Rectangular real matrix.
CDOC DOUBLE_PRECISION S      >  Square real matrix.
CDOC DOUBLE_PRECISION X      >  Rectangular real matrix.
CDOC DOUBLE_PRECISION SCAL   >  Real scalar.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, September 1996: Initial coding
CHST
CHST E.de Souza Neto & F.M.A.Pires , April 2002:
CHST                          Bug fix in skipping zero multiplication
CHST
      SUBROUTINE RTSX
     1(   AUXM       ,MODE       ,MROWQ      ,MROWR      ,NCOLR      ,
     2    NROWR      ,Q          ,R          ,S          ,X          ,
     3    SCAL       )  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION
     1    AUXM(NCOLR,NROWR)  ,Q(MROWQ,MROWQ)     ,R(MROWR,NCOLR)     ,
     2    S(MROWR,MROWR)     ,X(MROWR,NCOLR)
      DATA  R0   /0.0D0/
C***********************************************************************
C PERFORMS THE MATRIX PRODUCTS
C                                  T
C                     Q := SCAL * R  S X        (IF MODE=1)
C OR
C                                      T
C                     Q := Q + SCAL * R  S X    (OTHERWISE)
C
C WHERE 'R' AND 'X' ARE REAL RECTANGULAR MATRICES OF IDENTICAL
C DIMENSIONS, 'S' A REAL SQUARE MATRIX AND 'SCAL' A SCALAR.
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
C Construct the matrix Q
      DO 80 J=1,NCOLR
        DO 70 K=1,NROWR
          IF(X(K,J).NE.R0)THEN
            DO 60 I=1,NCOLR
              Q(I,J)=Q(I,J)+AUXM(I,K)*X(K,J)
   60       CONTINUE
          ENDIF
   70   CONTINUE
   80 CONTINUE
C
      RETURN
      END
CDOC END_SUBROUTINE RTSX
