CDOC BEGIN_DOUBLE_PRECISION_FUNCTION PLFUN
CDOC Returns the value of a piece-wise linear scalar function
CDOC
CDOC This function returns the value of a piece-wise linear scalar
CDOC function F(X) defined by a set of NPOINT pairs (X,F(X)) passed in
CDOC the matrix argument XFX.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION X      >  Point at which the function will be
CDOC C                          evaluated.
CDOC INTEGER          NPOINT >  Number of points defining the piece-wise
CDOC C                          linear function.
CDOC DOUBLE_PRECISION XFX    >  Matrix (dimension 2*NPOINT)
CDOC C                          containing the pairs (x,f(x)) which
CDOC C                          define the piece-wise linear function.
CDOC C                          evaluated. Each column of XFX
CDOC C                          contains a pair (xi,f(xi)). The pairs
CDOC C                          supplied in XFX
CDOC C                          must be ordered such that the x's are
CDOC C                          monotonically increasing. That is, the
CDOC C                          x [XFX(1,i+1)] of a column i+1 must be
CDOC C                          greater than XFX(1,i) (x of column i).
CDOC C                          If X < XFX(1,1) the
CDOC C                          piece-wise linear function is assumed
CDOC C                          constant equal to XFX(1,1).
CDOC C                          If X > XFX(1,NPOINT) the
CDOC C                          piece-wise linear function is assumed
CDOC C                          constant equal to XFX(1,NPOINT).
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, August 1992: Initial coding
CHST
      DOUBLE PRECISION FUNCTION PLFUN(X, NPOINT, XFX)
C
      INTEGER NPOINT, I
      DOUBLE PRECISION X, XFX(2,*)
C***********************************************************************
C PIECEWISE LINEAR FUNCTION DEFINED BY A SET OF NPOINT PAIRS
C {X,F(X)} STORED IN THE MATRIX XFX (DIM. 2*NPOINT).
C***********************************************************************
      DO 100 I=1,NPOINT 
        IF (X.GE.XFX(1,I)) THEN
          GOTO 100
        ELSE  
          IF (I.EQ.1) THEN
C           -- x < x1 --> f(x)=f(x1) --- 
            PLFUN=XFX(2,1)
            GOTO 999
          ELSE
C           -- x(i-1) <= x < x(i) ---
            PLFUN=XFX(2,I-1)+(X-XFX(1,I-1))*
     1                 (XFX(2,I)-XFX(2,I-1))/
     2                 (XFX(1,I)-XFX(1,I-1))
            GOTO 999
          ENDIF
        ENDIF
 100  CONTINUE
C     ----  x >= x(npoint) --> f(x) = f(x(npoint))  ---
      PLFUN=XFX(2,NPOINT)
 999  CONTINUE
      RETURN
      END
CDOC END_DOUBLE_PRECISION_FUNCTION PLFUN
