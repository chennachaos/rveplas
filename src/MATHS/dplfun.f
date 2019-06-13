CDOC BEGIN_DOUBLE_PRECISION_FUNCTION DPLFUN
CDOC Returns the derivative of piece-wise linear scalar function
CDOC
CDOC This procedure returns the derivative of the piece-wise linear
CDOC scalar function of procedure PLFUN.
CDOC The piece-wise linear function F(X) is defined by a set
CDOC of NPOINT pairs (X,F(X)) passed in the matrix argument XFX.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION X      >  Point at which the derivative will be
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
CDOC C                          If X > XFX(1,NPOINT) the piece-wise
CDOC C                          linear function is assumed constant
CDOC C                          equal to XFX(1,NPOINT).
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, August 1992: Initial coding
CHST
      DOUBLE PRECISION FUNCTION DPLFUN(X, NPOINT, XFX)
C
      INTEGER NPOINT, I
      DOUBLE PRECISION X, XFX(2,NPOINT), R0
      DATA R0 / 0.0D0 /
C***********************************************************************
C DERIVATIVE OF THE PIECEWISE LINEAR FUNCTION 'PLFUN' DEFINED BY A SET
C OF NPOINT PAIRS {X,F(X)} STORED IN THE MATRIX XFX (DIM. 2*NPOINT).
C***********************************************************************
      DO 100 I=1,NPOINT 
        IF (X.GE.XFX(1,I)) THEN
          GOTO 100
        ELSE
          IF (I.EQ.1) THEN
C           -- x < x1   --> f(x)=f(x1) --> df(x)/dx=0 --- 
            DPLFUN=R0
            GOTO 999
          ELSE
C           -- x(i-1) <= x < x(i) ---
            DPLFUN=(XFX(2,I)-XFX(2,I-1))/
     1             (XFX(1,I)-XFX(1,I-1))
            GOTO 999
          ENDIF
        ENDIF
 100  CONTINUE
C     ---- x >= x(npoint) --> f(x) = f(x(npoint)) --> df/dx=0 ---
      DPLFUN=R0 
 999  CONTINUE
      RETURN
      END
CDOC END_DOUBLE_PRECISION_FUNCTION DPLFUN
