CDOC BEGIN_SUBROUTINE SOLQUA
CDOC Solves a quadratic equation:  a x**2 + b x + c = 0
CDOC
CDOC Given the coeficients a, b and c, this routine computes the real
CDOC roots of the associated quadratic equation, a x**2 + b x + c = 0.
CDOC The return values of the arguments ROOT1 and ROOT2 (the roots) are
CDOC set only if real roots exist.
CDOC If the equation admits only one real solution, the return value
CDOC of the logical argument ONEROO is set to .TRUE. (set to .FALSE.
CDOC otherwise).
CDOC If the equation admits two real roots, the return value
CDOC of the logical argument TWOROO is set to .TRUE.  (set to .FALSE.
CDOC otherwise).
CDOC Consequently, if the roots are complex or no roots/infinite
CDOC number of roots exist (the equation is ill-defined), the return
CDOC value of the logical arguments ONEROO and TWOROO is set to .FALSE.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION A      >  Coefficient of the quadratic term.
CDOC DOUBLE_PRECISION B      >  Coefficient of the linear term.
CDOC DOUBLE_PRECISION C      >  Coefficient of the constant term.
CDOC LOGICAL          ONEROO <  Logical flag. Set to .TRUE. if
CDOC C                          there is only one (real) root.
CDOC C                          Set to .FALSE. otherwise.
CDOC LOGICAL          TWOROO <  Logical flag. Set to .TRUE. if
CDOC C                          there are two distinct (real) root.
CDOC C                          Set to .FALSE. otherwise.
CDOC DOUBLE_PRECISION ROOT1  <  One of the roots (the only root if there
CDOC C                          is only one real root - not set if there
CDOC C                          are no real roots).
CDOC DOUBLE_PRECISION ROOT2  <  The other root (not set if there is only
CDOC C                          one real root or no real roots).
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, June 1998: Initial coding
CHST
      SUBROUTINE SOLQUA
     1(   A          ,B          ,C          ,ONEROO     ,TWOROO     ,
     2    ROOT1      ,ROOT2      )
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL  ONEROO ,TWOROO
      DATA
     1    R0   ,R1   ,R2   ,R4   ,SMALL /
     2    0.D0 ,1.0D0,2.0D0,4.0D0,1.D-12/
C***********************************************************************
C FINDS THE REAL ROOTS OF A QUADRATIC EQUATION:  A X**2 + B X + C = 0.
C
C REFERENCE:
C W.H.Press, S.A.Teukolsky, W.T.Vetterling and B.P.Flannery. Numerical
C recipes in FORTRAN. The art of scientific computing. 2nd Ed.,
C Cambridge Univ. Press, 1992. (Section 5.6)
C***********************************************************************
C Initialises logical flags
      ONEROO=.FALSE.
      TWOROO=.FALSE.
      IF(A.NE.R0)THEN
C The equation is non-linear in fact
C ----------------------------------
        IF(B.NE.R0)THEN
          SIGNB=B/ABS(B)
        ELSE
          SIGNB=R1
        ENDIF
        B2=B*B
        R4AC=R4*A*C
        SQUAR=B2-R4AC
        IF(SQUAR.GT.R0)THEN
C there are two distinct real roots: uses formula which minimises
C round-off errors when the coefficients A and/or C are small
          TWOROO=.TRUE.
          SQUAR=SQRT(SQUAR)
          Q=-(B+SIGNB*SQUAR)/R2
          ROOT1=Q/A
          ROOT2=C/Q
        ELSEIF(SQUAR.EQ.R0.OR.
     1         (SQUAR/DMAX1(B2,ABS(R4AC))+SMALL).GE.R0)THEN
C there is only one root
          ONEROO=.TRUE.
          ROOT1=-B/(R2*A)
        ENDIF
      ELSE
C The equation is linear
C ----------------------
        IF(B.NE.R0)THEN
C and well defined -> (only) one root exists
          ONEROO=.TRUE.
          ROOT1=-C/B
        ENDIF
      ENDIF
      RETURN
      END
CDOC END_SUBROUTINE SOLQUA
