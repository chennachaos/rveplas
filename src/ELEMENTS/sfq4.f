CDOC BEGIN_SUBROUTINE SFQ4
CDOC Shape function and derivatives for element type QUAD4
CDOC
CDOC This routine computes the shape functions and shape function
CDOC derivatives for the element type QUAD4: Standard isoparametric
CDOC 4-noded quadrilateral for plane strain, plane stress and
CDOC axisymmetric analysis.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DERIV  <  Array of shape function derivatives.
CDOC DOUBLE_PRECISION ETASP  >  Isoparametric coordinate ETA of the
CDOC C                          point where the shape functions or their
CDOC C                          derivatives are to be evaluated.
CDOC DOUBLE_PRECISION EXISP  >  Isoparametric coordinate XI of the
CDOC C                          point where the shape functions or their
CDOC C                          derivatives are to be evaluated.
CDOC INTEGER          IBOUND >  Boundary interpolation flag. Entry
CDOC C                          must be 0 for domain interpolation.
CDOC C                          Boundary interpolation is assumed
CDOC C                          otherwise.
CDOC INTEGER          MDIME  >  Maximum permissible number of spatial
CDOC C                          dimensions. Used here only for array
CDOC C                          dimensioning.
CDOC DOUBLE_PRECISION SHAPE  <  Array of shape function values.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, September 1996: Initial coding
CHST
      SUBROUTINE SFQ4
     1(   DERIV      ,ETASP      ,EXISP      ,IBOUND     ,MDIME      ,
     2    SHAPE      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION
     1    DERIV(MDIME,*)     ,SHAPE(*)
      DATA RP25 ,RP5  ,R1   /0.25D0,0.5D0,1.0D0/
C***********************************************************************
C COMPUTES SHAPE FUNCTIONS AND SHAPE FUNCTION DERIVATIVES FOR
C ELEMENT 'QUAD_4':
C                         4         3
C                          o-------o
C                          |       |     STANDARD ISOPARAMETRIC
C                          |       |     BI-LINEAR 4-NODE QUADRILATERAL 
C                          |       |
C                          o-------o
C                         1         2
C
C REFERENCE: Expression (4.42)
C***********************************************************************
      IF(IBOUND.EQ.0)THEN
C Shape functions and derivatives on element DOMAIN
C -------------------------------------------------
        S=EXISP
        T=ETASP
        ST=S*T
C Shape functions
        SHAPE(1)=(R1-T-S+ST)*RP25
        SHAPE(2)=(R1-T+S-ST)*RP25
        SHAPE(3)=(R1+T+S+ST)*RP25
        SHAPE(4)=(R1+T-S-ST)*RP25
C Shape function derivatives
        DERIV(1,1)=(-R1+T)*RP25
        DERIV(1,2)=(+R1-T)*RP25
        DERIV(1,3)=(+R1+T)*RP25
        DERIV(1,4)=(-R1-T)*RP25
        DERIV(2,1)=(-R1+S)*RP25
        DERIV(2,2)=(-R1-S)*RP25
        DERIV(2,3)=(+R1+S)*RP25
        DERIV(2,4)=(+R1-S)*RP25
      ELSE
C Shape function and derivatives on element BOUNDARY (1-D)
C --------------------------------------------------------
        S=EXISP
C Shape functions
        SHAPE(1)=(-S+R1)*RP5
        SHAPE(2)=(+S+R1)*RP5
C Shape functions derivatives
        DERIV(1,1)=-RP5
        DERIV(1,2)=RP5
C
      ENDIF
C
      RETURN
      END
CDOC END_SUBROUTINE SFQ4
