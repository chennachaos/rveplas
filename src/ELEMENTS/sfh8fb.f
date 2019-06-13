CDOC BEGIN_SUBROUTINE SFH8FB
CDOC Shape function and derivatives for element type HEXA_8_FBAR
CDOC
CDOC This routine computes the shape functions and shape function
CDOC derivatives for the element type HEXA_8_FBAR: F-bar isoparametric
CDOC 8-noded hexahedron for three-dimensional analysis.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DERIV  <  Array of shape function derivatives.
CDOC DOUBLE_PRECISION ZETASP >  Isoparametric coordinate ZETA of the
CDOC C                          point where the shape functions or their
CDOC C                          derivatives are to be evaluated.
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
CHST D. de Bortoli, April 2015: Initial coding
CHST
      SUBROUTINE SFH8FB
     1(   DERIV      ,ZETASP      ,ETASP      ,EXISP      ,IBOUND     ,
     2    MDIME      ,SHAPE      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION
     1    DERIV(MDIME,*)     ,SHAPE(*)
      DATA RP25, RP125, R1   /0.25D0,0.125D0,1.0D0/
C***********************************************************************
C COMPUTES SHAPE FUNCTIONS AND SHAPE FUNCTION DERIVATIVES FOR
C ELEMENT 'HEXA_8_FBAR':
C                     
C                        4-------3
C     ETA               /|      /|     F-BAR ISOPARAMETRIC TRI-LINEAR
C     |                / |     / |     8-NODE HEXAHEDRON
C     |               /  |    /  |
C     o--- XI        /   1---/---2 
C    /              /   / o /   /      'o' denotes point where Fo is
C   /              8-------7   /       calculated (element centroid)
C ZETA             |  /    |  /
C                  | /     | /
C                  |/      |/
C                  5-------6
C
C
C REFERENCE: Expression (6.17), Zienkiewicz, O.C., Taylor, R.L., 
C            Zhu, J.Z., The Finite Element Method: Its Basis and 
C            Fundamentals. 7th ed. Butterworth-Heinemann, 2013.

C***********************************************************************
      IF(IBOUND.EQ.0)THEN
C Shape functions and derivatives on element DOMAIN
C -------------------------------------------------
        S=EXISP
        T=ETASP
        U=ZETASP
        ST=S*T
        SU=S*U
        TU=T*U
        STU=ST*U
C Shape functions
        SHAPE(1)=(-U+TU-T+SU-STU+ST-S+R1)*RP125
        SHAPE(2)=(-U+TU-T-SU+STU-ST+S+R1)*RP125
        SHAPE(3)=(-U-TU+T-SU-STU+ST+S+R1)*RP125
        SHAPE(4)=(-U-TU+T+SU+STU-ST-S+R1)*RP125
        SHAPE(5)=( U-TU-T-SU+STU+ST-S+R1)*RP125
        SHAPE(6)=( U-TU-T+SU-STU-ST+S+R1)*RP125
        SHAPE(7)=( U+TU+T+SU+STU+ST+S+R1)*RP125
        SHAPE(8)=( U+TU+T-SU-STU-ST-S+R1)*RP125
C Shape function derivatives
        DERIV(1,1)=( U-TU+T-R1)*RP125
        DERIV(1,2)=(-U+TU-T+R1)*RP125
        DERIV(1,3)=(-U-TU+T+R1)*RP125
        DERIV(1,4)=( U+TU-T-R1)*RP125
        DERIV(1,5)=(-U+TU+T-R1)*RP125
        DERIV(1,6)=( U-TU-T+R1)*RP125
        DERIV(1,7)=( U+TU+T+R1)*RP125
        DERIV(1,8)=(-U-TU-T-R1)*RP125
        DERIV(2,1)=( U-SU+S-R1)*RP125
        DERIV(2,2)=( U+SU-S-R1)*RP125
        DERIV(2,3)=(-U-SU+S+R1)*RP125
        DERIV(2,4)=(-U+SU-S+R1)*RP125
        DERIV(2,5)=(-U+SU+S-R1)*RP125
        DERIV(2,6)=(-U-SU-S-R1)*RP125
        DERIV(2,7)=( U+SU+S+R1)*RP125
        DERIV(2,8)=( U-SU-S+R1)*RP125
        DERIV(3,1)=( T-ST+S-R1)*RP125
        DERIV(3,2)=( T+ST-S-R1)*RP125
        DERIV(3,3)=(-T-ST-S-R1)*RP125
        DERIV(3,4)=(-T+ST+S-R1)*RP125
        DERIV(3,5)=(-T+ST-S+R1)*RP125
        DERIV(3,6)=(-T-ST+S+R1)*RP125
        DERIV(3,7)=( T+ST+S+R1)*RP125
        DERIV(3,8)=( T-ST-S+R1)*RP125
      ELSE
C Shape function and derivatives on element BOUNDARY (2-D)
C (these are the QUAD4 shape functions)
C --------------------------------------------------------
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
C
      ENDIF
C
      RETURN
      END
CDOC END_SUBROUTINE SFH8FB
