CDOC BEGIN_SUBROUTINE SFH8
CDOC Shape function and derivatives for element type HEXA_8
CDOC
CDOC This routine computes the shape functions and shape function
CDOC derivatives for the element type HEXA_8: Standard isoparametric
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
CHST D. de Bortoli, March 2015: Initial coding
CHST
      SUBROUTINE SFT7E
     1(   DERIV      ,ZETASP      ,ETASP      ,EXISP      ,IBOUND     ,
     2    MDIME      ,SHAPE      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION
     1    DERIV(MDIME,*)     ,SHAPE(*)
      DATA RP25, RP125, R1   /0.25D0,0.125D0,1.0D0/
C***********************************************************************
C COMPUTES SHAPE FUNCTIONS AND SHAPE FUNCTION DERIVATIVES FOR
C ELEMENT 'TETA_7':
C                     
C                        4-------3
C     ETA               /|      /|     STANDARD ISOPARAMETRIC
C     |                / |     / |     TRI-LINEAR 8-NODE HEXAHEDRON
C     |               /  |    /  |
C     o--- XI        /   1---/---2 
C    /              8-------7   /
C   /               |  /    |  /
C ZETA              | /     | /
C                   |/      |/
C                   5-------6
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
        STU=R1-S-T-U
C Shape functions
        SHAPE(1)=STU*(2.0*STU-R1)
        SHAPE(2)=S*(2.0*S-R1)
        SHAPE(3)=T*(2.0*T-R1)
        SHAPE(4)=U*(2.0*U-R1)
        SHAPE(5)=4.0*S*STU
        SHAPE(6)=4.0*S*U
        SHAPE(7)=4.0*T*STU
        SHAPE(8)=4.0*STU*U
        SHAPE(9)=4.0*S*U
        SHAPE(10)=4.0*T*U
C Shape function derivatives
        DERIV(1,1)=-4.0*STU+R1
        DERIV(1,2)=4.0*S-R1
        DERIV(1,3)=0.0
        DERIV(1,4)=0.0
        DERIV(1,5)=4.0*(STU-S)
        DERIV(1,6)=4.0*T
        DERIV(1,7)=-4.0*T
        DERIV(1,8)=-4*U
        DERIV(1,9)=4*U
        DERIV(1,10)=0.0
        
        DERIV(2,1)=-4.0*STU+R1
        DERIV(2,2)=0.0
        DERIV(2,3)=4.0*T-R1
        DERIV(2,4)=0
        DERIV(2,5)=-4.0*S
        DERIV(2,6)=4.0*S
        DERIV(2,7)=4.0*(STU-T) 
        DERIV(2,8)=-4.0*U
        DERIV(2,9)=0
        DERIV(2,10)=4.0*U        
        
        DERIV(3,1)=-4.0*STU+R1
        DERIV(3,2)= 0.0
        DERIV(3,3)=0.0
        DERIV(3,4)=4.0*U-R1
        DERIV(3,5)=-4.0*S
        DERIV(3,6)=0.0
        DERIV(3,7)=-4.0*T
        DERIV(3,8)=4.0*(STU-U)
        DERIV(3,9)=4.0*S
        DERIV(3,10)=4.0*T
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
CDOC END_SUBROUTINE SFT7E
