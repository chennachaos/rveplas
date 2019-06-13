CDOC BEGIN_SUBROUTINE SFT6
CDOC Shape function and derivatives for element type TRI6
CDOC
CDOC This routine computes the shape functions and shape function
CDOC derivatives for the element type TRI6: Standard isoparametric
CDOC 6-noded triangle for plane strain, plane stress and
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
CHST
CHST
      SUBROUTINE SFT6
     1(   DERIV      ,ETASP      ,EXISP      ,IBOUND     ,MDIME      ,
     2    SHAPE      )
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(OUT), DIMENSION(MDIME,*) :: DERIV
      DOUBLE PRECISION, INTENT(OUT), DIMENSION(*) :: SHAPE
      DOUBLE PRECISION, INTENT(IN) :: ETASP, EXISP
      INTEGER, INTENT(IN) :: IBOUND, MDIME
C
      DOUBLE PRECISION :: L1, L2, L3, L14, L24, L34, S, SS, S2
C
      DOUBLE PRECISION, PARAMETER :: R0=0.0D0, RP5=0.5D0, R1=1.0D0,
     1                               R2=2.0D0,  R4=4.0D0
C***********************************************************************
C COMPUTES SHAPE FUNCTIONS AND SHAPE FUNCTION DERIVATIVES FOR
C ELEMENT 'TRI_6':
C                        5 o
C                          |\
C                          | \       
C                          |  \     
C                        6 o   o 4       STANDARD QUADRATIC
C                          |    \        6 NODE TRIANGLE
C                          |     \  
C                          o--o---o
C                         1   2   3
C
C REFERENCE:
C Zienkiewicz, O.C., Taylor, R.L., Zhu, J.Z. The Finite Element Method: 
C Its Basis and Fundamentals. 7th ed. Butterworth-Heinemann, 2013.
C
C Section 6.2.1.3
C
C Note: the node ordering here is different from the reference above
C***********************************************************************
      IF(IBOUND.EQ.0)THEN
C Shape functions and derivatives on element DOMAIN
C -------------------------------------------------
C Area coordinates
        L1=R1-EXISP-ETASP
        L2=EXISP
        L3=ETASP
C
        L14=R4*L1
        L24=R4*L2
        L34=R4*L3
C Shape functions
        SHAPE(1)=L1*(R2*L1 - R1)
        SHAPE(2)=L14*L2
        SHAPE(3)=L2*(R2*L2 - R1)
        SHAPE(4)=L24*L3
        SHAPE(5)=L3*(R2*L3 - R1)
        SHAPE(6)=L34*L1
C Shape function derivatives
        DERIV(1,1)=R1-L14
        DERIV(1,2)=L14-L24
        DERIV(1,3)=L24-R1
        DERIV(1,4)=L34
        DERIV(1,5)=R0
        DERIV(1,6)=-L34
        DERIV(2,1)=R1-L14
        DERIV(2,2)=-L24
        DERIV(2,3)=R0
        DERIV(2,4)=L24
        DERIV(2,5)=L34-R1
        DERIV(2,6)=L14-L34
      ELSE
C Shape function and derivatives on element BOUNDARY (1-D)
C --------------------------------------------------------
        S=EXISP
        SS=S*S
        S2=S*R2
C Shape functions
        SHAPE(1)=(-S+SS)*RP5
        SHAPE(2)=R1-SS
        SHAPE(3)=(S+SS)*RP5
C Shape functions derivatives
        DERIV(1,1)=(-R1+S2)*RP5
        DERIV(1,2)=-S2
        DERIV(1,3)=(R1+S2)*RP5
C
      ENDIF
C
      RETURN
      END
CDOC END_SUBROUTINE SFT6