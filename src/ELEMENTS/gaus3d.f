CDOC BEGIN_SUBROUTINE GAUS3D
CDOC Sets Gaussian quadrature constants for 3-D domains
CDOC
CDOC Given the type of domain and the required number of integration
CDOC points, this routine sets the sampling point positions and the
CDOC corresponding weights following standard Gauss quadrature rules
CDOC for 3-D domains (hexahedral and tetrahedral domains).
CDOC
CDOC BEGIN_PARAMETERS
CDOC CHARACTER        DOMAIN >  Character string with domain type flag.
CDOC C                          Entry value of DOMAIN can be either
CDOC C                          'HEX' (for hexahedron domains or
CDOC C                          'TET' (for tetrahedron domains).
CDOC INTEGER          NGAUS  >  Number of integration points.
CDOC DOUBLE_PRECISION POSGP  <  Array containing the position of the
CDOC C                          integration points.
CDOC DOUBLE_PRECISION WEIGP  <  Array containing the weights of the
CDOC C                          integration points.
CDOC END_PARAMETERS
CHST
CHST D. de Bortoli, March 2015: Initial coding
CHST
      SUBROUTINE GAUS3D
     1(   DOMAIN     ,NGAUS      ,POSGP      ,WEIGP      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*3 DOMAIN
      DIMENSION
     1    POSGP(3,*)         ,WEIGP(*)
C***********************************************************************
C SET SAMPLING POINTS POSITIONS AND WEIGHTS FOR GAUSSIAN NUMERICAL
C INTEGRATION RULES IN 3-D
C
C REFERENCE: Expression (4.31)
C            OC Zienkiewicz & RL Taylor. The finite element method,
C            Volume 1: The basis. 5th Edn. Butterworth Heinemann, 2000.
C            J Fish & T Belytschko. A first course in finite element
C            analysis. Wiley, Chichester, 2007.
C***********************************************************************
      IF(DOMAIN.EQ.'HEX')THEN
C
C Integration over hexahedral domain with vertices
C                         {( 1,1,1),( 1,1,-1),( 1,-1,1),( 1,-1,-1),
C                          (-1,1,1),(-1,1,-1),(-1,-1,1),(-1,-1,-1)}
C
        IF(NGAUS.EQ.1)THEN
          POSGP(1,1)=0.0D0
          POSGP(2,1)=0.0D0
          POSGP(3,1)=0.0D0
          WEIGP(1)=8.0D0
        ELSEIF(NGAUS.EQ.8)THEN
          POSGP(1,1)=-0.577350269189626D0
          POSGP(2,1)=-0.577350269189626D0
          POSGP(3,1)=-0.577350269189626D0
          WEIGP(1)=1.0D0
          POSGP(1,2)=-0.577350269189626D0
          POSGP(2,2)=-0.577350269189626D0
          POSGP(3,2)=+0.577350269189626D0
          WEIGP(2)=1.0D0
          POSGP(1,3)=-0.577350269189626D0
          POSGP(2,3)=+0.577350269189626D0
          POSGP(3,3)=-0.577350269189626D0
          WEIGP(3)=1.0D0
          POSGP(1,4)=-0.577350269189626D0
          POSGP(2,4)=+0.577350269189626D0
          POSGP(3,4)=+0.577350269189626D0
          WEIGP(4)=1.0D0
          POSGP(1,5)=+0.577350269189626D0
          POSGP(2,5)=-0.577350269189626D0
          POSGP(3,5)=-0.577350269189626D0
          WEIGP(5)=1.0D0
          POSGP(1,6)=+0.577350269189626D0
          POSGP(2,6)=-0.577350269189626D0
          POSGP(3,6)=+0.577350269189626D0
          WEIGP(6)=1.0D0
          POSGP(1,7)=+0.577350269189626D0
          POSGP(2,7)=+0.577350269189626D0
          POSGP(3,7)=-0.577350269189626D0
          WEIGP(7)=1.0D0
          POSGP(1,8)=+0.577350269189626D0
          POSGP(2,8)=+0.577350269189626D0
          POSGP(3,8)=+0.577350269189626D0
          WEIGP(8)=1.0D0
        ELSE
          CALL ERRPRT('EI0014')
        ENDIF
      ELSEIF(DOMAIN.EQ.'TET')THEN
C
C Integration over tetrahedral domain with vertices 
C                                   {(0,0,0),(1,0,0),(0,1,0),(0,0,1)}
C
        IF(NGAUS.EQ.1)THEN
          POSGP(1,1)=0.25D0
          POSGP(2,1)=0.25D0
          POSGP(3,1)=0.25D0
          WEIGP(1)=1.0D0
        ELSEIF(NGAUS.EQ.4)THEN
C          POSGP(1,1)=0.58541020D0
C          POSGP(2,1)=0.13819660D0
C          POSGP(3,1)=0.13819660D0
C          WEIGP(1)=0.25D0
C          POSGP(1,2)=0.13819660D0
C          POSGP(2,2)=0.58541020D0
C          POSGP(3,2)=0.13819660D0
C          WEIGP(2)=0.25D0
C          POSGP(1,3)=0.13819660D0
C          POSGP(2,3)=0.13819660D0
C          POSGP(3,3)=0.58541020D0
C          WEIGP(3)=0.25D0
C          POSGP(1,4)=0.13819660D0
C          POSGP(2,4)=0.13819660D0
C          POSGP(3,4)=0.13819660D0
C          WEIGP(4)=0.25D0
          
          POSGP(1,1)=0.13819660D0
          POSGP(2,1)=0.13819660D0
          POSGP(3,1)=0.13819660D0
          WEIGP(1)=0.25D0
          POSGP(1,2)=0.58541020D0
          POSGP(2,2)=0.13819660D0
          POSGP(3,2)=0.13819660D0
          WEIGP(2)=0.25D0
          POSGP(1,3)=0.13819660D0
          POSGP(2,3)=0.58541020D0
          POSGP(3,3)=0.13819660D0
          WEIGP(3)=0.25D0
          POSGP(1,4)=0.13819660D0
          POSGP(2,4)=0.13819660D0
          POSGP(3,4)=0.58541020D0
          WEIGP(4)=0.25D0          
        ELSEIF(NGAUS.EQ.5)THEN
          POSGP(1,1)=0.25D0
          POSGP(2,1)=0.25D0
          POSGP(3,1)=0.25D0
          WEIGP(1)=-0.8D0
          POSGP(1,2)=0.333333333333333D0
          POSGP(2,2)=0.166666666666667D0
          POSGP(3,2)=0.166666666666667D0
          WEIGP(2)=0.45D0
          POSGP(1,3)=0.166666666666667D0
          POSGP(2,3)=0.333333333333333D0
          POSGP(3,3)=0.166666666666667D0
          WEIGP(3)=0.45D0
          POSGP(1,4)=0.166666666666667D0
          POSGP(2,4)=0.166666666666667D0
          POSGP(3,4)=0.333333333333333D0
          WEIGP(4)=0.45D0
          POSGP(1,5)=0.166666666666667D0
          POSGP(2,5)=0.166666666666667D0
          POSGP(3,5)=0.166666666666667D0
          WEIGP(5)=0.45D0
        ELSE
          CALL ERRPRT('EI0015')
        ENDIF
      ELSE
        CALL ERRPRT('EI0077')
      ENDIF
C
      RETURN
      END
CDOC END_SUBROUTINE GAUS3D
