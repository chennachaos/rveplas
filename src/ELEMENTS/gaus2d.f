CDOC BEGIN_SUBROUTINE GAUS2D
CDOC Sets Gaussian quadrature constants for 2-D domains
CDOC
CDOC Given the type of domain and the required number of integration
CDOC points, this routine sets the sampling point positions and the
CDOC corresponding weights following standard Gauss quadrature rules
CDOC for 2-D domains (quadrilateral and triangular domains).
CDOC
CDOC BEGIN_PARAMETERS
CDOC CHARACTER        DOMAIN >  Character string with domain type flag.
CDOC C                          Entry value of DOMAIN can be either
CDOC C                          'QUA' (for quadrilateral domains or
CDOC C                          'TRI' (for triangular domains).
CDOC INTEGER          NGAUS  >  Number of integration points.
CDOC DOUBLE_PRECISION POSGP  <  Array containing the position of the
CDOC C                          integration points.
CDOC DOUBLE_PRECISION WEIGP  <  Array containing the weights of the
CDOC C                          integration points.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, July 1996: Initial coding
CHST
      SUBROUTINE GAUS2D
     1(   DOMAIN     ,NGAUS      ,POSGP      ,WEIGP      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*3 DOMAIN
      DIMENSION
     1    POSGP(2,*)         ,WEIGP(*)
C***********************************************************************
C SET SAMPLING POINTS POSITIONS AND WEIGHTS FOR GAUSSIAN NUMERICAL
C INTEGRATION RULES IN 2-D
C
C REFERENCE: Expression (4.31)
C            OC Zienkiewicz & RL Taylor. The finite element method,
C            Volume 1: The basis. 5th Edn. Butterworth Heinemann, 2000.
C            J Fish & T Belytschko. A first course in finite element
C            analysis. Wiley, Chichester, 2007.
C***********************************************************************
      IF(DOMAIN.EQ.'QUA')THEN
C
C Integration over quadrilateral domain with vertices
C                                         {(1,1),(1,-1),(-1,-1),(-1,1)}
C
        IF(NGAUS.EQ.1)THEN
          POSGP(1,1)=0.0D0
          POSGP(2,1)=0.0D0
          WEIGP(1)=4.0D0
        ELSEIF(NGAUS.EQ.4)THEN
          POSGP(1,1)=-0.577350269189626D0
          POSGP(2,1)=-0.577350269189626D0
          WEIGP(1)=1.0D0
          POSGP(1,2)=-0.577350269189626D0
          POSGP(2,2)=+0.577350269189626D0
          WEIGP(2)=1.0D0
          POSGP(1,3)=+0.577350269189626D0
          POSGP(2,3)=-0.577350269189626D0
          WEIGP(3)=1.0D0
          POSGP(1,4)=+0.577350269189626D0
          POSGP(2,4)=+0.577350269189626D0
          WEIGP(4)=1.0D0
        ELSEIF(NGAUS.EQ.5)THEN
          POSGP(1,1)=-0.774596669241483D0
          POSGP(2,1)=-0.774596669241483D0
          WEIGP(1)=0.555555555555556D0
          POSGP(1,2)=-0.774596669241483D0
          POSGP(2,2)=+0.774596669241483D0
          WEIGP(2)=0.555555555555556D0
          POSGP(1,3)=+0.774596669241483D0
          POSGP(2,3)=-0.774596669241483D0
          WEIGP(3)=0.555555555555556D0
          POSGP(1,4)=+0.774596669241483D0
          POSGP(2,4)=+0.774596669241483D0
          WEIGP(4)=0.555555555555556D0
          POSGP(1,5)=+0.0D0
          POSGP(2,5)=+0.0D0
          WEIGP(5)=1.777777777777778D0
        ELSEIF(NGAUS.EQ.9)THEN
          POSGP(1,1)=-0.774596669241483D0
          POSGP(2,1)=-0.774596669241483D0
          WEIGP(1)=0.308641975308643D0
          POSGP(1,2)=-0.774596669241483D0
          POSGP(2,2)=+0.0D0
          WEIGP(2)=0.493827160493828D0
          POSGP(1,3)=-0.774596669241483D0
          POSGP(2,3)=+0.774596669241483D0
          WEIGP(3)=0.308641975308643D0
          POSGP(1,4)=+0.0D0
          POSGP(2,4)=-0.774596669241483D0
          WEIGP(4)=0.493827160493828D0
          POSGP(1,5)=+0.0D0
          POSGP(2,5)=+0.0D0
          WEIGP(5)=0.790123456790124D0
          POSGP(1,6)=+0.0D0
          POSGP(2,6)=+0.774596669241483D0
          WEIGP(6)=0.493827160493828D0
          POSGP(1,7)=+0.774596669241483D0
          POSGP(2,7)=-0.774596669241483D0
          WEIGP(7)=0.308641975308643D0
          POSGP(1,8)=+0.774596669241483D0
          POSGP(2,8)=+0.0D0
          WEIGP(8)=0.493827160493828D0
          POSGP(1,9)=+0.774596669241483D0
          POSGP(2,9)=+0.774596669241483D0
          WEIGP(9)=0.308641975308643D0
        ELSE
          CALL ERRPRT('EI0001')
        ENDIF
      ELSEIF(DOMAIN.EQ.'TRI')THEN
C
C Integration over triangular domain with vertices {(0,0),(1,0),(0,1)}
C
        IF(NGAUS.EQ.1)THEN
          POSGP(1,1)=0.333333333333333D0
          POSGP(2,1)=0.333333333333333D0
          WEIGP(1)=0.5D0
        ELSEIF(NGAUS.EQ.3)THEN
          POSGP(1,1)=0.166666666666667D0
          POSGP(2,1)=0.166666666666667D0
          WEIGP(1)=0.166666666666667D0
          POSGP(1,2)=0.666666666666667D0
          POSGP(2,2)=0.166666666666667D0
          WEIGP(2)=0.166666666666667D0
          POSGP(1,3)=0.166666666666667D0
          POSGP(2,3)=0.666666666666667D0
          WEIGP(3)=0.166666666666667D0
        ELSEIF(NGAUS.EQ.6)THEN
          POSGP(1,1)=0.816847572980459D0
          POSGP(2,1)=0.091576213509771D0
          WEIGP(1)=0.05497587182766099D0
          POSGP(1,2)=0.091576213509771D0
          POSGP(2,2)=0.816847572980459D0
          WEIGP(2)=0.05497587182766099D0
          POSGP(1,3)=0.091576213509771D0
          POSGP(2,3)=0.091576213509771D0
          WEIGP(3)=0.05497587182766099D0
          POSGP(1,4)=0.108103018168070D0
          POSGP(2,4)=0.445948490915965D0
          WEIGP(4)=0.1116907948390055D0
          POSGP(1,5)=0.445948490915965D0
          POSGP(2,5)=0.108103018168070D0
          WEIGP(5)=0.1116907948390055D0
          POSGP(1,6)=0.445948490915965D0
          POSGP(2,6)=0.445948490915965D0
          WEIGP(6)=0.1116907948390055D0
        ELSEIF(NGAUS.EQ.7)THEN
          POSGP(1,1)=0.333333333333333D0
          POSGP(2,1)=0.333333333333333D0
          WEIGP(1)=0.1125D0
          POSGP(1,2)=0.797426985353087D0
          POSGP(2,2)=0.101286507323456D0
          WEIGP(2)=0.062969590272414D0
          POSGP(1,3)=0.101286507323456D0
          POSGP(2,3)=0.797426985353087D0
          WEIGP(3)=0.062969590272414D0
          POSGP(1,4)=0.101286507323456D0
          POSGP(2,4)=0.101286507323456D0
          WEIGP(4)=0.062969590272414D0
          POSGP(1,5)=0.470142064105115D0
          POSGP(2,5)=0.470142064105115D0
          WEIGP(5)=0.066197076394253D0
          POSGP(1,6)=0.059715871789770D0
          POSGP(2,6)=0.470142064105115D0
          WEIGP(6)=0.066197076394253D0
          POSGP(1,7)=0.470142064105115D0
          POSGP(2,7)=0.059715871789770D0
          WEIGP(7)=0.066197076394253D0
        ELSE
          CALL ERRPRT('EI0002')
        ENDIF
      ELSE
        CALL ERRPRT('EI0003')
      ENDIF
C
      RETURN
      END
CDOC END_SUBROUTINE GAUS2D
