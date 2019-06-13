CDOC BEGIN_SUBROUTINE GAUS1D
CDOC Set Gauss point positions and weights for 1-D Gauss quadratures
CDOC
CDOC Given the required number of integration points,
CDOC this routine sets the sampling point positions and the
CDOC corresponding weights for 1-D Gauss quadratures for numerical
CDOC integration over the domain [-1,1].
CDOC
CDOC BEGIN_PARAMETERS
CDOC INTEGER          NGAUS  >  Number of Gauss points.
CDOC DOUBLE_PRECISION POSGP  <  Array containing the Gauss point
CDOC C                          positions in the standard domain.
CDOC DOUBLE_PRECISION WEIGP  <  Array containing the Gauss point
CDOC C                          weights.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, July 1996: Initial coding
CHST
      SUBROUTINE GAUS1D
     1(   NGAUS      ,POSGP      ,WEIGP      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION
     1    POSGP(*)           ,WEIGP(*)
C***********************************************************************
C SET SAMPLING POINTS POSITIONS AND WEIGHTS FOR GAUSSIAN NUMERICAL
C INTEGRATION RULES IN 1-D (INTEGRATION OVER LIMITS [-1,1]).
C
C REFERENCE: Expression (4.36)
C            OC Zienkiewicz & RL Taylor. The finite element method,
C            Volume 1: The basis. 5th Edn. Butterworth Heinemann, 2000.
C            J Fish & T Belytschko. A first course in finite element
C            analysis. Wiley, Chichester, 2007.
C***********************************************************************
      IF(NGAUS.EQ.1)THEN
        POSGP(1)=0.0D0
        WEIGP(1)=2.0D0
      ELSEIF(NGAUS.EQ.2)THEN
        POSGP(1)=-0.577350269189626D0
        WEIGP(1)=1.0D0
        POSGP(2)=+0.577350269189626D0
        WEIGP(2)=1.0D0
      ELSEIF(NGAUS.EQ.3)THEN
        POSGP(1)=-0.774596669241483D0
        WEIGP(1)=0.555555555555556D0
        POSGP(2)=+0.0D0
        WEIGP(2)=0.888888888888889D0
        POSGP(3)=+0.774596669241483D0
        WEIGP(3)=0.555555555555556D0
      ELSEIF(NGAUS.EQ.4)THEN
        POSGP(1)=-0.861136311594053D0
        WEIGP(1)=0.347854845137454D0
        POSGP(2)=-0.339981043584856D0
        WEIGP(2)=0.652145154862546D0
        POSGP(3)=+0.339981043584856D0
        WEIGP(3)=0.652145154862546D0
        POSGP(4)=+0.861136311594053D0
        WEIGP(4)=0.347854845137454D0
      ELSE
        CALL ERRPRT('EI0004')
      ENDIF
C
      RETURN
      END
CDOC END_SUBROUTINE GAUS1D
