CDOC BEGIN_SUBROUTINE LOGSTR
CDOC Logarithmic strain computation
CDOC
CDOC Given the left (right) Cauchy-Green strain tensor, this routine
CDOC computes the corresponding Eulerian (Lagrangian) logarithmic strain
CDOC tensor (engineering components).
CDOC Plane strain, plane stress, axisymmetric and 3-D implementations.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION B      >  Array of components of the Cauchy-Green
CDOC C                          strain tensor.
CDOC DOUBLE_PRECISION E      <  Array of (engineering) components of the
CDOC C                          logarithmic strain tensor.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, August 1996: Initial coding (as EETRIA)
CHST E.de Souza Neto, May    1998: Routine and some variables renamed
CHST D. de Bortoli  , March  2015: 3-D case added
CHST
      SUBROUTINE LOGSTR
     1(   B          ,E          ,NTYPE      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL  DLGD2
      LOGICAL   OUTOFP, IS2D
      DIMENSION
     1    B(*)               ,E(*)
      DATA  R2   /2.0D0/
C***********************************************************************
C COMPUTES THE LOGARITHMIC STRAIN TENSOR:
C
C                       E :=  1/2 ln[ B ]
C
C REFERENCE: Box 14.3, item (ii)
C***********************************************************************
      IF(NTYPE.EQ.2.OR.NTYPE.EQ.3)THEN
        OUTOFP=.TRUE.
        IS2D=.TRUE.
      ELSEIF(NTYPE.EQ.1)THEN
        OUTOFP=.FALSE.
        IS2D=.TRUE.
      ELSEIF(NTYPE.EQ.4)THEN
        OUTOFP=.FALSE.
        IS2D=.FALSE.
      ELSE
        CALL ERRPRT('EI0022')
      ENDIF
C Use isotropic tensor function to compute the logarithmic (physical)
C strain components, then convert physical components into engineering
C strain components...
      IF(IS2D)THEN
C ...plane strain/stress and axisymmetric cases
        CALL ISO2
     1(   DLGD2      ,OUTOFP     ,B          ,E          )
C
        E(3)=R2*E(3)
      ELSE
C ...three-dimensional case
        CALL ISO3
     1(   DLGD2      ,B          ,E          )
C
        E(4:6)=R2*E(4:6)
      ENDIF
C
      RETURN
      END
CDOC END_SUBROUTINE LOGSTR
