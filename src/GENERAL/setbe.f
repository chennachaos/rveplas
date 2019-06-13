CDOC BEGIN_SUBROUTINE SETBE
CDOC Obtains the left Cauchy-Green strain tensor from the log strain
CDOC
CDOC Given the Eulerian (Lagrangian) logarithmic strain tensor, this
CDOC routine computes the corresponding left (right) Cauchy-Green
CDOC strain tensor.
CDOC Plane strain, plane stress, axisymmetric and 3-D implementations.
CDOC In the present implementation, the out-of-plane component is
CDOC always computed for the plane strain/stress and axisymmetric cases.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION BE     <> Array of engineering logarithmic strain
CDOC C                          components on entry. Returns as the
CDOC C                          array of components of the corresponding
CDOC C                          Cauchy-Green strain tensor.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, August 1996: Initial coding
CHST E.de Souza Neto, June   2003: Out-of-plane always computed
CHST D. de Bortoli  , March  2015: 3-D case added
CHST
      SUBROUTINE SETBE
     1(   BE         ,NTYPE      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL  EXP2X
      LOGICAL   OUTOFP, IS2D
      DIMENSION BE(*)
      DATA  RP5  /0.5D0/
C***********************************************************************
C COMPUTES THE ELASTIC CAUCHY-GREEN TENSOR AS A FUNCTION OF
C THE ELASTIC LOGARITHMIC STRAIN TENSOR:
C
C                      Be   :=  exp[ 2 Ee  ]
C
C REFERENCE: Box 14.3, item (ii)
C***********************************************************************
      IF(NTYPE.EQ.2.OR.NTYPE.EQ.3)THEN
        OUTOFP=.TRUE.
        IS2D=.TRUE.
      ELSEIF(NTYPE.EQ.1)THEN
        OUTOFP=.TRUE.
        IS2D=.TRUE.
      ELSEIF(NTYPE.EQ.4)THEN
        OUTOFP=.FALSE.
        IS2D=.FALSE.
      ELSE
        CALL ERRPRT('EI0024')
      ENDIF
C Convert engineering elastic strain components into physical components
C and use isotropic tensor function to compute elastic Cauchy-Green
C strain tensor...
      IF(IS2D)THEN
C ...plane strain/stress and axisymmetric cases
        BE(3)=RP5*BE(3)
        CALL ISO2
     1(   EXP2X      ,OUTOFP     ,BE         ,BE         )
C ...three-dimensional case
      ELSE
        BE(4:6)=RP5*BE(4:6)
        CALL ISO3
     1(   EXP2X      ,BE         ,BE         )
      ENDIF
C
      RETURN
      END
CDOC END_SUBROUTINE SETBE
