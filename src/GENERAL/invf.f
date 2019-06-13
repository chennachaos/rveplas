CDOC BEGIN_SUBROUTINE INVF
CDOC Inverts the deformation gradient for 2-D and 3-D problems
CDOC
CDOC This routine inverts deformation gradient tensors for plane strain,
CDOC plane stress, axisymmetric and three-dimensional problems.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION F      >  Deformation gradient.
CDOC DOUBLE_PRECISION FINV   <  Inverse of the deformation gradient.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, August 1996: Initial coding
CHST D. de Bortoli  , March  2015: 3-D case added, changed name from
CHST                               INVF2 to INVF
CHST
      SUBROUTINE INVF
     1(   F          ,FINV       ,NTYPE      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION
     1    F(3,3)             ,FINV(3,3)
      LOGICAL IS2D
      DATA
     1   R0   ,R1   /
     2   0.0D0,1.0D0/
C***********************************************************************
C INVERTS DEFORMATION GRADIENT TENSORS FOR PLANE STRESS/STRAIN,
C AXISYMMETRIC AND THREE-DIMENSIONAL PROBLEMS
C***********************************************************************
      IF((NTYPE.EQ.1).OR.(NTYPE.EQ.2).OR.(NTYPE.EQ.3))THEN
        IS2D=.TRUE.
      ELSEIF(NTYPE.EQ.4)THEN
        IS2D=.FALSE.
      ELSE
        CALL ERRPRT('EI0072')
      ENDIF
C 2-D problems
C ------------
      IF(IS2D)THEN
        DETFPL=F(1,1)*F(2,2)-F(1,2)*F(2,1)
        IF(DETFPL.EQ.R0)CALL ERRPRT('EE0001')
        IF(NTYPE.EQ.3.AND.F(3,3).EQ.R0)CALL ERRPRT('EE0001')
C
        DETFIN=R1/DETFPL
        FINV(1,1)=F(2,2)*DETFIN
        FINV(2,2)=F(1,1)*DETFIN
        FINV(1,2)=-F(1,2)*DETFIN
        FINV(2,1)=-F(2,1)*DETFIN
        IF(NTYPE.EQ.2)THEN
          FINV(3,3)=R1
        ELSEIF(NTYPE.EQ.3)THEN
          FINV(3,3)=R1/F(3,3)
        ENDIF
C 3-D problems
C ------------
      ELSE
        CALL INVMT3
     1(   F          ,FINV       ,DETF3D       )
      ENDIF
C
      RETURN
      END
CDOC END_SUBROUTINE INVF
