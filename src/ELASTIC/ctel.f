CDOC BEGIN_SUBROUTINE CTEL
CDOC Computation of the elastic matrix for the linear elastic model.
CDOC
CDOC This routine returns the standard elasticity matrix of the linear
CDOC elastic material model. It contains the plane stress, plane strain
CDOC and axisymmetric implementations.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DMATX  <  Elastic matrix.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC DOUBLE_PRECISION RPROPS >  Array of real material properties.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, September 1996: Initial coding
CHST D. de Bortoli  , March     2015: 3-D case added (NTYPE=4)
CHST
      SUBROUTINE CTEL
     1(   DMATX      ,NTYPE      ,RPROPS     )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(MSTRE=6)
      DIMENSION
     1    DMATX(MSTRE,MSTRE),RPROPS(*)
      DIMENSION
     1    FOID(4,4) ,FOID3D(6,6), SOID(4), SOID3D(6)
C Fourth-order symmetric identity tensor (for plane stress/strain and 
C axisymmetric conditions)
      DATA
     1    FOID(1,1),FOID(1,2),FOID(1,3),FOID(1,4)/
     2    1.0D0    ,0.0D0    ,0.0D0    ,0.0D0    /
     3    FOID(2,1),FOID(2,2),FOID(2,3),FOID(2,4)/
     4    0.0D0    ,1.0D0    ,0.0D0    ,0.0D0    /
     5    FOID(3,1),FOID(3,2),FOID(3,3),FOID(3,4)/
     6    0.0D0    ,0.0D0    ,0.5D0    ,0.0D0    /
     7    FOID(4,1),FOID(4,2),FOID(4,3),FOID(4,4)/
     8    0.0D0    ,0.0D0    ,0.0D0    ,1.0D0    /
C Fourth-order symmetric identity tensor for the three-dimensional case
      DATA
     1    FOID3D(1,1),FOID3D(1,2),FOID3D(1,3),FOID3D(1,4),FOID3D(1,5)/
     2    1.0D0      ,0.0D0      ,0.0D0      ,0.0D0      ,0.0D0      /
     3    FOID3D(1,6)/
     4    0.0D0      /
     5    FOID3D(2,1),FOID3D(2,2),FOID3D(2,3),FOID3D(2,4),FOID3D(2,5)/
     6    0.0D0      ,1.0D0      ,0.0D0      ,0.0D0      ,0.0D0      /
     7    FOID3D(2,6)/
     8    0.0D0      /
     9    FOID3D(3,1),FOID3D(3,2),FOID3D(3,3),FOID3D(3,4),FOID3D(3,5)/
     O    0.0D0      ,0.0D0      ,1.0D0      ,0.0D0      ,0.0D0      /
     1    FOID3D(3,6)/
     2    0.0D0      /
     3    FOID3D(4,1),FOID3D(4,2),FOID3D(4,3),FOID3D(4,4),FOID3D(4,5)/
     4    0.0D0      ,0.0D0      ,0.0D0      ,0.5D0      ,0.0D0      /
     5    FOID3D(4,6)/
     6    0.0D0      /
     7    FOID3D(5,1),FOID3D(5,2),FOID3D(5,3),FOID3D(5,4),FOID3D(5,5)/
     8    0.0D0      ,0.0D0      ,0.0D0      ,0.0D0      ,0.5D0      /
     9    FOID3D(5,6)/
     O    0.0D0      /
     1    FOID3D(6,1),FOID3D(6,2),FOID3D(6,3),FOID3D(6,4),FOID3D(6,5)/
     2    0.0D0      ,0.0D0      ,0.0D0      ,0.0D0      ,0.0D0      /
     3    FOID3D(6,6)/
     4    0.5D0      /
C Second-order identity tensor (for plane stress/strain and axisymmetric
C cases) (array representation)
      DATA
     1    SOID(1)  ,SOID(2)  ,SOID(3)  ,SOID(4)  /
     2    1.0D0    ,1.0D0    ,0.0D0    ,1.0D0    /
C Second-order identity tensor (for the three-dimensional case)
C (array representation)
      DATA
     1    SOID3D(1),SOID3D(2),SOID3D(3),SOID3D(4),SOID3D(5),SOID3D(6)/
     2    1.0D0    ,1.0D0    ,1.0D0    ,0.0D0    ,0.0D0    ,0.0D0    /
      DATA
     1    R1   ,R2   ,R3   ,R4   /
     2    1.0D0,2.0D0,3.0D0,4.0D0/
C***********************************************************************
C COMPUTATION OF THE TANGENT MODULUS (ELASTICITY MATRIX) FOR THE LINEAR
C ELASTIC MATERIAL MODEL
C
C REFERENCE: Expression (4.44)
C***********************************************************************
C
C Set shear and bulk modulus
C --------------------------
C
      GMODU=RPROPS(2)
      BULK=RPROPS(3)
C
      R1D3=R1/R3
      R2G=R2*GMODU
      FACTOR=BULK-R2G*R1D3
C
C Compute elasticity matrix
C -------------------------
C
      IF(NTYPE.EQ.1)THEN
C plane stress
        NSTRE=3
        R4GD3=R4*GMODU*R1D3
        FACTOR=FACTOR*(R2G/(BULK+R4GD3))
      ELSEIF(NTYPE.EQ.2)THEN
C plane strain
        NSTRE=3
      ELSEIF(NTYPE.EQ.3)THEN
C axisymmetric
        NSTRE=4
C three-dimensional
      ELSEIF(NTYPE.EQ.4)THEN
        NSTRE=6
      ELSE
C stops program if other stress state
        CALL ERRPRT('EI0019')
      ENDIF
C
C Assemble matrix...
C
C ... plane stress/strain and axisymmetric cases
      IF((NTYPE.EQ.1.).OR.(NTYPE.EQ.2).OR.(NTYPE.EQ.3))THEN
        DO 20 I=1,NSTRE
          DO 10 J=I,NSTRE
            DMATX(I,J)=R2G*FOID(I,J)+FACTOR*SOID(I)*SOID(J)
   10     CONTINUE
   20   CONTINUE
C ... three-dimensional case
      ELSEIF(NTYPE.EQ.4)THEN
        DO 40 I=1,NSTRE
          DO 30 J=I,NSTRE
            DMATX(I,J)=R2G*FOID3D(I,J)+FACTOR*SOID3D(I)*SOID3D(J)
   30     CONTINUE
   40   CONTINUE
      ENDIF
C lower triangle (all cases)
      DO 60 J=1,NSTRE-1
        DO 50 I=J+1,NSTRE
          DMATX(I,J)=DMATX(J,I)
   50   CONTINUE
   60 CONTINUE
C
      RETURN
      END
CDOC END_SUBROUTINE CTEL
