CDOC BEGIN_SUBROUTINE DEFGRA
CDOC Deformation gradient for 2-D and 3-D isoparametric finite elements
CDOC
CDOC Given the element nodal displacements and the discrete gradient
CDOC operator, G-matrix, at a point, this routine computes the
CDOC corresponding deformation gradient at that point. This routine
CDOC contains the plane strain, plane stress, axisymmetric and three-
CDOC dimensional
CDOC implementations.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION ELDISP >  Array of nodal displacements of the
CDOC C                          finite element.
CDOC DOUBLE_PRECISION F      <  Deformation gradient.
CDOC DOUBLE_PRECISION GMATX  >  Discrete (full) gradient operator,
CDOC C                          G-matrix, at the point of interest.
CDOC INTEGER          MDOFN  >  Dimensioning parameter: number of
CDOC C                          rows of array ELDISP.
CDOC INTEGER          MGDIM  >  Dimensioning parameter: number of
CDOC C                          rows of array GMATX.
CDOC INTEGER          NDOFN  >  Number of degrees of freedom per node.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC INTEGER          NNODE  >  Number of nodes of the element.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, August 1996: Initial coding
CHST D. de Bortoli  , March  2015: 3-D case added
CHST
      SUBROUTINE DEFGRA
     1(   ELDISP     ,F          ,GMATX      ,MDOFN      ,MGDIM      ,
     2    NDOFN      ,NTYPE      ,NNODE      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER
     1(   MCOMP=9    )
      DIMENSION
     1    ELDISP(MDOFN,*)    ,F(3,3)             ,GMATX(MGDIM,*)
      LOGICAL OUTOFP
      DIMENSION
     1    FVEC(MCOMP)
      DATA  R0   ,R1    / 0.0D0,1.0D0 /
C***********************************************************************
C COMPUTES THE DEFORMATION GRADIENT TENSOR ASSOCIATED WITH THE ELEMENT
C DISPLACEMENT 'ELDISP'
C***********************************************************************
C Set total number of deformation gradient components
      IF(NTYPE.EQ.1.OR.NTYPE.EQ.2)THEN
        NCOMP=4
        NDIM=2
        OUTOFP=.FALSE.
      ELSEIF(NTYPE.EQ.3)THEN
        NCOMP=5
        NDIM=2
        OUTOFP=.TRUE.
      ELSEIF(NTYPE.EQ.4)THEN
        NCOMP=9
        NDIM=3
        OUTOFP=.FALSE.
      ELSE
        CALL ERRPRT('EI0021')
      ENDIF
C Evaluate the deformation gradient stored in vector form
      CALL RVZERO(FVEC,NCOMP)
      DO 30 ICOMP=1,NCOMP
        IEVAB=0
        DO 20 INODE=1,NNODE
          DO 10 IDOFN=1,NDOFN
            IEVAB=IEVAB+1
            FVEC(ICOMP)=FVEC(ICOMP)+
     1                  GMATX(ICOMP,IEVAB)*ELDISP(IDOFN,INODE)
   10     CONTINUE
   20   CONTINUE
   30 CONTINUE
C Store the deformation gradient in matrix form
      CALL ATASYM(FVEC, F, NDIM, OUTOFP)
C Add identity to F
      F(1,1)=F(1,1)+R1
      F(2,2)=F(2,2)+R1
      IF(NTYPE.NE.1)THEN
        F(3,3)=F(3,3)+R1
      ENDIF
C
      RETURN
      END
CDOC END_SUBROUTINE DEFGRA
