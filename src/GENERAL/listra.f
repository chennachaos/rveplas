CDOC BEGIN_SUBROUTINE LISTRA
CDOC Computes the infinitesimal strain components in 2-D and 3-D
CDOC
CDOC Given the nodal displacements of the element and the B-matrix
CDOC (discrete symmetric gradient operator) at a point in the element
CDOC domain, this routine computes the corresponding (engineering)
CDOC infinitesimal strain components at that point by performing
CDOC the standard operation: e = B u, where e is the array of
CDOC engineering strain components and u is the array of nodal
CDOC displacements. This routine contains the plane strain/stress,
CDOC axisymmetric and three-dimensional implementations.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION BMATX  >  The discrete symmetric gradient
CDOC C                          operator, B-matrix.
CDOC DOUBLE_PRECISION ELDISP >  Array containing the element nodal
CDOC C                          displacements.
CDOC INTEGER          MDOFN  >  Dimensioning parameter: Number of rows
CDOC C                          of ELDISP.
CDOC INTEGER          MBDIM  >  Dimensioning parameter: Number of rows
CDOC C                          of BMATX.
CDOC INTEGER          NDOFN  >  Number of degrees of freedom per node.
CDOC INTEGER          NNODE  >  Number of nodes of the element.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC DOUBLE_PRECISION STRAN  <  Array of engineering infinitesimal
CDOC C                          strain components.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto,   June 1996: Initial coding
CHST D. de Bortoli  , March  2015: 3-D case added
CHST
      SUBROUTINE LISTRA
     1(   BMATX      ,ELDISP     ,MDOFN      ,MBDIM      ,NDOFN      ,
     2    NNODE      ,NTYPE      ,STRAN      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION
     1    BMATX(MBDIM,*)     ,ELDISP(MDOFN,*)    ,STRAN(*)
C***********************************************************************
C COMPUTES THE SYMMETRIC GRADIENT (LINEAR STRAIN MEASURE) ASSOCIATED
C WITH THE ELEMENT DISPLACEMENT 'ELDISP' IN 3-D AND 2-D: PLANE STRAIN, 
C PLANE STRESS AND AXISYMMETRIC PROBLEMS
C
C REFERENCE: Expression (4.53)
C***********************************************************************
      IF(NTYPE.EQ.1)THEN
        NSTRE=3
        NBDIM=3
      ELSEIF(NTYPE.EQ.2)THEN
        NSTRE=4
        NBDIM=3
      ELSEIF(NTYPE.EQ.3)THEN
        NSTRE=4
        NBDIM=4
      ELSEIF(NTYPE.EQ.4)THEN
        NSTRE=6
        NBDIM=6
      ELSE
        CALL ERRPRT('EI0023')
      ENDIF
C
      CALL RVZERO(STRAN,NSTRE)
      DO 30 ISTRE=1,NBDIM
        IEVAB=0
        DO 20 INODE=1,NNODE
          DO 10 IDOFN=1,NDOFN
            IEVAB=IEVAB+1
            STRAN(ISTRE)=STRAN(ISTRE)+
     1                   BMATX(ISTRE,IEVAB)*ELDISP(IDOFN,INODE)
   10     CONTINUE
   20   CONTINUE
   30 CONTINUE
      RETURN
      END
CDOC END_SUBROUTINE LISTRA
