CDOC BEGIN_SUBROUTINE GETGCO
CDOC Gets coordinates of a point within an element by interpolation
CDOC
CDOC This routine computes the global cartesian coordinates of a point
CDOC within a finite element by interpolation of its nodal coordinates.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION CARTCO <  Cartesian coordinates of the point of
CDOC C                          interest.
CDOC DOUBLE_PRECISION ELCOD  >  Array of nodal coordinates of the
CDOC C                          element.
CDOC INTEGER          MDIME  >  Dimensioning parameter: Number of rows
CDOC C                          of ELCOD.
CDOC INTEGER          NDIME  >  Number of spatial dimensions.
CDOC INTEGER          NNODE  >  Number of nodes of the element.
CDOC DOUBLE_PRECISION SHAPE  >  Array containing the value of the shape
CDOC C                          function at the point of interest.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto  August  1996
CHST
      SUBROUTINE GETGCO
     1(   CARTCO     ,ELCOD      ,MDIME      ,NDIME      ,NNODE      ,
     2    SHAPE      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION
     1    CARTCO(NDIME)      ,ELCOD(MDIME,NNODE) ,SHAPE(NNODE)
      DATA R0/0.0D0/
C***********************************************************************
C EVALUATES THE GLOBAL CARTESIAN COORDINATES OF A POINT WITHIN AN
C ELEMENT BY INTERPOLATION OF THE ELEMENT NODAL COORDINATES
C
C REFERENCE: Expressions (4.39-40)
C***********************************************************************
      DO 20 IDIME=1,NDIME
        CARTCO(IDIME)=R0
        DO 10 INODE=1,NNODE
          CARTCO(IDIME)=CARTCO(IDIME)+ELCOD(IDIME,INODE)*SHAPE(INODE)
   10   CONTINUE
   20 CONTINUE
C
      RETURN
      END
CDOC END_SUBROUTINE GETGCO
