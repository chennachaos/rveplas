CDOC BEGIN_SUBROUTINE EXTNOD
CDOC Extrapolates Gauss point values of a given field to nodes
CDOC
CDOC This routine extrapolates the Gauss point values of a given field
CDOC to the nodes of a finite element. It simply multiplies the
CDOC coefficients matrix for extrapolation EXMATX by the matrix VARGP
CDOC which contains the value of each component of the given field in
CDOC all integration points of the element.
CDOC The results stored in the argument VARNOD are the extrapolated
CDOC (or locally smoothed) values of the components of the given field
CDOC at the nodes of the element.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION EXMATX >  Extrapolation matrix. Dimension
CDOC C                          NNODE x NGAUSP.
CDOC DOUBLE_PRECISION VARGP  >  Matrix containing, in each column,
CDOC C                          the Gauss point value of all
CDOC C                          components of the given field. 
CDOC C                          Dimension NVAR x NGAUSP.
CDOC DOUBLE_PRECISION VARNOD <  Matrix containing, in each column,
CDOC C                          the (extrapolated) nodal value of all
CDOC C                          components of the given field. 
CDOC C                          Dimension NVAR x NNODE.
CDOC INTEGER          NVAR   >  Number of components of the given field.
CDOC INTEGER          NGAUSP >  Number of Gauss points.
CDOC INTEGER          NNODE  >  Number of nodes of the element.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, August 1996: Initial coding
CHST
      SUBROUTINE EXTNOD
     1(   EXMATX     ,VARGP      ,VARNOD     ,NVAR       ,NGAUSP     ,
     2    NNODE      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE '../ELEMENTS.INC'
      DIMENSION
     1    EXMATX(NNODE,NGAUSP),VARGP(NVAR,NGAUSP),VARNOD(NVAR,NNODE)
C***********************************************************************
C EXTRAPOLATES GAUSS POINT VALUES OF A GIVEN FIELD TO NODES
C
C REFERENCE: Section 5.6.1
C            E Hinton & JS Campbel. Local and global Smoothing of
C            discontinuous finite element functions using a least
C            squares method. Int. J. Num. meth. Engng., 8:461-480, 1974.
C            E Hinton & DRJ Owen. An introduction to finite element
C            computations. Pineridge Press, Swansea, 1979.
C***********************************************************************
      CALL RVZERO(VARNOD,NVAR*NNODE)
      DO 30 IVAR=1,NVAR
        DO 20 INODE=1,NNODE
          DO 10 IGAUSP=1,NGAUSP
            VARNOD(IVAR,INODE)=VARNOD(IVAR,INODE)+
     1                         EXMATX(INODE,IGAUSP)*VARGP(IVAR,IGAUSP)
   10     CONTINUE
   20   CONTINUE
   30 CONTINUE
C
      RETURN
      END
CDOC END_SUBROUTINE EXTNOD
