CDOC BEGIN_SUBROUTINE CHKNDB
CDOC Checks whether a set of nodes matches one of the element boundaries
CDOC
CDOC This routine checks whether a given set of local node numbers match
CDOC one of the boundaries of an element. If the given set matches one
CDOC boundary, NODCHK returns with the appropriate node number ordering
CDOC for numerical integration on that boundary.
CDOC
CDOC BEGIN_PARAMETERS
CDOC LOGICAL          FOUND  <  Set to .TRUE. if the given set of local
CDOC C                          node numbers matches one of the
CDOC C                          boundaries of the element. Set to
CDOC C                          .FALSE. otherwise.
CDOC INTEGER          NNODE  >  Total number of nodes of the element.
CDOC INTEGER          NEDGEL >  Total number of edges (facets in 3-D) of
CDOC C                          the element.
CDOC INTEGER          NODCHK <> On entry, defines the set of local node
CDOC C                          numbers to be matched by a boundary: 
CDOC C                          If local node INODE is part of the set,
CDOC C                          then NODCHK(INODE)=1.
CDOC C                          NODCHK(INODE)=0 otherwise.
CDOC C                          This array returns with the sequence of
CDOC C                          local node numbers on the matching
CDOC C                          boundary.
CDOC INTEGER          NORDEB >  Array containing the sequence of local
CDOC C                          node numbers in each edge (facets in
CDOC C                          3-D) of the element.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, August 1996: Initial coding
CHST
      SUBROUTINE CHKNDB
     1(   FOUND      ,NNODE      ,NEDGEL     ,NODCHK     ,NORDEB     )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL FOUND
      DIMENSION
     1    NODCHK(NNODE)      ,NORDEB(NNODE,NEDGEL)
C***********************************************************************
C CHECKS WHETHER A GIVEN SET OF LOCAL ELEMENT NODE NUMBERS CORRESPOND TO
C ONE OF THE ELEMENT BOUNDARIES (EDGES IN 2-D AND FACETS IN 3-D). IF IT
C DOES, RETURNS (STORED IN NODCHK) THE LOCAL NODE NUMBERS ORDERED FOR
C NUMERICAL INTEGRATION ON BOUNDARY.
C***********************************************************************
      FOUND=.FALSE.
C Searches for the element boundary whose nodes coincide with the given
C set
      DO 20 IEDGEL=1,NEDGEL
        DO 10 INODE=1,NNODE
          IF((NODCHK(INODE).NE.0.AND.NORDEB(INODE,IEDGEL).EQ.0).OR.
     1       (NODCHK(INODE).EQ.0.AND.NORDEB(INODE,IEDGEL).NE.0))GOTO 20
   10   CONTINUE
        FOUND=.TRUE.
        GOTO 30
   20 CONTINUE
C
   30 CONTINUE
      IF(FOUND)THEN
C If the given node set corresponds indeed to one of the boundaries of
C the element, stores the node numbers in NODCHK ordered for numerical
C integration on the boundary.
        DO 40 INODE=1,NNODE
          INODEG=NORDEB(INODE,IEDGEL)
          IF(INODEG.NE.0)NODCHK(INODEG)=INODE
   40   CONTINUE 
      ENDIF 
C
      RETURN
      END
CDOC END_SUBROUTINE CHKNDB
