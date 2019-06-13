CDOC BEGIN_SUBROUTINE GETGMX
CDOC Computes the discrete (full) gradient operator for 2-D and 3-D 
CDOC elements.
CDOC
CDOC This routine assembles the discrete gradient operator, the 
CDOC G-matrix, for isoparametric 3-D and 2-D finite elements: plane 
CDOC strain, plane stress and axisymmetric cases.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION CARTCO >  Cartesian coordinates of the point where
CDOC C                          the G-matrix is to be computed.
CDOC DOUBLE_PRECISION CARTD  >  Array of cartesian derivatives of the
CDOC C                          element shape functions at the point of
CDOC C                          interest.
CDOC DOUBLE_PRECISION GMATX  <  The discrete gradient operator,
CDOC C                          G-matrix.
CDOC INTEGER          MDIME  >  Dimensioning parameter: Number of rows
CDOC C                          of CARTD.
CDOC INTEGER          MGDIM  >  Dimensioning parameter: Number of rows
CDOC C                          of GMATX.
CDOC INTEGER          NAXIS  >  Axis of symmetry flag. Used only for the
CDOC C                          axisymmetric case.
CDOC INTEGER          NNODE  >  Number of nodes of the element.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC DOUBLE_PRECISION SHAPE  >  Array containing the value of the shape
CDOC C                          functions at the point of interest.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, August 1996: Initial coding
CHST D. de Bortoli  , March  2015: 3-D case added
CHST
      SUBROUTINE GETGMX
     1(   CARTCO     ,CARTD      ,GMATX      ,MDIME      ,MGDIM      ,
     2    NAXIS      ,NNODE      ,NTYPE      ,SHAPE      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION
     1    CARTCO(MDIME)      ,CARTD(MDIME,*)     ,GMATX(MGDIM,*)     ,
     2    SHAPE(*)
      LOGICAL IS2D
      DATA R0/0.0D0/
C***********************************************************************
C EVALUATES THE DISCRETE (FULL) GRADIENT OPERATOR 'G' FOR:
C A) 2-D PROBLEMS - PLANE STRESS/STRAIN AND AXISYMMETRIC PROBLEMS;
C                   COMPONENT ORDERING (11,21,12,22,33).
C
C B) 3-D PROBLEMS - COMPONENT ORDERING (11,21,31,12,22,23,31,32,33).
C
C REFERENCE: Expression (4.97)
C***********************************************************************
      IF((NTYPE.EQ.1).OR.(NTYPE.EQ.2).OR.(NTYPE.EQ.3))THEN
        IS2D=.TRUE.
      ELSEIF(NTYPE.EQ.4)THEN
        IS2D=.FALSE.
      ELSE
        CALL ERRPRT('EI0069')
      ENDIF
C 2-D problems
C ============
      IF(IS2D)THEN
C Plane strain/stress
C -------------------
        IY=0
        DO 10 INODE=1,NNODE
          IX=IY+1
          IY=IX+1
          GMATX(1,IX)=CARTD(1,INODE)
          GMATX(1,IY)=R0
          GMATX(2,IX)=R0
          GMATX(2,IY)=CARTD(1,INODE)
          GMATX(3,IX)=CARTD(2,INODE)
          GMATX(3,IY)=R0
          GMATX(4,IX)=R0
          GMATX(4,IY)=CARTD(2,INODE)
   10   CONTINUE
        IF(NTYPE.EQ.3)THEN
C Axisymmetric problem
C --------------------
          IY=0
          DO 20 INODE=1,NNODE
            IX=IY+1
            IY=IX+1
            IF(NAXIS.EQ.1)THEN
C Axisymmetric about Y axis
              GMATX(5,IX)=SHAPE(INODE)/CARTCO(NAXIS)
              GMATX(5,IY)=R0
            ELSE IF(NAXIS.EQ.2)THEN
C Axisymmetric about X axis
              GMATX(5,IX)=R0
              GMATX(5,IY)=SHAPE(INODE)/CARTCO(NAXIS)
            ENDIF
   20     CONTINUE
        ENDIF
C 3-D problems
C ============
      ELSE
        IZ=0
        DO 30 INODE=1,NNODE
          IX=IZ+1
          IY=IX+1
          IZ=IY+1
C
          GMATX(1,IX)=CARTD(1,INODE)
          GMATX(1,IY)=R0
          GMATX(1,IZ)=R0
          GMATX(2,IX)=R0
          GMATX(2,IY)=CARTD(1,INODE)
          GMATX(2,IZ)=R0
          GMATX(3,IX)=R0
          GMATX(3,IY)=R0
          GMATX(3,IZ)=CARTD(1,INODE)
C
          GMATX(4,IX)=CARTD(2,INODE)
          GMATX(4,IY)=R0
          GMATX(4,IZ)=R0
          GMATX(5,IX)=R0
          GMATX(5,IY)=CARTD(2,INODE)
          GMATX(5,IZ)=R0
          GMATX(6,IX)=R0
          GMATX(6,IY)=R0
          GMATX(6,IZ)=CARTD(2,INODE)
C          
          GMATX(7,IX)=CARTD(3,INODE)
          GMATX(7,IY)=R0
          GMATX(7,IZ)=R0
          GMATX(8,IX)=R0
          GMATX(8,IY)=CARTD(3,INODE)
          GMATX(8,IZ)=R0
          GMATX(9,IX)=R0
          GMATX(9,IY)=R0
          GMATX(9,IZ)=CARTD(3,INODE)
   30   CONTINUE
      ENDIF
C
      RETURN
      END
CDOC END_SUBROUTINE GETGMX
