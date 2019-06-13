CDOC BEGIN_SUBROUTINE GETBMX
CDOC Computes the discrete symmetric gradient operator for 2-D and 3-D
CDOC elements
CDOC 
CDOC This routine assembles the discrete symmetric gradient operator
CDOC (strain-displacement matrix in small strains), the B-matrix, for
CDOC isoparametric 3-D and 2-D finite elements: plane strain, plane 
CDOC stress and axisymmetric cases.
CDOC 
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION BMATX  <  The discrete symmetric gradient
CDOC C                          operator, B-matrix.
CDOC DOUBLE_PRECISION CARTCO >  Cartesian coordinates of the point where
CDOC C                          the B-matrix is to be computed.
CDOC DOUBLE_PRECISION CARTD  >  Array of cartesian derivatives of the
CDOC C                          element shape functions at the point of
CDOC C                          interest.
CDOC INTEGER          NDIME  >  Dimensioning parameter: Number of rows
CDOC C                          of CARTCO and CARTD.
CDOC INTEGER          MBDIM  >  Dimensioning parameter: Number of rows
CDOC C                          of BMATX.
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
      SUBROUTINE GETBMX
     1(   BMATX      ,CARTCO     ,CARTD      ,NDIME      ,MBDIM      ,  
     2    NAXIS      ,NNODE      ,NTYPE      ,SHAPE      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION
     1    BMATX(MBDIM,*)     ,CARTCO(NDIME)      ,CARTD(NDIME,*)     ,
     2    SHAPE(*)
      LOGICAL IS2D
      DATA R0/0.0D0/
C***********************************************************************
C EVALUATES THE DISCRETE SYMMETRIC GRADIENT OPERATOR 'B' (SMALL
C STRAIN-DISPLACEMENT MATRIX) FOR:
C A) 2-D PROBLEMS - PLANE STRESS/STRAIN AND AXISYMMETRIC PROBLEMS;
C                   COMPONENT ORDERING (11,22,12,33).
C
C B) 3-D PROBLEMS - COMPONENT ORDERING (11,22,33,12,23,13).
C
C REFERENCE: Expression (4.30)
C***********************************************************************
      IF((NTYPE.EQ.1).OR.(NTYPE.EQ.2).OR.(NTYPE.EQ.3))THEN
        IS2D=.TRUE.
      ELSEIF(NTYPE.EQ.4)THEN
        IS2D=.FALSE.
      ELSE
        CALL ERRPRT('EI0070')
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
          BMATX(1,IX)=CARTD(1,INODE)
          BMATX(1,IY)=R0
          BMATX(2,IX)=R0
          BMATX(2,IY)=CARTD(2,INODE)
          BMATX(3,IX)=CARTD(2,INODE)
          BMATX(3,IY)=CARTD(1,INODE)
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
              BMATX(4,IX)=SHAPE(INODE)/CARTCO(NAXIS)
              BMATX(4,IY)=R0
            ELSE IF(NAXIS.EQ.2)THEN
C Axisymmetric about X axis
              BMATX(4,IX)=R0
              BMATX(4,IY)=SHAPE(INODE)/CARTCO(NAXIS)
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
          BMATX(1,IX)=CARTD(1,INODE)
          BMATX(1,IY)=R0
          BMATX(1,IZ)=R0
          BMATX(2,IX)=R0
          BMATX(2,IY)=CARTD(2,INODE)
          BMATX(2,IZ)=R0
          BMATX(3,IX)=R0
          BMATX(3,IY)=R0
          BMATX(3,IZ)=CARTD(3,INODE)
C
          BMATX(4,IX)=CARTD(2,INODE)
          BMATX(4,IY)=CARTD(1,INODE)
          BMATX(4,IZ)=R0
          BMATX(5,IX)=R0
          BMATX(5,IY)=CARTD(3,INODE)
          BMATX(5,IZ)=CARTD(2,INODE)
          BMATX(6,IX)=CARTD(3,INODE)
          BMATX(6,IY)=R0
          BMATX(6,IZ)=CARTD(1,INODE)
   30   CONTINUE
      ENDIF
C
      RETURN
      END
CDOC END_SUBROUTINE GETBMX
