CDOC BEGIN_SUBROUTINE EXH8
CDOC Sets Gauss point-node extrapolation matrix for element type HEXA_8
CDOC
CDOC This routine sets the coefficients matrix for extrapolation of
CDOC fields from Gauss point values to nodal values for element type
CDOC HEXA_8: Standard isoparametric 8-noded tri-linear hexahedron.
CDOC
CDOC BEGIN_PARAMETERS
CDOC INTEGER          NGAUSP >  Number of Gauss points.
CDOC DOUBLE_PRECISION EXMATX <  Extrapolation matrix.
CDOC END_PARAMETERS
CHST
CHST D. de Bortoli, March 2015: Initial coding, based on EXQ4
CHST
      SUBROUTINE EXH8
     1(   NGAUSP     ,EXMATX     )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NNODE=8)
      DIMENSION EXMATX(NNODE,NGAUSP)
      DATA R1    /
     1     1.0D0 /
      DATA
     1    A8                 ,B8                   ,C8                 ,
     2    D8                 /
     3    2.54903810567666D0 ,-0.683012701892219D0 ,0.183012701892219D0,
     4  -0.0490381056766580D0/
C***********************************************************************
C SETS COEFFICIENTS MATRIX (EXMATX) FOR EXTRAPOLATION FROM GAUSS POINTS
C TO NODES FOR ELEMENT TYPE 'HEXA_8' (STANDARD 8-NODED TRI-LINEAR
C HEXAHEDRON)
C
C REFERENCE: Section 5.6.1
C            E Hinton & JS Campbel. Local and global Smoothing of
C            discontinuous finite element functions using a least
C            squares method. Int. J. Num. meth. Engng., 8:461-480, 1974.
C            E Hinton & DRJ Owen. An introduction to finite element
C            computations. Pineridge Press, Swansea, 1979.
C***********************************************************************
      IF(NGAUSP.EQ.1)THEN
        EXMATX(1,1)=R1
        EXMATX(2,1)=R1
        EXMATX(3,1)=R1
        EXMATX(4,1)=R1
        EXMATX(5,1)=R1
        EXMATX(6,1)=R1
        EXMATX(7,1)=R1
        EXMATX(8,1)=R1
      ELSEIF(NGAUSP.EQ.8)THEN
        EXMATX(1,1)=A8
        EXMATX(1,2)=B8
        EXMATX(1,3)=B8
        EXMATX(1,4)=C8
        EXMATX(1,5)=B8
        EXMATX(1,6)=C8
        EXMATX(1,7)=C8
        EXMATX(1,8)=D8
C
        EXMATX(2,1)=B8
        EXMATX(2,2)=C8
        EXMATX(2,3)=C8
        EXMATX(2,4)=D8
        EXMATX(2,5)=A8
        EXMATX(2,6)=B8
        EXMATX(2,7)=B8
        EXMATX(2,8)=C8
C
        EXMATX(3,1)=C8
        EXMATX(3,2)=D8
        EXMATX(3,3)=B8
        EXMATX(3,4)=C8
        EXMATX(3,5)=B8
        EXMATX(3,6)=C8
        EXMATX(3,7)=A8
        EXMATX(3,8)=B8
C
        EXMATX(4,1)=B8
        EXMATX(4,2)=C8
        EXMATX(4,3)=A8
        EXMATX(4,4)=B8
        EXMATX(4,5)=C8
        EXMATX(4,6)=D8
        EXMATX(4,7)=B8
        EXMATX(4,8)=C8
C
        EXMATX(5,1)=B8
        EXMATX(5,2)=A8
        EXMATX(5,3)=C8
        EXMATX(5,4)=B8
        EXMATX(5,5)=C8
        EXMATX(5,6)=B8
        EXMATX(5,7)=D8
        EXMATX(5,8)=C8
C
        EXMATX(6,1)=C8
        EXMATX(6,2)=B8
        EXMATX(6,3)=D8
        EXMATX(6,4)=C8
        EXMATX(6,5)=B8
        EXMATX(6,6)=A8
        EXMATX(6,7)=C8
        EXMATX(6,8)=B8
C
        EXMATX(7,1)=D8
        EXMATX(7,2)=C8
        EXMATX(7,3)=C8
        EXMATX(7,4)=B8
        EXMATX(7,5)=C8
        EXMATX(7,6)=B8
        EXMATX(7,7)=B8
        EXMATX(7,8)=A8
C
        EXMATX(8,1)=C8
        EXMATX(8,2)=B8
        EXMATX(8,3)=B8
        EXMATX(8,4)=A8
        EXMATX(8,5)=D8
        EXMATX(8,6)=C8
        EXMATX(8,7)=C8
        EXMATX(8,8)=B8
      ENDIF
C
      RETURN
      END
CDOC END_SUBROUTINE EXH8
