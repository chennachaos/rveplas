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
      SUBROUTINE EXT7
     1(   NGAUSP     ,EXMATX     )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NNODE=10)
      DIMENSION EXMATX(NNODE,NGAUSP)
      DATA R1    /
     1     1.0D0 /
      DATA
     1    w1                 ,w2                  ,w3                 ,
     2    w4                /
     3   0.2500000000000000 ,0.2500000000000000,0.2500000000000000,
     4  0.2500000000000000/
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
        EXMATX(9,1)=R1
        EXMATX(10,1)=R1
      ELSEIF(NGAUSP.EQ.4)THEN
        EXMATX(1,1)=1.5
        EXMATX(1,2)=-0.5
        EXMATX(1,3)=-0.5
        EXMATX(1,4)=-0.5
        EXMATX(1,5)=0.5
        EXMATX(1,6)=-0.5
        EXMATX(1,7)=0.5
        EXMATX(1,8)=0.5
        EXMATX(1,9)=-0.5
        EXMATX(1,10)=-0.5
C
        EXMATX(2,1)=-0.16666667
        EXMATX(2,2)=1.83333333
        EXMATX(2,3)=-0.16666667
        EXMATX(2,4)=-0.16666667
        EXMATX(2,5)=0.83333333
        EXMATX(2,6)=0.83333333
        EXMATX(2,7)=-0.16666667
        EXMATX(2,8)=-0.16666667
        EXMATX(2,9)=0.83333333
        EXMATX(2,10)=-0.16666667
C
        EXMATX(3,1)=-0.16666667
        EXMATX(3,2)=-0.16666667
        EXMATX(3,3)=1.83333333
        EXMATX(3,4)=-0.16666667
        EXMATX(3,5)=-0.16666667
        EXMATX(3,6)=0.83333333
        EXMATX(3,7)=0.83333333
        EXMATX(3,8)=-0.16666667
        EXMATX(3,9)=-0.16666667
        EXMATX(3,10)=0.83333333
C
        EXMATX(4,1)=-0.16666667
        EXMATX(4,2)=-0.16666667
        EXMATX(4,3)=-0.16666667
        EXMATX(4,4)=1.83333333
        EXMATX(4,5)=-0.16666667
        EXMATX(4,6)=-0.16666667
        EXMATX(4,7)=-0.16666667
        EXMATX(4,8)=0.83333333
        EXMATX(4,9)=0.83333333
        EXMATX(4,10)=0.83333333
      ENDIF
C
      RETURN
      END
CDOC END_SUBROUTINE EXT7
