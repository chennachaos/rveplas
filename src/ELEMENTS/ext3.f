CDOC BEGIN_SUBROUTINE EXT3
CDOC Gauss point-node extrapolation matrix for element type TRI3
CDOC
CDOC This routine sets the coefficients matrix for extrapolation of
CDOC fields from Gauss point values to nodal values for element type
CDOC TRI3: Standard isoparametric 3-noded linear triangle.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION EXMATX <  Extrapolation matrix.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, February 1997: Initial coding
CHST
      SUBROUTINE EXT3
     1(   EXMATX     )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NNODE=3,NGAUSP=1)
      DIMENSION EXMATX(NNODE,NGAUSP)
      DATA R1    /
     1     1.0D0 /
C***********************************************************************
C SET COEFFICIENTS MATRIX (EXMATX) FOR EXTRAPOLATION FROM GAUSS POINTS
C TO NODES FOR ELEMENT TYPE 'TRI_3' (STANDARD 3-NODED LINEAR TRIANGLE)
C
C REFERENCE: Section 5.6.1
C            E Hinton & JS Campbel. Local and global Smoothing of
C            discontinuous finite element functions using a least
C            squares method. Int. J. Num. meth. Engng., 8:461-480, 1974.
C            E Hinton & DRJ Owen. An introduction to finite element
C            computations. Pineridge Press, Swansea, 1979.
C***********************************************************************
      EXMATX(1,1)=R1
      EXMATX(2,1)=R1
      EXMATX(3,1)=R1
C
      RETURN
      END
CDOC END_SUBROUTINE EXT3
