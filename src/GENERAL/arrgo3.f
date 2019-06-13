CDOC BEGIN_SUBROUTINE ARRGO3
CDOC Arrange a fourth order tensor in matrix form with G matrix ordering
CDOC
CDOC This routine re-arranges a given fourth order tensor, stored as a
CDOC 4-index array, in matrix form (2-index array) using G matrix
CDOC component ordering.
CDOC Implemented only for 3-D tensors.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION A4TH   >  Fourth order tensor stored as a 4-index
CDOC C                          array.
CDOC DOUBLE_PRECISION AMATX  <  2-index array containing the components
CDOC C                          of the given 4th order tensor stored
CDOC C                          using G matrix ordering.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, January  1999: Initial coding
CDOC
      SUBROUTINE ARRGO3
     1(   A4TH       ,AMATX      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER
     1(   NDIM=3     ,NGDIM=9    )
C Arguments
      DIMENSION
     1    A4TH(NDIM,NDIM,NDIM,NDIM)              ,AMATX(NGDIM,NGDIM)
C***********************************************************************
C RE-ARRANGES A FOURTH ORDER TENSOR, STORED AS A 4-INDEX ARRAY, IN
C MATRIX FORM (2-INDEX ARRAY) USING G MATRIX COMPONENT ORDERING
C C (11,21,31,12,22,32,13,23,33). FOR 3-D ONLY.
C***********************************************************************
C
      AMATX(1,1)=A4TH(1,1,1,1) 
      AMATX(1,2)=A4TH(1,1,2,1) 
      AMATX(1,3)=A4TH(1,1,3,1)
      AMATX(1,4)=A4TH(1,1,1,2) 
      AMATX(1,5)=A4TH(1,1,2,2) 
      AMATX(1,6)=A4TH(1,1,3,2) 
      AMATX(1,7)=A4TH(1,1,1,3) 
      AMATX(1,8)=A4TH(1,1,2,3) 
      AMATX(1,9)=A4TH(1,1,3,3) 
C
      AMATX(2,1)=A4TH(2,1,1,1) 
      AMATX(2,2)=A4TH(2,1,2,1) 
      AMATX(2,3)=A4TH(2,1,3,1)
      AMATX(2,4)=A4TH(2,1,1,2) 
      AMATX(2,5)=A4TH(2,1,2,2) 
      AMATX(2,6)=A4TH(2,1,3,2) 
      AMATX(2,7)=A4TH(2,1,1,3) 
      AMATX(2,8)=A4TH(2,1,2,3) 
      AMATX(2,9)=A4TH(2,1,3,3) 
C
      AMATX(3,1)=A4TH(3,1,1,1) 
      AMATX(3,2)=A4TH(3,1,2,1) 
      AMATX(3,3)=A4TH(3,1,3,1)
      AMATX(3,4)=A4TH(3,1,1,2) 
      AMATX(3,5)=A4TH(3,1,2,2) 
      AMATX(3,6)=A4TH(3,1,3,2) 
      AMATX(3,7)=A4TH(3,1,1,3) 
      AMATX(3,8)=A4TH(3,1,2,3) 
      AMATX(3,9)=A4TH(3,1,3,3) 
C
      AMATX(4,1)=A4TH(1,2,1,1) 
      AMATX(4,2)=A4TH(1,2,2,1) 
      AMATX(4,3)=A4TH(1,2,3,1)
      AMATX(4,4)=A4TH(1,2,1,2) 
      AMATX(4,5)=A4TH(1,2,2,2) 
      AMATX(4,6)=A4TH(1,2,3,2) 
      AMATX(4,7)=A4TH(1,2,1,3) 
      AMATX(4,8)=A4TH(1,2,2,3) 
      AMATX(4,9)=A4TH(1,2,3,3) 
C
      AMATX(5,1)=A4TH(2,2,1,1) 
      AMATX(5,2)=A4TH(2,2,2,1) 
      AMATX(5,3)=A4TH(2,2,3,1)
      AMATX(5,4)=A4TH(2,2,1,2) 
      AMATX(5,5)=A4TH(2,2,2,2) 
      AMATX(5,6)=A4TH(2,2,3,2) 
      AMATX(5,7)=A4TH(2,2,1,3) 
      AMATX(5,8)=A4TH(2,2,2,3) 
      AMATX(5,9)=A4TH(2,2,3,3) 
C
      AMATX(6,1)=A4TH(3,2,1,1) 
      AMATX(6,2)=A4TH(3,2,2,1) 
      AMATX(6,3)=A4TH(3,2,3,1)
      AMATX(6,4)=A4TH(3,2,1,2) 
      AMATX(6,5)=A4TH(3,2,2,2) 
      AMATX(6,6)=A4TH(3,2,3,2) 
      AMATX(6,7)=A4TH(3,2,1,3) 
      AMATX(6,8)=A4TH(3,2,2,3) 
      AMATX(6,9)=A4TH(3,2,3,3) 
C
      AMATX(7,1)=A4TH(1,3,1,1) 
      AMATX(7,2)=A4TH(1,3,2,1) 
      AMATX(7,3)=A4TH(1,3,3,1)
      AMATX(7,4)=A4TH(1,3,1,2) 
      AMATX(7,5)=A4TH(1,3,2,2) 
      AMATX(7,6)=A4TH(1,3,3,2) 
      AMATX(7,7)=A4TH(1,3,1,3) 
      AMATX(7,8)=A4TH(1,3,2,3) 
      AMATX(7,9)=A4TH(1,3,3,3) 
C
      AMATX(8,1)=A4TH(2,3,1,1) 
      AMATX(8,2)=A4TH(2,3,2,1) 
      AMATX(8,3)=A4TH(2,3,3,1)
      AMATX(8,4)=A4TH(2,3,1,2) 
      AMATX(8,5)=A4TH(2,3,2,2) 
      AMATX(8,6)=A4TH(2,3,3,2) 
      AMATX(8,7)=A4TH(2,3,1,3) 
      AMATX(8,8)=A4TH(2,3,2,3) 
      AMATX(8,9)=A4TH(2,3,3,3) 
C
      AMATX(9,1)=A4TH(3,3,1,1) 
      AMATX(9,2)=A4TH(3,3,2,1) 
      AMATX(9,3)=A4TH(3,3,3,1)
      AMATX(9,4)=A4TH(3,3,1,2) 
      AMATX(9,5)=A4TH(3,3,2,2) 
      AMATX(9,6)=A4TH(3,3,3,2) 
      AMATX(9,7)=A4TH(3,3,1,3) 
      AMATX(9,8)=A4TH(3,3,2,3) 
      AMATX(9,9)=A4TH(3,3,3,3) 
C
C
      RETURN
      END
CDOC END_SUBROUTINE ARRGO3
