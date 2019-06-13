CDOC BEGIN_SUBROUTINE ARRGAX
CDOC Arrange a fourth order tensor in matrix form with G matrix ordering
CDOC
CDOC This routine re-arranges a given fourth order tensor, stored as a
CDOC 4-index array, in matrix form (2-index array) using G matrix
CDOC component ordering.
CDOC Implemented only for axisymmetric problems.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION A4TH   >  Fourth order tensor stored as a 4-index
CDOC C                          array.
CDOC DOUBLE_PRECISION AMATX  <  2-index array containing the components
CDOC C                          of the given 4th order tensor stored
CDOC C                          using G matrix ordering.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, May  1999: Initial coding
CDOC
      SUBROUTINE ARRGAX
     1(   A4TH       ,AMATX      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER
     1(   NDIM=3     ,NGDIM=5    )
C Arguments
      DIMENSION
     1    A4TH(NDIM,NDIM,NDIM,NDIM)              ,AMATX(NGDIM,NGDIM)
C***********************************************************************
C RE-ARRANGES A FOURTH ORDER TENSOR, STORED AS A 4-INDEX ARRAY, IN
C MATRIX FORM (2-INDEX ARRAY) USING G MATRIX COMPONENT ORDERING
C (11,21,12,22,33). FOR AXISYMMETRIC CASE ONLY.
C
C REFERENCE: Section D.2.1
C
C***********************************************************************
C
      AMATX(1,1)=A4TH(1,1,1,1) 
      AMATX(1,2)=A4TH(1,1,2,1) 
      AMATX(1,3)=A4TH(1,1,1,2)
      AMATX(1,4)=A4TH(1,1,2,2) 
      AMATX(1,5)=A4TH(1,1,3,3) 
C
      AMATX(2,1)=A4TH(2,1,1,1) 
      AMATX(2,2)=A4TH(2,1,2,1) 
      AMATX(2,3)=A4TH(2,1,1,2) 
      AMATX(2,4)=A4TH(2,1,2,2) 
      AMATX(2,5)=A4TH(2,1,3,3) 
C
      AMATX(3,1)=A4TH(1,2,1,1) 
      AMATX(3,2)=A4TH(1,2,2,1) 
      AMATX(3,3)=A4TH(1,2,1,2) 
      AMATX(3,4)=A4TH(1,2,2,2) 
      AMATX(3,5)=A4TH(1,2,3,3) 
C
      AMATX(4,1)=A4TH(2,2,1,1) 
      AMATX(4,2)=A4TH(2,2,2,1) 
      AMATX(4,3)=A4TH(2,2,1,2) 
      AMATX(4,4)=A4TH(2,2,2,2) 
      AMATX(4,5)=A4TH(2,2,3,3) 
C
      AMATX(5,1)=A4TH(3,3,1,1) 
      AMATX(5,2)=A4TH(3,3,2,1) 
      AMATX(5,3)=A4TH(3,3,1,2) 
      AMATX(5,4)=A4TH(3,3,2,2) 
      AMATX(5,5)=A4TH(3,3,3,3) 
C
      RETURN
      END
CDOC END_SUBROUTINE ARRGAX
