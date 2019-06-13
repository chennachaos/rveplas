CDOC BEGIN_SUBROUTINE ARRGO2
CDOC Arrange a fourth order tensor in matrix form with G matrix ordering
CDOC
CDOC This routine re-arranges a given fourth order tensor, stored as a
CDOC 4-index array, in matrix form (2-index array) using G matrix
CDOC component ordering.
CDOC Implemented only for 2-D tensors.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION A4TH   >  Fourth order tensor stored as a 4-index
CDOC C                          array.
CDOC DOUBLE_PRECISION AMATX  <  2-index array containing the components
CDOC C                          of the given 4th order tensor stored
CDOC C                          using G matrix ordering.
CDOC END_PARAMETERS
CHST 
CHST E.de Souza Neto, January  1999: Initial coding
CHST
      SUBROUTINE ARRGO2
     1(   A4TH       ,AMATX      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER
     1(   NDIM=2     ,NGDIM=4    )
C Arguments
      DIMENSION
     1    A4TH(NDIM,NDIM,NDIM,NDIM)              ,AMATX(NGDIM,NGDIM)
C***********************************************************************
C RE-ARRANGES A FOURTH ORDER TENSOR, STORED AS A 4-INDEX ARRAY, IN
C MATRIX FORM (2-INDEX ARRAY) USING G-MATRIX COMPONENT ORDERING
C (11,21,12,22). FOR 2-D ONLY.
C
C REFERENCE: Section D.2.1
C***********************************************************************
C
      AMATX(1,1)=A4TH(1,1,1,1) 
      AMATX(1,2)=A4TH(1,1,2,1) 
      AMATX(1,3)=A4TH(1,1,1,2)
      AMATX(1,4)=A4TH(1,1,2,2) 
C
      AMATX(2,1)=A4TH(2,1,1,1) 
      AMATX(2,2)=A4TH(2,1,2,1) 
      AMATX(2,3)=A4TH(2,1,1,2) 
      AMATX(2,4)=A4TH(2,1,2,2) 
C
      AMATX(3,1)=A4TH(1,2,1,1) 
      AMATX(3,2)=A4TH(1,2,2,1) 
      AMATX(3,3)=A4TH(1,2,1,2) 
      AMATX(3,4)=A4TH(1,2,2,2) 
C
      AMATX(4,1)=A4TH(2,2,1,1) 
      AMATX(4,2)=A4TH(2,2,2,1) 
      AMATX(4,3)=A4TH(2,2,1,2) 
      AMATX(4,4)=A4TH(2,2,2,2) 
C
      RETURN
      END
CDOC END_SUBROUTINE ARRGO2
      