CDOC BEGIN_SUBROUTINE ISO3
CDOC Computes the value of isotropic tensor functions of one tensor.
CDOC
CDOC This subroutine evaluates isotropic tensor functions Y(X), of one
CDOC tensor belonging to the class described below.
CDOC This implementation is restricted to 3-D.
CDOC The class of symmetric tensor functions Y(X) is assumed to be
CDOC defined as Y(X)= Sum[y(xi) ei(x)ei], where the scalar function
CDOC y(xi) defines the eigenvalues of the tensor Y and xi the
CDOC eigenvalues of the tensor X. ei are the egenvectors of X (which by
CDOC definition of Y(X), coincide with those of tensor Y) and "(x)"
CDOC denotes the tensor product.
CDOC
CDOC BEGIN_PARAMETERS
CDOC SYMBOLIC_NAME    FUNC   >  Symbolic name of the double
CDOC C                          precision function defining y(xi).
CDOC DOUBLE_PRECISION X      >  Array of components of the tensor at
CDOC C                          which the function is to be evaluated.
CDOC DOUBLE_PRECISION Y      <  Array of components of the tensor
CDOC C                          function at X.
CDOC END_PARAMETERS
CHST
CHST D. de Bortoli, March 2015: Initial coding, based on ISO2
CHST
      SUBROUTINE ISO3
     1(   FUNC       ,X          ,Y          )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER
     1(   MCOMP=6    ,NDIM=3     )
      LOGICAL REPEAT(5)
      DIMENSION
     1    X(MCOMP)                  ,Y(MCOMP)
      DIMENSION
     1    EIGPRJ(MCOMP,NDIM)        ,EIGX(NDIM)                ,
     1    EIGY(NDIM)
C***********************************************************************
C COMPUTE THE TENSOR Y (STORED IN VECTOR FORM) AS AN ISOTROPIC
C FUNCTION OF THE TYPE:
C
C                     Y(X) = sum{ y(x_i) E_i }
C
C WHERE Y AND X ARE SYMMETRIC TENSORS, x_i AND E_i ARE, RESPECTIVELY
C THE EIGENVALUES AND EIGENPROJECTIONS OF X, AND y(.) IS A SCALAR
C FUNCTION. THIS ROUTINE IS RESTRICTED TO 3-D.
C
C REFERENCE: Section A.5
C***********************************************************************
C Performs the spectral decomposition of X
      CALL SPDEC3
     1(   EIGPRJ     ,EIGX       ,REPEAT     ,X          )
C Computes the eigenvalues of Y
      DO 10 IDIR=1,NDIM
        EIGY(IDIR)=FUNC(EIGX(IDIR))
   10 CONTINUE
C Assembles Y (in vector form)
      CALL RVZERO(Y,MCOMP)
      DO 30 ICOMP=1,MCOMP
        DO 20 IDIR=1,NDIM
          Y(ICOMP)=Y(ICOMP)+EIGY(IDIR)*EIGPRJ(ICOMP,IDIR)
   20   CONTINUE
   30 CONTINUE
      RETURN
      END
CDOC END_SUBROUTINE ISO3
