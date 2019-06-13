CDOC BEGIN_SUBROUTINE ISO2
CDOC Computes the value of isotropic tensor functions of one tensor.
CDOC
CDOC This subroutine evaluates isotropic tensor functions Y(X), of one
CDOC tensor belonging to the class described below.
CDOC This implementation is restricted to 2-D with one possible
CDOC out-of-plane component (normally needed in axisymmetric problems).
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
CDOC LOGICAL          OUTOFP >  Out-of-plane component flag. If set to
CDOC C                          .TRUE. the out-of-plane
CDOC C                          component (normally required in
CDOC C                          axisymmetric problems) is computed.
CDOC C                          The out-of-plane component is not
CDOC C                          computed otherwise.
CDOC DOUBLE_PRECISION X      >  Array of components of the tensor at
CDOC C                          which the function is to be evaluated.
CDOC DOUBLE_PRECISION Y      <  Array of components of the tensor
CDOC C                          function at X.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, August   1996: Initial coding
CHST E.de Souza Neto, February 2004: External declaration for argument
CHST                                 FUNC removed
CHST
      SUBROUTINE ISO2
     1(   FUNC       ,OUTOFP     ,X          ,Y          )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER
     1(   MCOMP=4    ,NDIM=2     )
      LOGICAL OUTOFP ,REPEAT
      DIMENSION
     1    X(*)                      ,Y(*)
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
C FUNCTION. THIS ROUTINE IS RESTRICTED TO 2-D TENSORS WITH ONE
C POSSIBLE (TRANSVERSAL) OUT-OF-PLANE COMPONENT.
C
C REFERENCE: Section A.5
C***********************************************************************
C Performs the spectral decomposition of X
      CALL SPDEC2
     1(   EIGPRJ     ,EIGX       ,REPEAT     ,X          )
C Computes the in-plane eigenvalues of Y
      DO 10 IDIR=1,2
        EIGY(IDIR)=FUNC(EIGX(IDIR))
   10 CONTINUE
C Assembles in-plane component of Y (in vector form)
      CALL RVZERO(Y,3)
      DO 30 ICOMP=1,3
        DO 20 IDIR=1,2
          Y(ICOMP)=Y(ICOMP)+EIGY(IDIR)*EIGPRJ(ICOMP,IDIR)
   20   CONTINUE
   30 CONTINUE
C Out-of-plane component required
      IF(OUTOFP)Y(4)=FUNC(X(4))
C
      RETURN
      END
CDOC END_SUBROUTINE ISO2
