CDOC BEGIN_SUBROUTINE DISO3
CDOC Derivative of a class of isotropic tensor function of one tensor
CDOC
CDOC This function computes the derivative, dY(X)/dX, of a particular
CDOC class of symmetric isotropic tensor valued function of one
CDOC symmetric tensor, Y(X).
CDOC This implementation is restricted to 3-D.
CDOC The subroutine for the general derivative of isotropic tensor
CDOC functions of one tensor, DGISO3, is particularised to this class
CDOC of tensor functions.
CDOC The class of tensor functions Y(X) is assumed to be defined as
CDOC Y(X)= Sum[y(xi) ei(x)ei], where the scalar function y(xi) defines
CDOC the eigenvalues of the tensor Y and xi the eigenvalues of the
CDOC tensor X. ei are the egenvectors of X (which by definition of
CDOC Y(X), coincide with those of tensor Y) and "(x)" denotes the
CDOC tensor product.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DYDX   <  Matrix of components of the derivative
CDOC C                          (fourth order tensor) dY/dX.
CDOC SYMBOLIC_NAME    DFUNC  >  Symbolic name of the double precision
CDOC C                          function defining the derivative
CDOC C                          dy(x)/dx of the eigenvalues of the
CDOC C                          tensor function.
CDOC SYMBOLIC_NAME    FUNC   >  Symbolic name of the double precision
CDOC C                          function defining the eigenvalues y(x)
CDOC C                          of the tensor fumction.
CDOC DOUBLE_PRECISION X      >  Point (argument) at which the
CDOC C                          derivative is to be computed.
CDOC END_PARAMETERS
CHST
CHST D. de Bortoli, March 2015: Initial coding, based on DISO2
CHST
      SUBROUTINE DISO3
     1(   DYDX       ,DFUNC      ,FUNC       ,X          )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL
     1    DFUNC      ,FUNC
      PARAMETER
     1(   MCOMP=6    ,NDIM=3     )
C Arguments
      DIMENSION
     1    DYDX(MCOMP,MCOMP) ,X(MCOMP)
C Local variables and arrays
      LOGICAL REPEAT(5)
      DIMENSION
     1    DEIGY(NDIM,NDIM)  ,EIGPRJ(MCOMP,NDIM)       ,EIGX(NDIM)      ,
     2    EIGY(NDIM)
C***********************************************************************
C COMPUTE (AND STORE IN MATRIX FORM) THE DERIVATIVE dY/dX OF AN
C ISOTROPIC TENSOR FUNCTION OF THE TYPE:
C
C                        Y(X) = sum{ y(x_i) E_i }
C
C WHERE Y AND X ARE SYMMETRIC TENSORS, x_i AND E_i ARE, RESPECTIVELY
C THE EIGENVALUES AND EIGENPROJECTIONS OF X, AND y(.) IS A SCALAR
C FUNCTION. THIS ROUTINE IS RESTRICTED TO 3-D TENSORS.
C
C REFERENCE: Section A.5.2
C***********************************************************************
C Spectral decomposition of X
      CALL SPDEC3
     1(   EIGPRJ     ,EIGX       ,REPEAT     ,X          )
C Evaluate eigenvalues of Y and their derivatives
C -----------------------------------------------
C DEIGY is a diagonal matrix containing the derivatives of the 
C eigenvalues of Y with respect to each eigenvalue of X
      CALL RVZERO(DEIGY, NDIM*NDIM)
      DO 10 IDIR=1,3
        EIGY(IDIR)=FUNC(EIGX(IDIR))
        DEIGY(IDIR,IDIR)=DFUNC(EIGX(IDIR))
   10 CONTINUE
C
C Calculate components of dY/dX using DGISO3
C ------------------------------------------
      CALL DGISO3
     1(   DEIGY      ,DYDX       ,EIGPRJ     ,EIGX       ,EIGY       ,
     2    REPEAT     ,X          )
C
      RETURN
      END
CDOC END_SUBROUTINE DISO3
