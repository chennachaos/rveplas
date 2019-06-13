CDOC BEGIN_SUBROUTINE DGISO2
CDOC Derivative of a general isotropic tensor function of one tensor
CDOC
CDOC This function computes the derivative, dY(X)/dX, of a general
CDOC isotropic tensor function of one tensor, Y(X). This implementation
CDOC is restricted to 2-D with one possible out-of-plane component
CDOC (normally needed in axisymmetric problems). The tensor function
CDOC Y(X) is assumed to be defined as Y(X)= Sum[yi(x1,x2,x3) ei(x)ei],
CDOC where yi are the eigenvalues of the tensor Y and xi the eigenvalues
CDOC of the tensor X. ei are the egenvectors of X (which by definition
CDOC of Y(X), coincide with those of tensor Y) and "(x)" denotes the
CDOC tensor product.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DEIGY  >  Matrix containing the derivatives
CDOC C                          dyi/dxj of the eigenvalues of Y(X) with
CDOC C                          respect to the eigenvalues of X.
CDOC DOUBLE_PRECISION DYDX   <  Matrix of components of the derivative
CDOC C                          (fourth order tensor) dY/dX.
CDOC DOUBLE_PRECISION EIGPRJ >  Matrix with each column containing the
CDOC C                          components of one eigenprojection tensor
CDOC C                          ei(x)ei.
CDOC DOUBLE_PRECISION EIGX   >  Array of eigenvalues of X.
CDOC DOUBLE_PRECISION EIGY   >  Array of eigenvalues of Y.
CDOC LOGICAL          OUTOFP >  Out-of-plane component flag. If set to
CDOC C                          .TRUE. the out-of-plane
CDOC C                          component (normally required in
CDOC C                          axisymmetric problems) is computed.
CDOC C                          The out-of-plane component is not
CDOC C                          computed otherwise.
CDOC LOGICAL          REPEAT >  Repeated in-plane eigenvalues flag.
CDOC C                          If the in-plane eigenvalues are repeated
CDOC C                          this argument must be set to
CDOC C                          .TRUE. on entry, so that the
CDOC C                          appropriate limit expression for the
CDOC C                          the derivative is employed.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, May 1996: Initial coding
CHST
      SUBROUTINE DGISO2
     1(   DEIGY      ,DYDX       ,EIGPRJ     ,EIGX       ,EIGY       ,
     2    OUTOFP     ,REPEAT     )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER
     1(   MCOMP=4    ,MDIM=3     ,NDIM=2     )
C Arguments
      LOGICAL OUTOFP, REPEAT
      DIMENSION
     1    DEIGY(MDIM,MDIM)  ,DYDX(MCOMP,MCOMP) ,EIGPRJ(MCOMP,NDIM),
     2    EIGX(NDIM)        ,EIGY(NDIM)
C Local arrays
      DIMENSION
     1    EIGPR3(MCOMP)     ,FOID(MCOMP,MCOMP) ,SOPID(MCOMP)
      DATA
     1    FOID(1,1)     ,FOID(1,2)     ,FOID(1,3)     /
     2    1.0D0         ,0.0D0         ,0.0D0         /
     3    FOID(2,1)     ,FOID(2,2)     ,FOID(2,3)     /
     4    0.0D0         ,1.0D0         ,0.0D0         /
     5    FOID(3,1)     ,FOID(3,2)     ,FOID(3,3)     /
     6    0.0D0         ,0.0D0         ,0.5D0         /
      DATA
     1    SOPID(1)      ,SOPID(2)      ,SOPID(3)      ,SOPID(4)        /
     2    1.0D0         ,1.0D0         ,0.0D0         ,0.0D0           /
      DATA
     1    EIGPR3(1)     ,EIGPR3(2)     ,EIGPR3(3)     ,EIGPR3(4)       /
     2    0.0D0         ,0.0D0         ,0.0D0         ,1.0D0           /
C***********************************************************************
C COMPUTE THE DERIVATIVE OF A GENERAL ISOTROPIC TENSOR FUNCTION OF ONE
C TENSOR IN 2-D (WITH ONE POSSIBLE OUT-OF-PLANE COMPONENT)
C
C REFERENCE: Sections A.3.1-2
C            Box A.2
C***********************************************************************
      CALL RVZERO(DYDX,MCOMP*MCOMP)
      IF(REPEAT)THEN
C Derivative dY/dX for repeated in-plane eigenvalues of X
C -------------------------------------------------------
C In-plane component
        DO 20 I=1,3
          DO 10 J=1,3
            DYDX(I,J)=(DEIGY(1,1)-DEIGY(1,2))*FOID(I,J)+
     1                 DEIGY(1,2)*SOPID(I)*SOPID(J)
   10     CONTINUE
   20   CONTINUE
        IF(OUTOFP)THEN
C out-of-plane components required
          DO 40 I=1,4
            DO 30 J=1,4
              IF(I.EQ.4.OR.J.EQ.4)DYDX(I,J)=
     1                    DEIGY(1,3)*SOPID(I)*EIGPR3(J)+
     2                    DEIGY(3,1)*EIGPR3(I)*SOPID(J)+
     3                    DEIGY(3,3)*EIGPR3(I)*EIGPR3(J)
   30       CONTINUE
   40     CONTINUE
        ENDIF
      ELSE
C Derivative dY/dX for distinct in-plane eigenvalues of X
C -------------------------------------------------------
C Assemble in-plane DYDX
        A1=(EIGY(1)-EIGY(2))/(EIGX(1)-EIGX(2))
        DO 70 I=1,3
          DO 60 J=1,3
            DYDX(I,J)=A1*(FOID(I,J)-EIGPRJ(I,1)*EIGPRJ(J,1)-
     1                EIGPRJ(I,2)*EIGPRJ(J,2))+
     2                DEIGY(1,1)*EIGPRJ(I,1)*EIGPRJ(J,1)+
     3                DEIGY(1,2)*EIGPRJ(I,1)*EIGPRJ(J,2)+
     4                DEIGY(2,1)*EIGPRJ(I,2)*EIGPRJ(J,1)+
     5                DEIGY(2,2)*EIGPRJ(I,2)*EIGPRJ(J,2)
   60     CONTINUE
   70   CONTINUE
        IF(OUTOFP) THEN
C out-of-plane components required
          DO 90 I=1,4
            DO 80 J=1,4
              IF(I.EQ.4.OR.J.EQ.4)DYDX(I,J)=
     1                DEIGY(1,3)*EIGPRJ(I,1)*EIGPR3(J)+
     2                DEIGY(2,3)*EIGPRJ(I,2)*EIGPR3(J)+
     3                DEIGY(3,1)*EIGPR3(I)*EIGPRJ(J,1)+
     4                DEIGY(3,2)*EIGPR3(I)*EIGPRJ(J,2)+
     5                DEIGY(3,3)*EIGPR3(I)*EIGPR3(J)
   80       CONTINUE
   90     CONTINUE
        ENDIF
      ENDIF
      RETURN
      END
CDOC END_SUBROUTINE DGISO2
