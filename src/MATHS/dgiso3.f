CDOC BEGIN_SUBROUTINE DGISO3
CDOC Derivative of a general isotropic tensor function of one tensor
CDOC
CDOC This function computes the derivative, dY(X)/dX, of a general
CDOC isotropic tensor function of one tensor, Y(X). This implementation
CDOC is restricted to 3-D. The tensor function Y(X) is assumed to be 
CDOC defined as Y(X)= Sum[yi(x1,x2,x3) ei(x)ei], where yi are the 
CDOC eigenvalues of the tensor Y and xi the eigenvalues of the tensor X.
CDOC ei are the eigenvectors of X (which by definition of Y(X), coincide 
CDOC with those of tensor Y) and "(x)" denotes the tensor product.
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
CDOC LOGICAL          REPEAT >  Array of repeated eigenvalues flags.
CDOC C                          If the in-plane eigenvalues are repeated
CDOC C                          this argument must be set to
CDOC C                          .TRUE. on entry, so that the
CDOC C                          appropriate limit expression for the
CDOC C                          the derivative is employed.
CDOC DOUBLE_PRECISION X      >  Function argument (in vector form)
CDOC END_PARAMETERS     
CHST
CHST E.de Souza Neto, October 1994: Initial coding as DKIOG3 (Ogden 
CHST                                material specific function)
CHST E.de Souza Neto,   March 1999: Coded as a general function
CHST
      SUBROUTINE DGISO3
     1(   DEIGY      ,DYDX       ,EIGPRJ     ,EIGX       ,EIGY       ,
     2    REPEAT     ,X          )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER
     1(   MCOMP=6    ,NDIM=3     )
      LOGICAL   REPEAT(5)
      DIMENSION
     1    DEIGY(NDIM,NDIM)   ,DYDX(MCOMP,MCOMP)  ,EIGPRJ(MCOMP,NDIM) ,
     2    EIGX(NDIM)         ,EIGY(NDIM)         ,X(MCOMP)
C
      LOGICAL DIF123, REP123,  REP12,  REP13,  REP23
      DIMENSION
     1    DEIGDX(MCOMP,MCOMP,NDIM)               ,DX2DX(MCOMP,MCOMP) ,
     2    SOID(MCOMP)        ,FOID(MCOMP,MCOMP)
      DATA
     1    FOID(1,1),FOID(1,2),FOID(1,3),FOID(1,4),FOID(1,5),FOID(1,6)  /
     2    1.0D0    ,0.0D0    ,0.0D0    ,0.0D0    ,0.0D0    ,0.0D0      /
     3    FOID(2,1),FOID(2,2),FOID(2,3),FOID(2,4),FOID(2,5),FOID(2,6)  /
     4    0.0D0    ,1.0D0    ,0.0D0    ,0.0D0    ,0.0D0    ,0.0D0      /
     5    FOID(3,1),FOID(3,2),FOID(3,3),FOID(3,4),FOID(3,5),FOID(3,6)  /
     6    0.0D0    ,0.0D0    ,1.0D0    ,0.0D0    ,0.0D0    ,0.0D0      /
     7    FOID(4,1),FOID(4,2),FOID(4,3),FOID(4,4),FOID(4,5),FOID(4,6)  /
     8    0.0D0    ,0.0D0    ,0.0D0    ,0.5D0    ,0.0D0    ,0.0D0      /
     9    FOID(5,1),FOID(5,2),FOID(5,3),FOID(5,4),FOID(5,5),FOID(5,6)  /
     O    0.0D0    ,0.0D0    ,0.0D0    ,0.0D0    ,0.5D0    ,0.0D0      /
     1    FOID(6,1),FOID(6,2),FOID(6,3),FOID(6,4),FOID(6,5),FOID(6,6)  /
     2    0.0D0    ,0.0D0    ,0.0D0    ,0.0D0    ,0.0D0    ,0.5D0      /
      DATA
     1    SOID(1)  ,SOID(2)  ,SOID(3)  ,SOID(4)  ,SOID(5)  ,SOID(6)    /
     2    1.0D0    ,1.0D0    ,1.0D0    ,0.0D0    ,0.D0     ,0.D0       /
      DATA
     1    R0   ,R2   ,R3   /
     2    0.0D0,2.0D0,3.0D0/
C***********************************************************************
C COMPUTE THE DERIVATIVE OF A GENERAL ISOTROPIC TENSOR FUNCTION OF ONE
C TENSOR IN 3-D
C
C REFERENCE: Sections A.4.1-2
C            Box A.6
C***********************************************************************
C Retrieve information on repeated eigenvalues
C ============================================
      DIF123=REPEAT(1)
      REP123=REPEAT(2)
      REP12=REPEAT(3)
      REP13=REPEAT(4)
      REP23=REPEAT(5)
C Evaluate DYDX
C =============
      CALL RVZERO(DYDX,MCOMP*MCOMP)
      IF(REP123)THEN
C Three repeated eigenvalues
C --------------------------
        A0=DEIGY(1,1)-DEIGY(1,2)
        A1=DEIGY(1,2)
        DO 110 I=1,MCOMP
          DO 100 J=1,MCOMP
            DYDX(I,J)=DYDX(I,J)+A0*FOID(I,J)+A1*SOID(I)*SOID(J)
 100      CONTINUE
 110    CONTINUE
      ELSE
C Set components of the derivative of X**2 with respect to X
        DX2DX(1,1)=R2*X(1)
        DX2DX(1,2)=R0
        DX2DX(1,3)=R0
        DX2DX(1,4)=X(4)
        DX2DX(1,5)=R0
        DX2DX(1,6)=X(6)
        DX2DX(2,2)=R2*X(2)
        DX2DX(2,3)=R0
        DX2DX(2,4)=X(4)
        DX2DX(2,5)=X(5)
        DX2DX(2,6)=R0
        DX2DX(3,3)=R2*X(3)
        DX2DX(3,4)=R0
        DX2DX(3,5)=X(5)
        DX2DX(3,6)=X(6)
        DX2DX(4,4)=(X(2)+X(1))/R2
        DX2DX(4,5)=X(6)/R2
        DX2DX(4,6)=X(5)/R2
        DX2DX(5,5)=(X(3)+X(2))/R2
        DX2DX(5,6)=X(4)/R2
        DX2DX(6,6)=(X(3)+X(1))/R2
        DO 130 I=2,MCOMP
          DO 120 J=1,I-1
            DX2DX(I,J)=DX2DX(J,I)
 120      CONTINUE
 130    CONTINUE
        IF(DIF123) THEN
C Three distinct eigenvalues
C --------------------------
C Set components of the derivatives of the eigenprojections
          CALL RVZERO(DEIGDX,MCOMP*MCOMP*NDIM)
          DO 160 IDIR=1,3
            IF(IDIR.EQ.1)THEN
              IA=2
              IB=3
            ELSEIF(IDIR.EQ.2)THEN
              IA=3
              IB=1
            ELSEIF(IDIR.EQ.3)THEN
              IA=1
              IB=2
            ENDIF
            A0=EIGX(IDIR)
            B1=X(1)+X(2)+X(3)
            A1=B1-A0
            A2=B1-R3*A0
            A3=(A0-EIGX(IA))*(A0-EIGX(IB))
            DO 150 I=1,MCOMP
              DO 140 J=1,MCOMP
              DEIGDX(I,J,IDIR)=(DX2DX(I,J)-A1*FOID(I,J)+
     1                         A2*EIGPRJ(I,IDIR)*EIGPRJ(J,IDIR)-
     2                         (EIGX(IA)-EIGX(IB))*
     3         (EIGPRJ(I,IA)*EIGPRJ(J,IA)-EIGPRJ(I,IB)*EIGPRJ(J,IB)))/A3
 140          CONTINUE
 150        CONTINUE
 160      CONTINUE
C Compute components of DYDX
          DO 190 IDIR=1,3
            DO 180 I=1,MCOMP
              DO 170 J=1,MCOMP
                DYDX(I,J)=DYDX(I,J)+EIGY(IDIR)*DEIGDX(I,J,IDIR)
 170          CONTINUE
 180        CONTINUE
 190      CONTINUE
          DO 230 IDIR=1,3
            DO 220 JDIR=1,3
              A1=DEIGY(IDIR,JDIR)
              DO 210 I=1,MCOMP
                DO 200 J=1,MCOMP
                  DYDX(I,J)=DYDX(I,J)+A1*EIGPRJ(I,IDIR)*EIGPRJ(J,JDIR)
 200            CONTINUE
 210          CONTINUE
 220        CONTINUE
 230      CONTINUE
        ELSE
C Two repeated eigenvalues
C ------------------------
C Compute components of DYDX
          IF(REP12)THEN
            IA=3
            IB=1
            IC=2
          ELSEIF(REP13)THEN
            IA=2
            IB=3
            IC=1
          ELSEIF(REP23)THEN
            IA=1
            IB=2
            IC=3
          ENDIF
          F1=EIGY(IA)
          F3=EIGY(IC)
          F11=DEIGY(IA,IA)
          F13=DEIGY(IA,IC)
          F31=DEIGY(IC,IA)
          F32=DEIGY(IC,IB)
          F33=DEIGY(IC,IC)
          X1=EIGX(IA)
          X3=EIGX(IC)
          X3SQ=X3*X3
          A1=X1-X3
          A1SQ=A1*A1
          A0=R2/(A1SQ*A1)
          A2=X1+X3
          A3=X3/A1SQ
          S1=(F1-F3)/A1SQ+(F32-F33)/A1
          S2=R2*A3*(F1-F3)+A2/A1*(F32-F33)
          S3=-(A0*(F3-F1)+(F33+F11-F31-F13)/A1SQ)
          S4=A0*X3*(F1-F3)+(F13-F32)/A1+A3*(F13+F31-F11-F33)
          S5=A0*X3*(F1-F3)+(F31-F32)/A1+A3*(F13+F31-F11-F33)
          S6=A0*X3SQ*(F3-F1)-X1*A3*(F31+F13)+
     1                         A2/A1*F32+X3SQ/A1SQ*(F11+F33)
          DO 250 I=1,MCOMP
            DO 240 J=1,MCOMP
             DYDX(I,J)=S1*DX2DX(I,J)-S2*FOID(I,J)-S3*X(I)*X(J)+
     1                 S4*X(I)*SOID(J)+S5*SOID(I)*X(J)+
     2                 S6*SOID(I)*SOID(J)
 240        CONTINUE
 250      CONTINUE
        ENDIF
      ENDIF
      RETURN
      END
CDOC END_SUBROUTINE DGISO3