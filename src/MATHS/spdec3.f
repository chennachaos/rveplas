CDOC BEGIN_SUBROUTINE SPDEC3
CDOC Spectral decomposition of 3-D symmetric tensors
CDOC
CDOC This routine performs the spectral decomposition of 3-D symmetric
CDOC tensors. The tensor is passed as argument (stored in vector form).
CDOC 
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION EIGPRJ <  Matrix with one eigenprojection tensor
CDOC C                          of \smparm{X} stored in each column.
CDOC DOUBLE_PRECISION EIGX   <  Array containing the eigenvalues of
CDOC C                          \smparm{X}.
CDOC LOGICAL          REPEAT <  Array of repeated (within a small
CDOC C                          tolerance) eigenvalues flag.
CDOC DOUBLE_PRECISION X      >  Array containing the components of a
CDOC C                          symmetric tensor.
CDOC END_PARAMETERS
CDOC
CDOC E.de Souza Neto, March 1999: Initial coding
CDOC
      SUBROUTINE SPDEC3
     1(   EIGPRJ     ,EIGX       ,REPEAT     ,X          )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER
     1(   MCOMP=6    ,NDIM=3     )
C Arguments
      LOGICAL REPEAT(5)
      DIMENSION
     1    EIGPRJ(MCOMP,NDIM)        ,EIGX(NDIM)                ,
     2    X(MCOMP)
C Local arrays and variables
      LOGICAL DIF123 ,REP123 ,REP12 ,REP13 ,REP23
      DIMENSION
     1    AUXMTX(NDIM,NDIM)         ,EIGVEC(NDIM,NDIM)
      DATA
     1    R0   ,R1   ,SMALL  /
     2    0.0D0,1.0D0,1.D-8  /
C***********************************************************************
C PERFORMS THE SPECTRAL DECOMPOSITION OF A SYMMETRIC 3-D TENSOR STORED
C IN VECTOR FORM
C***********************************************************************
      AUXMTX(1,1)=X(1)
      AUXMTX(2,2)=X(2)
      AUXMTX(3,3)=X(3)
      AUXMTX(1,2)=X(4)
      AUXMTX(2,3)=X(5)
      AUXMTX(1,3)=X(6)
      AUXMTX(2,1)=AUXMTX(1,2)
      AUXMTX(3,2)=AUXMTX(2,3)
      AUXMTX(3,1)=AUXMTX(1,3)
      CALL JACOB(AUXMTX,EIGX,EIGVEC,3)
      DO IDIR=1,3
        EIGPRJ(1,IDIR)=EIGVEC(1,IDIR)*EIGVEC(1,IDIR)
        EIGPRJ(2,IDIR)=EIGVEC(2,IDIR)*EIGVEC(2,IDIR)
        EIGPRJ(3,IDIR)=EIGVEC(3,IDIR)*EIGVEC(3,IDIR)
        EIGPRJ(4,IDIR)=EIGVEC(1,IDIR)*EIGVEC(2,IDIR)
        EIGPRJ(5,IDIR)=EIGVEC(2,IDIR)*EIGVEC(3,IDIR)
        EIGPRJ(6,IDIR)=EIGVEC(1,IDIR)*EIGVEC(3,IDIR)
      END DO
C Identify possible repetitions of principal stretches
C ----------------------------------------------------
      DIF123=.FALSE.
      REP123=.FALSE.
      REP12=.FALSE.
      REP13=.FALSE.
      REP23=.FALSE.
      AMXEIG=ABS(EIGX(1))
      IF(ABS(EIGX(2)).GT.AMXEIG)AMXEIG=ABS(EIGX(2))
      IF(ABS(EIGX(3)).GT.AMXEIG)AMXEIG=ABS(EIGX(3))
      FACTOR=R1
      IF(AMXEIG.NE.R0)THEN
        FACTOR=AMXEIG
      ELSE
        FACTOR=R1
      ENDIF
      IF(ABS(EIGX(1)-EIGX(2))/FACTOR.LT.SMALL.AND.
     1   ABS(EIGX(1)-EIGX(3))/FACTOR.LT.SMALL)THEN
        REP123=.TRUE.
      ELSEIF(ABS(EIGX(1)-EIGX(2))/FACTOR.LT.SMALL.AND.
     1       ABS(EIGX(1)-EIGX(3))/FACTOR.GE.SMALL)THEN
        REP12=.TRUE.
      ELSEIF(ABS(EIGX(1)-EIGX(3))/FACTOR.LT.SMALL.AND.
     1       ABS(EIGX(1)-EIGX(2))/FACTOR.GE.SMALL)THEN
        REP13=.TRUE.
      ELSEIF(ABS(EIGX(2)-EIGX(3))/FACTOR.LT.SMALL.AND.
     1       ABS(EIGX(2)-EIGX(1))/FACTOR.GE.SMALL)THEN
        REP23=.TRUE.
      ELSE
        DIF123=.TRUE.
      ENDIF
      REPEAT(1)=DIF123
      REPEAT(2)=REP123
      REPEAT(3)=REP12
      REPEAT(4)=REP13
      REPEAT(5)=REP23
      RETURN
      END
CDOC END_SUBROUTINE SPDEC3
