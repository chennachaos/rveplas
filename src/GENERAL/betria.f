CDOC BEGIN_SUBROUTINE BETRIA
CDOC Computes the left Cauchy-Green strain tensor
CDOC
CDOC Given the previous elastic logarithmic strains and the incremental
CDOC deformation gradient between the previous and current
CDOC configuration, this routine computes the current elastic trial left
CDOC Cauchy-Green strain tensor.
CDOC This routine contains the plane strain, plane stress and
CDOC axisymmetric implementations.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION BETRL  <  Array of components of the current
CDOC C                          elastic trial (or total) left
CDOC C                          Cauchy-Green strain tensor.
CDOC DOUBLE_PRECISION EEN    >  Array of previous elastic engineering
CDOC C                          logarithmic strains.
CDOC DOUBLE_PRECISION FINCR  >  Incremental deformation gradient between
CDOC C                          the previous and current configuration.
CDOC INTEGER          NTYPE  >  Stress state type. Present routine is
CDOC C                          compatible with NTYPE=1 (plane stress),
CDOC C                          NTYPE=2 (plane strain) and NTYPE=3
CDOC C                          (axisymmetric condition).
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, August 1996: Initial coding
CHST E.de Souza Neto, June   2003: List of arguments changed. Now the
CHST                               previous elastic logarithmic strain
CHST                               replaces the previous elastic left
CHST                               Cauchy-Green tensor.
CHST D. de Bortoli  , March  2015: 3-D case added
CHST
      SUBROUTINE BETRIA
     1(   BETRL      ,EEN        ,FINCR      ,NTYPE      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION
     1    BETRL(*)           ,EEN(*)             ,FINCR(3,3)
      DIMENSION
     1    AUXM(3,3)          ,BEN(6)              ,BENMTX(3,3)        ,
     2    BETRLM(3,3)
C***********************************************************************
C COMPUTES THE "ELASTIC TRIAL" LEFT CAUCHY-GREEN STRAIN TENSOR FOR
C HYPERELASTIC-BASED LARGE STRAIN ELASTO-PLASTIC MODELS ACCORDING TO
C THE FORMULA:
C
C                     e trial             e   T
C                    B        :=   F     B   F
C                     n+1           incr  n   incr
C         e
C WHERE  B   IS OBTAINED AS:
C         n
C                     e              e
C                    B   :=  exp[ 2 E  ]
C                     n              n
C        e
C WITH  E   DENOTING THE ELASTIC LOGARITHMIC STRAIN AT  t .
C        n                                               n
C
C REFERENCE: Box 14.3, item (ii)
C***********************************************************************
      IF((NTYPE.EQ.1).OR.(NTYPE.EQ.2).OR.(NTYPE.EQ.3))THEN
        NDIM=2
        NSTRA=4
      ELSEIF(NTYPE.EQ.4)THEN
        NDIM=3
        NSTRA=6
      ELSE
        CALL ERRPRT('EI0066')
      ENDIF
C
      BEN(1:NSTRA)=EEN(1:NSTRA)
C Convert engineering elastic logarithmic strain components into the
C corresponding elastic left Cauchy-Green strain tensor components
      CALL SETBE(BEN,NTYPE)
C Convert left Cauchy-Green strain tensor from vector array to matrix
C form
      CALL ATSYM(BEN, BENMTX, NDIM, .FALSE.)
C
C Compute elastic trial left Cauchy-Green tensor (for the plane stress/
C strain and axisymmetric cases, compute the in-plane components only)
C
      CALL RVZERO(AUXM,9)
      DO 30 I=1,NDIM
        DO 20 J=1,NDIM
          DO 10 K=1,NDIM
            AUXM(I,J)=AUXM(I,J)+FINCR(I,K)*BENMTX(K,J)
   10     CONTINUE
   20   CONTINUE
   30 CONTINUE
      CALL RVZERO(BETRLM,9)
      DO 60 I=1,NDIM
        DO 50 J=1,NDIM
          DO 40 K=1,NDIM
            BETRLM(I,J)=BETRLM(I,J)+AUXM(I,K)*FINCR(J,K)
   40     CONTINUE
   50   CONTINUE
   60 CONTINUE
C
C        e trial
C Store B        in array form
C        n+1 
C
      CALL SYMTA(BETRLM, BETRL, NDIM, .FALSE.)
C set out-of-plane component for plane strain and axisymmetric analyses
      IF(NTYPE.EQ.2)THEN
        BETRL(4)=BEN(4)
      ELSEIF(NTYPE.EQ.3)THEN
        BETRL(4)=BEN(4)*FINCR(3,3)**2
      ENDIF
C
      RETURN
      END
CDOC END_SUBROUTINE BETRIA
