CDOC BEGIN_SUBROUTINE LENGTH
CDOC Adjusts step length for the Arc-Length Method.
CDOC
CDOC This subroutine adjusts the step length for the Arc-Length Method
CDOC according to the prescribed desired number of iterations for
CDOC equilibrium convergence and the actual number of iterations for
CDOC convergence in the previous load step.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DLEN   <  Calculated step length.
CDOC DOUBLE_PRECISION DLENM  >  Prescribed maximum permissible step
CDOC C                          length.
CDOC INTEGER          ITACT  >  Actual number of iterations required for
CDOC C                          convergence in the previous load step.
CDOC INTEGER          ITDES  >  Desired number of iterations for
CDOC C                          convergence in the iterative solution of
CDOC C                          the non-linear finite element
CDOC C                          equilibrium equations.
CDOC END_PARAMETERS
CDOC
      SUBROUTINE LENGTH(DLENG ,DLENM ,ITACT ,ITDES )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C***********************************************************************
C ADJUSTS STEP LENGTH ACCORDING TO THE DESIRED NUMBER OF ITERATIONS AND
C THE NUMBER OF ITERATIONS REQUIRED FOR CONVERGENCE IN THE PREVIOUS
C LOAD STEP (USED FOR ARC-LENGTH METHOD ONLY)
C
C REFERENCE: Expression (5.3)
C***********************************************************************
      DLENG=DLENG*DBLE(ITDES)/DBLE(ITACT)
      DLENG=MIN(DLENG,DLENM)
      RETURN
      END
CDOC END_SUBROUTINE LENGTH
