CDOC BEGIN_SUBROUTINE SWOGD
CDOC Initialise/switch state variables for the Ogden hyperelastic model
CDOC
CDOC This initialises and switches state variables (between current and
CDOC previous values) for the Ogden hyperelastic material model.
CDOC
CDOC BEGIN_PARAMETERS
CDOC INTEGER          MODE   >  Initialisation/Switching mode. See
CDOC C                          source code comments for details.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC DOUBLE_PRECISION RSTAVC <> Array of real state variables at Gauss
CDOC C                          point. Current values.
CDOC DOUBLE_PRECISION RSTAVL <> Array of real state variables at Gauss
CDOC C                          point. Last converged (equilibrium)
CDOC C                          values.
CDOC DOUBLE_PRECISION STRESC <> Array of stress (Cauchy in large strain)
CDOC C                          components. Current values.
CDOC DOUBLE_PRECISION STRESL <> Array of stress (Cauchy in large strain)
CDOC C                          components. Last converged (equilibrium)
CDOC C                          values.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, July      1999: Initial coding
CHST
      SUBROUTINE SWOGD
     1(   MODE       ,NTYPE      ,RSTAVC     ,RSTAVL     ,STRESC     ,
     2    STRESL     )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C Arguments
      DIMENSION
     1    RSTAVC(*)          ,RSTAVL(*)          ,STRESC(*)          ,
     2    STRESL(*)
C Local numerical constants
      DATA R0   ,R1   /
     1     0.0D0,1.0D0/
C***********************************************************************
C INITIALISE/SWITCH DATA FOR THE OGDEN MATERIAL MODEL
C
C    MODE=0:   Initialises the relevant data.
C
C    MODE=1:   Assigns current values of the state variables to
C              converged solution (when the current iteration
C              satisfies the convergence criterion).
C
C    MODE=2:   Assigns the last converged solution to current state
C              variables values (when a new iteration is required by
C              the iterative process).
C
C    MODE=3:   Assigns the last converged solution to current state
C              variables values (when increment cutting is required).
C***********************************************************************
C
      IF(NTYPE.EQ.1.OR.NTYPE.EQ.2.OR.NTYPE.EQ.3)THEN
        NSTRE=4
        NRSTAV=4
      ENDIF
C
      IF(MODE.EQ.0)THEN
C Initialisation mode
C ===================
C Zero stress component array
        CALL RVZERO(STRESC,NSTRE)
C RSTAVA stores the left Cauchy-Green strain tensor components.
C Initialised as identity
        IF(NTYPE.EQ.1.OR.NTYPE.EQ.2.OR.NTYPE.EQ.3)THEN
          RSTAVC(1)=R1
          RSTAVC(2)=R1
          RSTAVC(3)=R0
          RSTAVC(4)=R1
        ENDIF
C Switching modes
C ===============
      ELSEIF(MODE.EQ.1)THEN
        DO 10 I=1,NSTRE
          STRESL(I)=STRESC(I)
   10   CONTINUE
        DO 20 I=1,NRSTAV
          RSTAVL(I)=RSTAVC(I)
   20   CONTINUE
      ELSEIF(MODE.EQ.2.OR.MODE.EQ.3)THEN
        DO 30 I=1,NSTRE
          STRESC(I)=STRESL(I)
   30   CONTINUE
        DO 40 I=1,NRSTAV
          RSTAVC(I)=RSTAVL(I)
   40   CONTINUE
      ENDIF
      RETURN
      END
CDOC END_SUBROUTINE SWOGD
