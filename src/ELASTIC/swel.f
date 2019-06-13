CDOC BEGIN_SUBROUTINE SWEL
CDOC Initialise/switch state variables for the elastic material model
CDOC
CDOC This initialises and switches state variables (between current and
CDOC previous values) for the elastic material model (Hencky material
CDOC in large strain analysis).
CDOC
CDOC BEGIN_PARAMETERS
CDOC INTEGER          MODE   >  Initialisation/Switching mode.
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
CHST E.de Souza Neto, October   2008: Unused arguments removed
CHST D. de Bortoli  , March     2015: 3-D case added (NTYPE=4)
CHST
      SUBROUTINE SWEL
     1(   MODE       ,NTYPE      ,RSTAVC     ,RSTAVL     ,STRESC     ,
     2    STRESL     )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C Arguments
      DIMENSION
     1    RSTAVC(*)          ,RSTAVL(*)          ,STRESC(*)          ,
     2    STRESL(*)
C***********************************************************************
C INITIALISE/SWITCH DATA FOR THE ELASTIC MATERIAL MODEL
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
      IF((NTYPE.EQ.1).OR.(NTYPE.EQ.2).OR.(NTYPE.EQ.3))THEN
        NSTRE=4
      ELSEIF(NTYPE.EQ.4)THEN
        NSTRE=6
      ENDIF
C
      IF(MODE.EQ.0)THEN
C Initialisation mode
C ===================
        CALL RVZERO(STRESC,NSTRE)
C RSTAVA stores the infinitesimal engineering strain tensor components
C in small strains and the logarithmic eng. strains in large strain
C analysis
        CALL RVZERO(RSTAVC,NSTRE)
      ELSE
C Switching mode
C ==============
        IF(MODE.EQ.1)THEN
          DO 10 ISTRE=1,NSTRE
            STRESL(ISTRE)=STRESC(ISTRE)
            RSTAVL(ISTRE)=RSTAVC(ISTRE)
   10     CONTINUE
        ELSEIF(MODE.EQ.2.OR.MODE.EQ.3)THEN
          DO 20 ISTRE=1,NSTRE
            STRESC(ISTRE)=STRESL(ISTRE)
            RSTAVC(ISTRE)=RSTAVL(ISTRE)
   20     CONTINUE
        ENDIF
      ENDIF
      RETURN
      END
CDOC END_SUBROUTINE SWEL
