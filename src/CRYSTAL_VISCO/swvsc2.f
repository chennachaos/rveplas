CDOC BEGIN_SUBROUTINE SWVPCR
CDOC Initialise/switch state variables for the viscoplastic single
CDOC crystal model (plane strain and 3D implementations).
CDOC
CDOC BEGIN_PARAMETERS
CDOC INTEGER          MODE   >  Initialisation/Switching mode.
CDOC DOUBLE_PRECISION LALGVC <> Array of logical algorithmic variables
CDOC C                          at Gauss point. Current values.
CDOC DOUBLE_PRECISION LALGVL <> Array of logical algorithmic variables
CDOC C                          at Gauss point. Last converged
CDOC C                          (equilibrium) values.
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
CHST D. de Bortoli,   May       2016: Updated real state variable
CHST                                  numbers and initialisation for 3D
CHST                                  implementation;
CHST                                  removed unused arguments
CHST
C1      SUBROUTINE SWVPCR
      SUBROUTINE SWVSC2
     1(  MODE, LALGVC, LALGVL, RSTAVC, RSTAVL, STRESC, STRESL  )
      IMPLICIT NONE
C
      INTEGER, PARAMETER :: NSTRE=6, NRSTAV=19, NLALGV=2, NRALGV=0
C Arguments
      INTEGER, INTENT(IN) :: MODE
      LOGICAL, DIMENSION(NLALGV), INTENT(INOUT) :: LALGVC, LALGVL
      DOUBLE PRECISION, DIMENSION(NRSTAV), INTENT(INOUT) :: RSTAVC,
     1                                                      RSTAVL
      DOUBLE PRECISION, DIMENSION(NSTRE), INTENT(INOUT) :: STRESC,
     1                                                     STRESL
C Local parameters
      DOUBLE PRECISION, PARAMETER :: R0=0.0D0, R1=1.0D0
C***********************************************************************
C INITIALISE/SWITCH DATA FOR THE VISCO-PLASTIC SINGLE CRYSTAL MATERIAL
C MODEL (PLANE STRAIN AND 3D IMPLEMENTATIONS)
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
C Initialisation mode
C ===================
      IF(MODE==0)THEN
        STRESC=R0
        LALGVC=.FALSE.
        RSTAVC=R0
C RSTAVA 1-9: elastic deformation gradient
        RSTAVC([1,5,9])=R1 ! set components (1,1), (2,2) and (3,3) to 1
C RSTAVA 10-18: plastic deformation gradient
        RSTAVC([10,14,18])=R1
C RSTAVA 19: accumulated plastic slip (hardening variable)
C
C Switching modes
C ===============
      ELSEIF(MODE==1)THEN
        STRESL=STRESC
        LALGVL=LALGVC
        RSTAVL=RSTAVC
      ELSEIF((MODE==2).OR.(MODE==3))THEN
        STRESC=STRESL
        LALGVC=LALGVL
        RSTAVC=RSTAVL
      ENDIF
C
      RETURN
      END
CDOC END_SUBROUTINE SWVPCR
