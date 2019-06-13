CDOC BEGIN_SUBROUTINE SWMTSC
      SUBROUTINE SWMTSC
     1(   MODE       ,LALGVC     ,LALGVL     ,RALGVC     ,RSTAVC     ,
     2    RSTAVL     ,STRESC     ,STRESL     )
      
CC      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      IMPLICIT NONE
      INTEGER, PARAMETER :: NSTRE=6, NRSTAV=52, NLALGV=6, NRALGV=1
C Arguments
      
CC      LOGICAL
CC     1    LALGVC             ,LALGVL
CC      DIMENSION
CC     1    LALGVC(*)          ,LALGVL(*)          ,RALGVC(*)          ,
CC     2    RSTAVC(*)          ,RSTAVL(*)          ,STRESC(*)          ,
CC     3    STRESL(*)
      INTEGER, INTENT(IN) :: MODE
      LOGICAL, DIMENSION(NLALGV), INTENT(INOUT) :: LALGVC, LALGVL
      DOUBLE PRECISION, DIMENSION(NRSTAV), INTENT(INOUT) :: RSTAVC,
     1                                                      RSTAVL
      DOUBLE PRECISION, DIMENSION(NRALGV), INTENT(INOUT) :: RALGVC
      DOUBLE PRECISION, DIMENSION(NSTRE), INTENT(INOUT) :: STRESC,
     1                                                     STRESL
      
CC      DATA R0   ,R1   /
CC     1     0.0D0,1.0D0/
      
C Local parameters
      DOUBLE PRECISION, PARAMETER :: R0=0.0D0, R1=1.0D0
C***********************************************************************
C INITIALISE/SWITCH DATA FOR THE MARTENSITIC TRANSFORMAITON MATERIAL 
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
      
CC      NSTRE=6
CC      NRSTAV=52
CC      NLALGV=6
CC      NRALGV=1
      
C
C Initialisation mode
C ===================
      IF(MODE==0)THEN
CC        CALL RVZERO(STRESC,NSTRE)
CC        CALL RVZERO(RALGVC,NRALGV)
CC        DO 10 I=1,NLALGV
CC          LALGVC(I)=.FALSE.
CC   10   CONTINUE
        STRESC=R0
        RALGVC=R0
        LALGVC=.FALSE.
CC        RSTAVC(1:NRSTAV)=R0
        RSTAVC=R0
C RSTAVA 1-9: elastic deformation gradient
        RSTAVC([1,5,9])=R1 !Components (1,1), (2,2) and (3,3)
C RSTAVA 10-18: austenite plastic deformation gradient
        RSTAVC([10,14,18])=R1 !Components (1,1), (2,2) and (3,3)
C RSTAVA 19: accumulated plastic slip
C RSTAVA 20-43: volume fraction of transformation systems
C RSTAVA 44: total martensite volume fraction
C Switching modes
C ===============
      ELSE
        IF(MODE==1)THEN
CC          DO 20 I=1,NSTRE
CC            STRESL(I)=STRESC(I)
CC   20     CONTINUE
CC          DO 30 I=1,NRSTAV
CC            RSTAVL(I)=RSTAVC(I)
CC   30     CONTINUE
CC          DO 40 I=1,NLALGV
CC            LALGVL(I)=LALGVC(I)
CC   40     CONTINUE
CC          CALL RVZERO(RALGVC,NRALGV)
          STRESL=STRESC
          RSTAVL=RSTAVC
          LALGVL=LALGVC
          
C Zero plastic multipliers before starting a new increment
          RALGVC=R0
          
        ELSEIF(MODE==2.OR.MODE==3)THEN
CC          DO 50 I=1,NSTRE
CC            STRESC(I)=STRESL(I)
CC   50     CONTINUE
CC          DO 60 I=1,NRSTAV
CC            RSTAVC(I)=RSTAVL(I)
CC   60     CONTINUE
CC          DO 70 I=1,NLALGV
CC            LALGVC(I)=LALGVL(I)
CC   70     CONTINUE
          STRESC=STRESL
          RSTAVC=RSTAVL
          LALGVC=LALGVL
          
          
          
          
C Set once a transformation is identified, it is irreversible 
C including in the increment step (!!under review)  
C ***************************************	            
C8          IFTRAN=0
C8          IF(RSTAVC(52).EQ.R1)IFTRAN=1         
C ***************************************
C Set a flag identifying whether or not a transformation has been  
C initiated during the previous converged solution         
C8          IF(RSTAVL(39).NE.R0)THEN
C8            RSTAVC(40)=R1
C8          ENDIF
C Set a flag identifying the state of transformation in the previous  
C converged solution          
C8          IF(LALGVL(1))THEN
C8            RSTAVC(51)=R1
C8          ENDIF
C Set once a transformation is identified, it is irreversible 
C including in the increment step (!!under review)	    
C ****************************
C8          RSTAVC(52)=R0 	    
C8          IF(IFTRAN.EQ.1)THEN
C8            RSTAVC(52)=R1
C8          ELSE
C8            RSTAVC(52)=R0
C8          ENDIF
C ****************************
          

C Zero plastic multipliers before starting a new increment
          IF(MODE==3)THEN
CC            CALL RVZERO(RALGVC,NRALGV)
            RALGVC=R0
          ENDIF
          
        ENDIF
      ENDIF
      RETURN
      END
CDOC END_SUBROUTINE SWMTSC