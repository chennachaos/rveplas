CDOC BEGIN_SUBROUTINE CONVRV
CDOC Checks for convergence of the equilibrium iterations of RVE
CDOC computational homogenisation problem
CDOC
CDOC This routine computes the global residual (out-of-balance) force
CDOC vector and its relative norm and checks whether the convergence
CDOC criterion has been satisfied.
CDOC If convergence has been attained the flag logical CONVRG returns as
CDOC .TRUE., otherwise it returns as .FALSE.
CDOC Also, if divergence is detected, the flag DIVERG returns as .TRUE.
CDOC (and .FALSE. otherwise).
CDOC
CDOC BEGIN_PARAMETERS
CDOC LOGICAL          CONVRG <  Convergence flag.
CDOC LOGICAL          DIVERG <  Divergence flag.
CDOC DOUBLE_PRECISION FBD    <  Boundary nodal reaction vector. 
CDOC INTEGER          IITER  >  Current equilibrium iteration number.
CDOC DOUBLE_PRECISION RESVEC <  Reduced vector of global internal force 
CDOC                            specific for the chosen kinematical
CDOC                            constraint
CDOC DOUBLE_PRECISION TOLER  >  Equilibrium convergence tolerance.
CDOC DOUBLE_PRECISION TFACT  >  Current total load factor.
CDOC END_PARAMETERS
CDOC
CHST M.F. Adziman, July 2013: Initially referring to CONVER.F (HYPLAS). 
CHST                          Added additional operations to the global  
CHST                          internal force vectors accommodating 
CHST                          micro-to-macro transitions (homogenisa-
CHST                          tion). Residual vector and residual norm 
CHST                          are specific for RVEPLAS.
CHST
      SUBROUTINE CONVRV
     1(   CONVRG     ,DIVERG     ,IITER      ,TOLER     ,TFACT       ,
     2    RESVEC     ,FBD)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C Hyplas database: Global parameters and common blocks
      INCLUDE '../MAXDIM.INC'
      INCLUDE '../MATERIAL.INC'
      INCLUDE '../ELEMENTS.INC'
      INCLUDE '../GLBDBASE.INC'
      INCLUDE '../RVE.INC'
C
      LOGICAL CONVRG,DIVERG
      SAVE RATOLD ,REMOLD
C
C Local arrays and variables
      DIMENSION 
     1   RESVEC(MTOTV)
      DIMENSION 
     1   FIN(MTOTV)    ,FPL(MTOTV)    ,FMI(MTOTV)     ,
     2   FCO(MTOTV)    ,FBD(MTOTV)    ,
     3   FFR(MTOTV)    ,FDP(MTOTV)    ,RTFD(MTOTV)   ,RAUX(MTOTV,4)
      DATA R0   ,R20   ,R100   ,R1000   /
     1     0.0D0,20.0D0,100.0D0,1000.0D0/
C***********************************************************************
C COMPUTE GLOBAL RESIDUAL (OUT-OF-BALANCE) FORCE VECTOR AND ITS RELATIVE
C NORM AND SET EQUILIBRIUM CONVERGENCE FLAG
C
C REFERENCE: Expressions (4.72) and (4.77)
C***********************************************************************
 1000 FORMAT(6X,I3,19X,G14.6,15X,G14.6)
 1010 FORMAT(6X,I3,11X,G14.6,6X,G14.6,6X,G14.6)
C Initialise relevant variables
C -----------------------------
      CONVRG=.FALSE.
      DIVERG=.FALSE.
      RESID=R0
      RETOT=R0
      REMAX=R0
C
C Evaluate global nodal internal and external forces
C --------------------------------------------------
      CALL RVZERO(STFOR,NTOTV)
      CALL RVZERO(TOFOR,NTOTV)
      DO 30 IELEM=1,NELEM
        IGRUP =IGRPID(IELEM)
        IELIDN=IELTID(IGRUP)
        NNODE =IELPRP(3,IELIDN)
        KEVAB=0
        DO 20 INODE=1,NNODE
          LOCNO=IABS(LNODS(IELEM,INODE))
          DO 10 IDOFN=1,NDOFN
            KEVAB=KEVAB+1
            NPOSI=MASTER((LOCNO-1)*NDOFN+IDOFN)
C current internal force
            STFOR(NPOSI)=STFOR(NPOSI)+ELOAD(KEVAB,IELEM)
C current external force
            TOFOR(NPOSI)=TOFOR(NPOSI)+TFACT*RLOAD(KEVAB,IELEM)
   10     CONTINUE
   20   CONTINUE
   30 CONTINUE
C Compute residual vector 
C -----------------------
C RESVEC needs to be in reduced order, like RHS in solrve
      RESVEC=R0
C Interior dofs
      DO I=IRV_DFSPPT(1),IRV_DFSPPT(2)
        IGL=IRV_DFSP(I)
        IRD=IRV_IFILT(IGL)
        IF(IRD/=0)RESVEC(IRD)=RESVEC(IRD)+STFOR(IGL)
      ENDDO
      IF(IRV_RVEOPT==2)THEN
C Periodic kinematic assumption...
C ...plus dofs
        DO I=IRV_DFSPPT(3),IRV_DFSPPT(4)
          IGL=IRV_DFSP(I)
          IRD=IRV_IFILT(IGL)
          IF(IRD/=0)RESVEC(IRD)=RESVEC(IRD)+STFOR(IGL)
        ENDDO
C ...minus dofs
        DO I=IRV_DFSPPT(5),IRV_DFSPPT(6)
          IGL=IRV_DFSP(I)
          IRD=IRV_IFILT(IGL)
          IF(IRD/=0)RESVEC(IRD)=RESVEC(IRD)+STFOR(IGL)
        ENDDO
C ...corner dofs
        DO I=IRV_DFSPPT(7),IRV_DFSPPT(8)
          IGL=IRV_DFSP(I)
          IRD=IRV_IFILT(IGL)
          IF(IRD/=0)RESVEC(IRD)=RESVEC(IRD)+STFOR(IGL)
        ENDDO
      ELSEIF(IRV_RVEOPT==3)THEN
C Uniform traction kinematic assumption...
C ...free dofs
        DO I=IRV_DFSPPT(9),IRV_DFSPPT(10)
          IGL=IRV_DFSP(I)
          IRD=IRV_IFILT(IGL)
          IF(IRD/=0)RESVEC(IRD)=RESVEC(IRD)+STFOR(IGL)
        ENDDO
C ...dependent dofs
        DO I=IRV_DFSPPT(11),IRV_DFSPPT(12)
          IGLDP=IRV_DFSP(I)
          IRDDP=IRV_IFILT(IGLDP)
          IF(IRDDP/=0)THEN
            DO J=IRV_DFSPPT(9),IRV_DFSPPT(10)
C DRV_BNDR's second index is ordered from 1 to NDEFDP
              JTEMP=J+1-IRV_DFSPPT(9)
              JGLFR=IRV_DFSP(J)
              JRDFR=IRV_IFILT(JGLFR)
              IF(JRDFR/=0)THEN
C Ff + R^T * Fd
                RESVEC(JRDFR)=RESVEC(JRDFR)
     1                       +DRV_BNDR(IRDDP,JTEMP)*STFOR(IGLDP)
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDIF
C
C Normalise residual norm
C -----------------------
C
C Number of dofs in reduced system for each assumption
      NDFIN=IRV_DFSPPT(2)-IRV_DFSPPT(1)+1
      IF(IRV_RVEOPT.EQ.1)THEN
        NDFSOL=NDFIN
      ELSEIF(IRV_RVEOPT.EQ.2)THEN
        NDFPL=IRV_DFSPPT(4)-IRV_DFSPPT(3)+1
        NDFSOL=NDFIN+NDFPL
      ELSEIF(IRV_RVEOPT.EQ.3)THEN
        NDFFR=IRV_DFSPPT(10)-IRV_DFSPPT(9)+1
        NDFSOL=NDFIN+NDFFR
      ENDIF
C
C Euclidean norm and maximum value of residual vector
      RVNORM=NORM2(RESVEC(1:NDFSOL))
      REMAX=MAXVAL(ABS(RESVEC(1:NDFSOL)))
C
C Evaluate norm of boundary dof reactive forces (used for normalising
C residual norm)
      FPNORM=R0
      IF((IRV_RVEOPT==1).OR.(IRV_RVEOPT==2))THEN
C Linear displacement and periodic assumptions...
C ...plus dofs
        DO I=IRV_DFSPPT(3),IRV_DFSPPT(4)
          IGL=IRV_DFSP(I)
          FPNORM=FPNORM+STFOR(IGL)*STFOR(IGL)
        ENDDO
C ...minus dofs
        DO I=IRV_DFSPPT(5),IRV_DFSPPT(6)
          IGL=IRV_DFSP(I)
          FPNORM=FPNORM+STFOR(IGL)*STFOR(IGL)
        ENDDO
C ...corner dofs
        DO I=IRV_DFSPPT(7),IRV_DFSPPT(8)
          IGL=IRV_DFSP(I)
          FPNORM=FPNORM+STFOR(IGL)*STFOR(IGL)
        ENDDO
      ELSEIF(IRV_RVEOPT==3)THEN
C Uniform traction assumption...
C ...free dofs
        DO I=IRV_DFSPPT(9),IRV_DFSPPT(10)
          IGL=IRV_DFSP(I)
          FPNORM=FPNORM+STFOR(IGL)*STFOR(IGL)
        ENDDO
C ...dependent dofs
        DO I=IRV_DFSPPT(11),IRV_DFSPPT(12)
          IGL=IRV_DFSP(I)
          FPNORM=FPNORM+STFOR(IGL)*STFOR(IGL)
        ENDDO
C ...prescribed dofs
        DO I=IRV_DFSPPT(13),IRV_DFSPPT(14)
          IGL=IRV_DFSP(I)
          FPNORM=FPNORM+STFOR(IGL)*STFOR(IGL)
        ENDDO
      ENDIF
      FPNORM=SQRT(FPNORM)
C
C Normalise residual norm
      IF(FPNORM.NE.R0)THEN
        RATIO=R100*RVNORM/FPNORM
      ELSE
        RATIO=R100*RVNORM
      ENDIF     
      IF(NALGO.GT.0)THEN
c        WRITE(16,1000)IITER,RATIO,REMAX
        WRITE(*,1000) IITER,RATIO,REMAX
      ELSE
c        WRITE(16,1010)IITER,RATIO,REMAX,TFACT
        WRITE(*,1010) IITER,RATIO,REMAX,TFACT
      ENDIF
C Set convergence/divergence flags
C --------------------------------
      IF(RATIO.LE.TOLER.OR.ABS(REMAX).LE.(TOLER/R1000))CONVRG=.TRUE.
      IF(IITER.NE.1.AND.(RATIO.GT.R20*RATOLD.OR.REMAX.GT.R20*REMOLD))
     1    DIVERG=.TRUE.
      RATOLD=RATIO
      REMOLD=REMAX
C Evaluate element residual forces before exit -> store it in ELOAD
C -----------------------------------------------------------------
      DO 100 IELEM=1,NELEM
        IGRUP=IGRPID(IELEM)
        IELIDN=IELTID(IGRUP)
        NEVAB=IELPRP(5,IELIDN)
C Current element residual (out-of-balance force) = current element
C external load - current element internal force
        DO 90 IEVAB=1,NEVAB
          ELOAD(IEVAB,IELEM)=TFACT*RLOAD(IEVAB,IELEM)-ELOAD(IEVAB,IELEM)
   90   CONTINUE
  100 CONTINUE
C
      RETURN
      END
CDOC END_SUBROUTINE CONVRV
