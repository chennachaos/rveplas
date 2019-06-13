CDOC BEGIN_SUBROUTINE SUDAMA
CDOC State update procedure for Lemaitre's ductile damage model.
CDOC
CDOC This routine uses the fully implicit elastic predictor/return
CDOC mapping algorithm as the state update procedure for Lemaitre's
CDOC ductile damage elasto-plastic model with general non-linear
CDOC (piece-wise linear) isotropic hardening under plane strain and
CDOC axisymmetric conditions.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DGAMA  <  Incremental plastic multiplier.
CDOC INTEGER          IPROPS >  Array of integer material properties.
CDOC C                          The number of points on the piece-wise
CDOC C                          linear hardening curve is the only
CDOC C                          element of this array used here.
CDOC C                          This array is set in routines
CDOC C                          INDATA and RDDAMA.
CDOC LOGICAL          LALGVA <  Array of logical algorithmic flags.
CDOC C                          For the present material model, this
CDOC C                          array contains the plastic yielding
CDOC C                          flag and the return algorithm failure
CDOC C                          flag. The plastic yielding flag is set
CDOC C                          to .TRUE. if plastic yielding has
CDOC C                          occurred and to .FALSE. if the step is
CDOC C                          elastic. The algorithm failure flag is
CDOC C                          set to .FALSE. if the state update
CDOC C                          algorithm has been successful and to
CDOC C                          .TRUE. if the return mapping algorithm
CDOC C                          has failed to converge.
CDOC INTEGER          NTYPE  >  Stress state type flag. The present
CDOC C                          routine is compatible only with NTYPE=2
CDOC C                          (plane strain) and NTYPE=3 (axisymmetric
CDOC C                          condition).
CDOC DOUBLE_PRECISION RPROPS >  Array of real material properties.
CDOC C                          This array contains the density
CDOC C                          (not used in this routine), the elastic
CDOC C                          properties: Young's modulus and
CDOC C                          Poisson's ratio, the two damage
CDOC C                          evolution constants and the pairs
CDOC C                          ``hardening variable strain-uniaxial
CDOC C                          yield stress'' defining the (user
CDOC C                          supplied) piece-wise linear hardening
CDOC C                          curve. This array is set in routine
CDOC C                          RDDAMA.
CDOC DOUBLE_PRECISION RSTAVA <> Array of real state variables other
CDOC C                          than the stress tensor components.
CDOC C                          Previous converged values on entry,
CDOC C                          updated values on exit.
CDOC C                          The state variables stored in
CDOC C                          this array are the (engineering)
CDOC C                          elastic strain components, the
CDOC C                          hardening variable and the
CDOC C                          damage internal variable.
CDOC DOUBLE_PRECISION STRAT  >  Array of elastic trial (engineering)
CDOC C                          strain components.
CDOC DOUBLE_PRECISION STRES  <  Array of updated stress tensor
CDOC C                          components.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, Jan 2000: Initial coding
CHST E.de Souza Neto, Jul 2008: Small changes in Newton-Raphson loop
CHST
      SUBROUTINE SUDAMA
     1(   DGAMA      ,IPROPS     ,LALGVA     ,NTYPE      ,RPROPS     ,
     2    RSTAVA     ,STRAT      ,STRES      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(IPHARD=6  ,MSTRE=4)
      LOGICAL IFPLAS, LALGVA(2), SUFAIL
      DIMENSION
     1    IPROPS(*)          ,RPROPS(*)          ,RSTAVA(MSTRE+2)    ,
     2    STRAT(MSTRE)       ,STRES(MSTRE)
      DIMENSION
     1    EET(MSTRE)
      DATA
     1    R0   ,RP5  ,R1   ,R2   ,R3   ,R6   ,SMALL ,TOL   / 
     2    0.0D0,0.5D0,1.0D0,2.0D0,3.0D0,6.0D0,1.D-20,1.D-08/
      DATA MXITER / 25 /
C***********************************************************************
C STATE UPDATE PROCEDURE FOR LEMAITRE'S DUCTILE DAMAGE MATERIAL MODEL
C WITH NON-LINEAR (PIECEWISE LINEAR) ISOTROPIC HARDENING:
C IMPLICIT ELASTIC PREDICTOR/RETURN MAPPING ALGORITHM.
C PLANE STRAIN AND AXISYMMETRIC IMPLEMENTATIONS.
C
C REFERENCE: Box 12.4
C***********************************************************************
C Stop program if neither plane strain nor axisymmetric state
      IF(NTYPE.NE.2.AND.NTYPE.NE.3)CALL ERRPRT('EI0051')
C Initialise some algorithmic and internal variables
      DGAMA=R0
      IFPLAS=.FALSE.
      SUFAIL=.FALSE.
C Retrieve hardening and damage internal variables
      HVARN=RSTAVA(MSTRE+1)
      DAMAGN=RSTAVA(MSTRE+2)
C... integrity
      OMEGAN=R1-DAMAGN
C Retrieve some material properties
      YOUNG=RPROPS(2)
      POISS=RPROPS(3)
      DAMEXP=RPROPS(4)
      DAMDEN=RPROPS(5)
      NHARD=IPROPS(3)
C Shear and bulk moduli and other necessary constants
      GMODU=YOUNG/(R2*(R1+POISS))
      BULK=YOUNG/(R3*(R1-R2*POISS))
      R2BULK=R2*BULK
      R2G=R2*GMODU
      R3G=R3*GMODU
      R6G=R6*GMODU
C Elastic predictor: Compute elastic trial state
C ==============================================
C Volumetric strain and (undamaged) pressure stress
      EEV=STRAT(1)+STRAT(2)+STRAT(4)
      PTILDE=BULK*EEV
C Elastic trial deviatoric strain
      EEVD3=EEV/R3
      EET(1)=STRAT(1)-EEVD3
      EET(2)=STRAT(2)-EEVD3
      EET(4)=STRAT(4)-EEVD3
C Convert engineering shear component into physical component
      EET(3)=STRAT(3)/R2
C Compute trial (undamaged) von Mises effective stress and uniaxial
C yield stress
      VARJ2T=R2G*R2G*(EET(3)*EET(3)+RP5*(EET(1)*EET(1)+
     1                     EET(2)*EET(2)+EET(4)*EET(4)))
      QTILTR=SQRT(R3*VARJ2T)
      SIGMAY=PLFUN(HVARN,NHARD,RPROPS(IPHARD))
C Check for plastic admissibility
C ===============================
      PHI=QTILTR-SIGMAY
      IF(PHI/SIGMAY.GT.TOL)THEN
C Plastic step: Apply return mapping - use Newton-Raphson algorithm
C               to solve the return mapping equation for DGAMA
C =================================================================
C Reset plastic flag
        IFPLAS=.TRUE.
C Initial guess for DGAMA: Use perfectly plastic solution with frozen
C yield surface at the beginning of the load increment
        DGAMA=OMEGAN*PHI/R3G
C Initialise hardening variable
        HVAR=HVARN+DGAMA
C Start N-R iterations
C --------------------
        PTILD2=PTILDE**2
        DO 10 NRITER=1,MXITER
C yield stress
          SIGMAY=PLFUN(HVAR,NHARD,RPROPS(IPHARD))
C integrity
          AUX2=R3G/(QTILTR-SIGMAY)
          OMEGA=AUX2*DGAMA
C stress triaxiality and damage energy release rate
          Y=-SIGMAY**2/R6G-PTILD2/R2BULK
C Compute residual function
          RES=OMEGA-OMEGAN+(-Y/DAMDEN)**DAMEXP/AUX2
C Check for convergence
          IF(ABS(RES).LE.TOL)THEN
C... update hardening and damage variables
            IF(OMEGA.LT.SMALL)THEN
C... check if converged damage variable is acceptable
              SUFAIL=.TRUE.
              CALL ERRPRT('WE0019')
              GOTO 999
            ENDIF
            DAMAGE=R1-OMEGA
            RSTAVA(MSTRE+1)=HVAR
            RSTAVA(MSTRE+2)=DAMAGE
C... update stress components
            P=OMEGA*PTILDE
            Q=SIGMAY*OMEGA
            FACTOR=R2G*Q/QTILTR
            STRES(1)=FACTOR*EET(1)+P
            STRES(2)=FACTOR*EET(2)+P
            STRES(3)=FACTOR*EET(3)
            STRES(4)=FACTOR*EET(4)+P
C... compute and store converged elastic (engineering) strain components
            FACTOR=R1-R3G*DGAMA/(OMEGA*QTILTR)
            RSTAVA(1)=FACTOR*EET(1)+EEVD3
            RSTAVA(2)=FACTOR*EET(2)+EEVD3
            RSTAVA(3)=FACTOR*EET(3)*R2
            RSTAVA(4)=FACTOR*EET(4)+EEVD3
            GOTO 999
          ENDIF
C Compute derivative of residual function
C... slope of hardening function
          HSLOPE=DPLFUN(HVAR,NHARD,RPROPS(IPHARD))
C... derivative of -Y
          DY=-HSLOPE*SIGMAY/R3G
C... residual derivative
          AUX=-Y/DAMDEN
          DRES=AUX2+AUX2*DGAMA*HSLOPE/(QTILTR-SIGMAY)-
     1         HSLOPE/R3G*AUX**DAMEXP-
     2         DAMEXP*DY/(AUX2*DAMDEN)*AUX**(DAMEXP-R1)
C Apply N-R correction to DGAMA
          DDGAMA=-RES/DRES
          DGAMA=DGAMA+DDGAMA
C... update hardening variable
          HVAR=HVARN+DGAMA
   10   CONTINUE
C N-R loop failed to converge: Reset failure flag and issue warning
C                              message before exiting
        SUFAIL=.TRUE.
        CALL ERRPRT('WE0018')
      ELSE
C Elastic step: Update stress using damaged elastic law
C =====================================================
        FACTOR=R2G*OMEGAN
        P=OMEGAN*PTILDE
        STRES(1)=FACTOR*EET(1)+P
        STRES(2)=FACTOR*EET(2)+P
        STRES(3)=FACTOR*EET(3)
        STRES(4)=FACTOR*EET(4)+P
C elastic engineering strain
        RSTAVA(1)=STRAT(1)
        RSTAVA(2)=STRAT(2)
        RSTAVA(3)=STRAT(3)
        RSTAVA(4)=STRAT(4)
      ENDIF
  999 CONTINUE
C Update some algorithmic variables before exit
      LALGVA(1)=IFPLAS
      LALGVA(2)=SUFAIL
      RETURN
      END
CDOC END_SUBROUTINE SUDAMA
