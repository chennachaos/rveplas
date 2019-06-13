CDOC BEGIN_SUBROUTINE SUVMPS
CDOC State update procedure for the von Mises model in plane stress.
CDOC
CDOC This routine uses the fully implicit elastic predictor/return
CDOC mapping algorithm as the state update procedure for the von Mises
CDOC elasto-plastic material model with general non-linear (piece-wise
CDOC linear) isotropic hardening under plane stress condition.
CDOC The algorithm used here is derived from the plane stress-projected
CDOC von Mises model, in which only in-plane stress/strain components
CDOC take part in the evolution equations.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DGAMA  <  Incremental plastic multiplier.
CDOC INTEGER          IPROPS >  Array of integer material properties.
CDOC C                          The number of points on the piece-wise
CDOC C                          linear hardening curve is the only
CDOC C                          element of this array used here.
CDOC C                          This array is set in routines INDATA
CDOC C                          and RDVM.
CDOC LOGICAL          LALGVA <  Array of logical algorithmic flags.
CDOC C                          For the present material model, this
CDOC C                          array contains the plastic yielding
CDOC C                          flag and the return algorithm failure
CDOC C                          flag. The plastic yielding flag is set
CDOC C                          to .TRUE. if plastic yielding
CDOC C                          has occurred and to .FALSE. if
CDOC C                          the step is elastic. The algorithm
CDOC C                          failure flag is set to .FALSE. if
CDOC C                          the state update algorithm has been
CDOC C                          successful and to .TRUE. if the
CDOC C                          return mapping algorithm has failed
CDOC C                          to converge.
CDOC INTEGER          NTYPE  >  Stress state type flag. The present
CDOC C                          routine is compatible only with NTYPE=1
CDOC C                          (plane stress).
CDOC DOUBLE_PRECISION RPROPS >  Array of real material properties.
CDOC C                          This array contains the density
CDOC C                          (not used in this routine), the elastic
CDOC C                          properties: Young's modulus and
CDOC C                          Poisson's ratio, and the pairs
CDOC C                          ``accumulated plastic strain-uniaxial
CDOC C                          yield stress'' defining the (user
CDOC C                          supplied) piece-wise linear hardening
CDOC C                          curve. This array is set in routine
CDOC C                          RDVM.
CDOC DOUBLE_PRECISION RSTAVA <> Array of real state variables other
CDOC C                          than the stress tensor components.
CDOC C                          Previous converged values on entry,
CDOC C                          updated values on exit.
CDOC C                          The state variables stored in
CDOC C                          this array are the (engineering)
CDOC C                          elastic strain components and the
CDOC C                          accumulated plastic strain.
CDOC DOUBLE_PRECISION STRAT  >  Array of in-plane elastic trial
CDOC                            (engineering) strain components.
CDOC DOUBLE_PRECISION STRES  <  Array of updated stress tensor
CDOC C                          components.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, Aug 1998: Initial coding
CHST
      SUBROUTINE SUVMPS
     1(   DGAMA      ,IPROPS     ,LALGVA     ,NTYPE      ,RPROPS     ,
     2    RSTAVA     ,STRAT      ,STRES      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER( IPHARD=4  ,MSTRE=4  ,NSTRE=3 )
      LOGICAL IFPLAS, LALGVA(2), SUFAIL
      DIMENSION
     1    IPROPS(*)          ,RPROPS(*)          ,RSTAVA(MSTRE+1)    ,
     2    STRAT(MSTRE)       ,STRES(MSTRE)
      DIMENSION
     1    EET(MSTRE)         ,STREST(NSTRE)
      DATA
     1    R0   ,RP5  ,R1   ,R2   ,R3   ,R4   ,R6   ,TOL   / 
     2    0.0D0,0.5D0,1.0D0,2.0D0,3.0D0,4.0D0,6.0D0,1.D-08/
      DATA MXITER / 50 /
C***********************************************************************
C STATE UPDATE PROCEDURE FOR THE VON MISES ELASTO-PLASTIC MODEL WITH
C NON-LINEAR (PIECEWISE LINEAR) ISOTROPIC HARDENING IN PLANE STRESS:
C IMPLICIT PLANE STRESS-PROJECTED ELASTIC PREDICTOR/RETURN MAPPING
C ALGORITHM (BOXES 9.4-5).
C
C REFERENCE: Section 9.4.3
C            Boxes 9.4-5
C***********************************************************************
C Stop program if not plane stress
      IF(NTYPE.NE.1)CALL ERRPRT('EI0031')
C Initialise some algorithmic and internal variables
      DGAMA=R0
      IFPLAS=.FALSE.
      SUFAIL=.FALSE.
C...set previously (equilibrium) converged accumulated plastic strain
      EPBARN=RSTAVA(MSTRE+1)
C Set some material properties
      YOUNG=RPROPS(2)
      POISS=RPROPS(3)
      NHARD=IPROPS(3)
C Shear and bulk moduli and other necessary constants
      GMODU=YOUNG/(R2*(R1+POISS))
      BULK=YOUNG/(R3*(R1-R2*POISS))
      R2G=R2*GMODU
      R4G=R4*GMODU
      R1D3=R1/R3
      R1D6=R1/R6
      R2D3=R2*R1D3
      SQR2D3=SQRT(R2D3)
      R4GD3=R4G*R1D3
C Elastic predictor: Compute elastic trial state
C ----------------------------------------------
C Volumetric strain
      FACTOR=R2G/(BULK+R4GD3)
      EEV=(STRAT(1)+STRAT(2))*FACTOR
C Elastic trial deviatoric strain
      EEVD3=EEV/R3
      EET(1)=STRAT(1)-EEVD3
      EET(2)=STRAT(2)-EEVD3
C Convert engineering shear component into physical component
      EET(3)=STRAT(3)*RP5
C Elastic trial stress components
      PT=BULK*EEV
      STREST(1)=R2G*EET(1)+PT
      STREST(2)=R2G*EET(2)+PT
      STREST(3)=R2G*EET(3)
C Compute yield function value at trial state
      A1=(STREST(1)+STREST(2))*(STREST(1)+STREST(2))
      A2=(STREST(2)-STREST(1))*(STREST(2)-STREST(1))
      A3=STREST(3)*STREST(3)
      XI=R1D6*A1+RP5*A2+R2*A3
      SIGMAY=PLFUN(EPBARN,NHARD,RPROPS(IPHARD))
C...yield function
      PHI=RP5*XI-R1D3*SIGMAY*SIGMAY
C Check for plastic admissibility
C -------------------------------
      IF(PHI/SIGMAY.GT.TOL)THEN
C Plastic step: Apply return mapping - use Newton-Raphson algorithm
C               to solve the plane stress-projected return mapping
C               equation for the plastic multiplier (Box 9.5)
C -----------------------------------------------------------------
        IFPLAS=.TRUE.
        EPBAR=EPBARN
        SQRTXI=SQRT(XI)
        B1=R1
        B2=R1
        FMODU=YOUNG/(R3*(R1-POISS))
        DO 10 NRITER=1,MXITER
C Compute residual derivative
          HSLOPE=DPLFUN(EPBAR,NHARD,RPROPS(IPHARD))
          DXI=-A1*FMODU/(R3*B1*B1*B1)-R2G*(A2+R4*A3)/(B2*B2*B2)
          HBAR=R2*SIGMAY*HSLOPE*SQR2D3*(SQRTXI+DGAMA*DXI/(R2*SQRTXI))
          DPHI=RP5*DXI-R1D3*HBAR
C Compute Newton-Raphson increment and update equation variable DGAMA
          DGAMA=DGAMA-PHI/DPHI
C Compute new residual (yield function value)
          B1=R1+FMODU*DGAMA
          B2=R1+R2G*DGAMA
          XI=R1D6*A1/(B1*B1)+(RP5*A2+R2*A3)/(B2*B2)
          SQRTXI=SQRT(XI)
          EPBAR=EPBARN+DGAMA*SQR2D3*SQRTXI
          SIGMAY=PLFUN(EPBAR,NHARD,RPROPS(IPHARD))
          PHI=RP5*XI-R1D3*SIGMAY*SIGMAY
C Check for convergence
          RESNOR=ABS(PHI/SIGMAY)
          IF(RESNOR.LE.TOL)THEN
C update accumulated plastic strain
            RSTAVA(MSTRE+1)=EPBAR
C update stress components:   sigma := A sigma^trial
            ASTAR1=R3*(R1-POISS)/(R3*(R1-POISS)+YOUNG*DGAMA)
            ASTAR2=R1/(R1+R2G*DGAMA)
            A11=RP5*(ASTAR1+ASTAR2)
            A22=A11
            A12=RP5*(ASTAR1-ASTAR2)
            A21=A12
            A33=ASTAR2
            STRES(1)=A11*STREST(1)+A12*STREST(2)
            STRES(2)=A21*STREST(1)+A22*STREST(2)
            STRES(3)=A33*STREST(3)
C compute corresponding elastic (engineering) strain components
            FACTG=R1/R2G
            P=R1D3*(STRES(1)+STRES(2))
            EEV=P/BULK
            EEVD3=R1D3*EEV
            RSTAVA(1)=FACTG*(R2D3*STRES(1)-R1D3*STRES(2))+EEVD3
            RSTAVA(2)=FACTG*(R2D3*STRES(2)-R1D3*STRES(1))+EEVD3
            RSTAVA(3)=FACTG*STRES(3)*R2
            RSTAVA(4)=-POISS/(R1-POISS)*(RSTAVA(1)+RSTAVA(2))
            GOTO 999
          ENDIF
   10   CONTINUE
C reset failure flag and print warning message if N-R algorithm fails
        SUFAIL=.TRUE.
        CALL ERRPRT('WE0013')
      ELSE
C Elastic step: Update stress using linear elastic law
C ----------------------------------------------------
        STRES(1)=STREST(1)
        STRES(2)=STREST(2)
        STRES(3)=STREST(3)
C elastic engineering strain
        RSTAVA(1)=STRAT(1)
        RSTAVA(2)=STRAT(2)
        RSTAVA(3)=STRAT(3)
        RSTAVA(4)=-POISS/(R1-POISS)*(STRAT(1)+STRAT(2))
      ENDIF
  999 CONTINUE
C Update some algorithmic variables before exit
      LALGVA(1)=IFPLAS
      LALGVA(2)=SUFAIL
      RETURN
      END
CDOC END_SUBROUTINE SUVMPS
