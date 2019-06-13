CDOC BEGIN_SUBROUTINE ORDAMA
CDOC Output results Lemaitre's ductile damage elasto-plastic model
CDOC
CDOC This routine writes to the results file the internal and
CDOC algorithmic variables of Lemaitre's ductile damage elasto-plastic
CDOC material with non-linear (piece-wise linear) isotropic hardening.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DGAMA  >  Incremental plastic multiplier. Computed
CDOC C                          in routine SUDAMA.
CDOC INTEGER          NOUTF  >  Results file unit identifier.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC DOUBLE_PRECISION RSTAVA >  Array of real state variables other than
CDOC C                          the stress tensor components.
CDOC DOUBLE_PRECISION STRES  >  Array of stress tensor components.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, January 2000: Initial coding
CHST
      SUBROUTINE ORDAMA
     1(   DGAMA      ,NOUTF      ,NTYPE      ,RSTAVA     ,STRES      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(IPHARD=6  ,MSTRE=4)
      DIMENSION  RSTAVA(MSTRE+2), STRES(*)
      DATA   R2   ,R3    / 2.0D0,3.0D0 /
C***********************************************************************
C OUTPUT RESULTS (INTERNAL AND ALGORITHMIC VARIABLES) FOR LEMAITRE'S
C DUCTILE DAMAGE ELASTO-PLASTIC MODEL WITH NON-LINEAR ISOTROPIC
C HARDENING
C***********************************************************************
 1000 FORMAT(' S-eff = ',G12.4,' R     = ',G12.4,' dgama = ',G12.4)
 2000 FORMAT(' Damage= ',G12.4)
C
C Retrieve current values of hardening variable and damage
      HVAR=RSTAVA(MSTRE+1)
      DAMAGE=RSTAVA(MSTRE+2)
      IF(NTYPE.EQ.1)THEN
C Plane stress
        P=(STRES(1)+STRES(2))/R3
        EFFST=SQRT(R3/R2*((STRES(1)-P)**2+(STRES(2)-P)**2+
     1                     R2*STRES(3)**2+P**2))
      ELSEIF(NTYPE.EQ.2.OR.NTYPE.EQ.3)THEN
C Plane strain and axisymmetric
        P=(STRES(1)+STRES(2)+STRES(4))/R3
        EFFST=SQRT(R3/R2*((STRES(1)-P)**2+(STRES(2)-P)**2+
     1                     R2*STRES(3)**2+(STRES(4)-P)**2))
      ENDIF
C Write to output file
      WRITE(NOUTF,1000)EFFST,HVAR,DGAMA
      WRITE(NOUTF,2000)DAMAGE
      RETURN
      END
CDOC END_SUBROUTINE ORDAMA
