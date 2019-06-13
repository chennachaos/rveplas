CDOC BEGIN_SUBROUTINE ORMC
CDOC Output results for the Mohr-Coulomb elasto-plastic material model
CDOC
CDOC This routine writes to the results file the internal and
CDOC algorithmic variables of the Mohr-Coulomb elasto-plastic material
CDOC with non-linear isotropic hardening. The results printed here
CDOC are obtained by the return mapping algorithm implemented in routine
CDOC SUMC.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DGAM   >  Array of incremental plastic
CDOC C                          multipliers.
CDOC C                          Computed in routine SUMC.
CDOC INTEGER          NOUTF  >  Results file unit identifier.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC DOUBLE_PRECISION RSTAVA >  Array of real state variables other than
CDOC C                          the stress tensor components.
CDOC DOUBLE_PRECISION STRES  >  Array of stress tensor components.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, July 1996: Initial coding
CHST
      SUBROUTINE ORMC
     1(   DGAM       ,NOUTF      ,NTYPE      ,RSTAVA     ,STRES      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(IPHARD=7  ,MSTRE=4)
      DIMENSION  DGAM(2), RSTAVA(MSTRE+1), STRES(*)
      DATA  R2   ,R3    / 2.0D0,3.0D0 /
C***********************************************************************
C OUTPUT RESULTS (INTERNAL AND ALGORITHMIC VARIABLES) FOR THE
C MOHR-COULOMB TYPE ELASTO-PLASTIC MATERIAL WITH ASSOCIATIVE/NON-
C ASSOCIATIVE FLOW RULE AND NON-LINEAR ISOTROPIC HARDENING
C***********************************************************************
 1000 FORMAT(' S-eff = ',G12.4,' Press.= ',G12.4,' Eps.  = ',G12.4,
     1       ' dgama = ',G12.4,' dgamb = ',G12.4)
C
      EPBAR=RSTAVA(MSTRE+1)
      IF(NTYPE.EQ.1)THEN
        P=(STRES(1)+STRES(2))/R3
        EFFST=SQRT(R3/R2*((STRES(1)-P)**2+(STRES(2)-P)**2+
     1                     R2*STRES(3)**2+P**2))
      ELSEIF(NTYPE.EQ.2.OR.NTYPE.EQ.3)THEN
        P=(STRES(1)+STRES(2)+STRES(4))/R3
        EFFST=SQRT(R3/R2*((STRES(1)-P)**2+(STRES(2)-P)**2+
     1                     R2*STRES(3)**2+(STRES(4)-P)**2))
      ENDIF
C Write to output file
      WRITE(NOUTF,1000)EFFST,P,EPBAR,DGAM(1),DGAM(2)
      RETURN
      END
CDOC END_SUBROUTINE ORMC
