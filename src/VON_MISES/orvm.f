CDOC BEGIN_SUBROUTINE ORVM
CDOC Output results for the von Mises elasto-plastic material model
CDOC
CDOC This routine writes to the results file the internal and
CDOC algorithmic variables of the von Mises elasto-plastic material
CDOC with non-linear isotropic hardening. The results printed here
CDOC are the von Mises effective stress, the equivalent plastic strain
CDOC and the incremental plastic multiplier obtained by the return
CDOC mapping algorithm of routine SUVM or SUVMPS.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DGAMA  >  Incremental plastic multiplier.
CDOC C                          Computed in routine SUVM or SUVMPS.
CDOC INTEGER          NOUTF  >  Results file unit identifier.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC DOUBLE_PRECISION RSTAVA >  Array of real state variables other than
CDOC C                          the stress tensor components.
CDOC DOUBLE_PRECISION STRES  >  Array of stress tensor components.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, June 1996: Initial coding
CHST
      SUBROUTINE ORVM
     1(   DGAMA      ,NOUTF      ,NTYPE      ,RSTAVA     ,STRES      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(IPHARD=4  ,MSTRE=4)
      DIMENSION  RSTAVA(MSTRE+1), STRES(*)
      DATA   R2   ,R3    / 2.0D0,3.0D0 /
C***********************************************************************
C OUTPUT RESULTS (INTERNAL AND ALGORITHMIC VARIABLES) FOR VON MISES
C TYPE ELASTO-PLASTIC MATERIAL WITH NON-LINEAR ISOTROPIC HARDENING
C***********************************************************************
 1000 FORMAT(' S-eff = ',G12.4,' Eps.  = ',G12.4,' dgama = ',G12.4)
C
      EPBAR=RSTAVA(MSTRE+1)
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
      WRITE(NOUTF,1000)EFFST,EPBAR,DGAMA
      RETURN
      END
CDOC END_SUBROUTINE ORVM
