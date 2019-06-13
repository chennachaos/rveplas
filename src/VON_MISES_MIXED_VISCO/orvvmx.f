CDOC BEGIN_SUBROUTINE ORVVMX
CDOC Output results for a mixed hardening von Mises viscoplastic model
CDOC
CDOC This routine writes to the results file the internal and
CDOC algorithmic variables of the von Mises elasto-viscoplastic material
CDOC with non-linear mixed hardening and Peric's power law for
CDOC viscoplastic flow.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DGAMA  >  Incremental plastic multiplier.
CDOC C                          Computed in routine SUVMMX.
CDOC INTEGER          NOUTF  >  Results file unit identifier.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC DOUBLE_PRECISION RSTAVA >  Array of real state variables other than
CDOC C                          the stress tensor components.
CDOC DOUBLE_PRECISION STRES  >  Array of stress tensor components.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto & F.Adziman, September 2012: Initial coding
CHST
      SUBROUTINE ORVVMX
     1(   DGAMA      ,NOUTF      ,NTYPE      ,RSTAVA     ,STRES      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(IPHARD=6  ,MSTRE=4)
C Arguments
      DIMENSION  RSTAVA(2*MSTRE+1), STRES(*)
C Local arrays and variables
      DIMENSION  BACSTR(MSTRE)
      DATA   R2   ,R3    / 2.0D0,3.0D0 /
C***********************************************************************
C OUTPUT RESULTS (INTERNAL AND ALGORITHMIC VARIABLES) FOR VON MISES
C ELASTO-VISCOPLASTIC MATERIAL WITH NON-LINEAR MIXED HARDENING
C AND PERIC'S POWER LAW FOR VISCOPLASTIC FLOW
C***********************************************************************
 1000 FORMAT(' b-xx  = ',G12.4,' b-yy  = ',G12.4,' b-xy  = ',G12.4,
     1       ' b-zz  = ',G12.4)
 1100 FORMAT(' S-eff = ',G12.4,' Eps.  = ',G12.4,' dgama = ',G12.4)
C
C Plane strain and axisymmetric only
      IF(NTYPE.NE.2.AND.NTYPE.NE.3)CALL ERRPRT('EI0063')
C Print backstress tensor
      BACSTR(1)=RSTAVA(6)
      BACSTR(2)=RSTAVA(7)
      BACSTR(3)=RSTAVA(8)
      BACSTR(4)=RSTAVA(9)
      WRITE(NOUTF,1000)BACSTR(1),BACSTR(2),BACSTR(3),BACSTR(4)
C Compute effective stress
      IF(NTYPE.EQ.2.OR.NTYPE.EQ.3)THEN
        P=(STRES(1)+STRES(2)+STRES(4))/R3
        EFFST=SQRT(R3/R2*((STRES(1)-P)**2+(STRES(2)-P)**2+
     1                     R2*STRES(3)**2+(STRES(4)-P)**2))
      ENDIF
      EPBAR=RSTAVA(MSTRE+1)
      WRITE(NOUTF,1100)EFFST,EPBAR,DGAMA
      RETURN
      END
CDOC END_SUBROUTINE ORVVMX
