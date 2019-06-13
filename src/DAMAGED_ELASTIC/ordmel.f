CDOC BEGIN_SUBROUTINE ORDMEL
CDOC Output results for the damaged elastic/crack closure model.
CDOC
CDOC This routine writes to the results file some state variables for
CDOC the isotropically damaged isotropic elastic model accounting for
CDOC partial microcrack/void closure effects (quasi-unilateral
CDOC conditions) under compressive stresses.
CDOC
CDOC BEGIN_PARAMETERS
CDOC INTEGER          NOUTF  >  Results file unit identifier.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC DOUBLE_PRECISION STRES  >  Array of stress tensor components.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto  July  2001
CHST
      SUBROUTINE ORDMEL
     1(   NOUTF      ,NTYPE      ,STRES      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(MSTRE=4)
      DIMENSION  STRES(*)
      DATA   R2   ,R3    / 2.0D0,3.0D0 /
C***********************************************************************
C OUTPUT RESULTS FOR ISOTROPICALLY DAMAGED ISOTROPIC ELASTIC MODEL
C ACCOUNTING FOR PARTIAL MICROCRACK/VOID CLOSURE EFFECTS
C***********************************************************************
 1000 FORMAT(' S-eff = ',G12.4,' Press.= ',G12.4)
C
      IF(NTYPE.EQ.2.OR.NTYPE.EQ.3)THEN
        P=(STRES(1)+STRES(2)+STRES(4))/R3
        EFFST=SQRT(R3/R2*((STRES(1)-P)**2+(STRES(2)-P)**2+
     1                     R2*STRES(3)**2+(STRES(4)-P)**2))
      ELSE
        CALL ERRPRT('EI0055')
      ENDIF
C Write to output file
      WRITE(NOUTF,1000)EFFST,P
      RETURN
      END
CDOC END_SUBROUTINE ORDMEL
