CDOC BEGIN_SUBROUTINE OREL
CDOC Output results for the linear elastic material model
CDOC
CDOC This routine writes to the results file some state variables for
CDOC the linear elastic material model.
CDOC
CDOC BEGIN_PARAMETERS
CDOC INTEGER          NOUTF  >  Results file unit identifier.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC DOUBLE_PRECISION STRES  >  Array of stress tensor components.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, September 1996: Initial coding
CHST D. de Bortoli  , March     2015: 3-D case added (NTYPE=4)
CHST
      SUBROUTINE OREL
     1(   NOUTF      ,NTYPE      ,STRES      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION  STRES(*)
      DATA   R2   ,R3    / 2.0D0,3.0D0 /
C***********************************************************************
C OUTPUT RESULTS FOR LINEAR ELASTIC MATERIAL MODEL
C***********************************************************************
 1000 FORMAT(' S-eff = ',G12.4,' Press.= ',G12.4)
C
      IF(NTYPE.EQ.1)THEN
        P=(STRES(1)+STRES(2))/R3
        EFFST=SQRT(R3/R2*((STRES(1)-P)**2+(STRES(2)-P)**2+
     1                     R2*STRES(3)**2+P**2))
      ELSEIF(NTYPE.EQ.2.OR.NTYPE.EQ.3)THEN
        P=(STRES(1)+STRES(2)+STRES(4))/R3
        EFFST=SQRT(R3/R2*((STRES(1)-P)**2+(STRES(2)-P)**2+
     1                     R2*STRES(3)**2+(STRES(4)-P)**2))
      ELSEIF(NTYPE.EQ.4)THEN
        P=(STRES(1)+STRES(2)+STRES(3))/R3
        EFFST=SQRT(R3/R2*((STRES(1)-P)**2+(STRES(2)-P)**2+
     1                    (STRES(3)-P)**2+R2*STRES(4)**2+
     2                     R2*STRES(5)**2+R2*STRES(6)**2))
      ENDIF
C Write to output file
      WRITE(NOUTF,1000)EFFST,P
      RETURN
      END
CDOC END_SUBROUTINE OREL
