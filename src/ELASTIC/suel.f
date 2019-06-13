CDOC BEGIN_SUBROUTINE SUEL
CDOC State update procedure for the linear elastic material model.
CDOC
CDOC Given the total strain, this routine computes the corresponding
CDOC stress using the standard generalised Hooke's law for linear
CDOC elastic materials. This routine contains the plane stress, plane
CDOC strain and axisymmetric implementations of the linear elastic
CDOC model.
CDOC
CDOC BEGIN_PARAMETERS
CDOC INTEGER          NTYPE  >  Stress state type flag. The present
CDOC C                          routine is compatible with NTYPE=1
CDOC C                          (plane stress), NTYPE=2 (plane strain)
CDOC C                          and NTYPE=3 (axisymmetric condition).
CDOC DOUBLE_PRECISION RPROPS >  Array of real material properties.
CDOC DOUBLE_PRECISION RSTAVA <  Array of real state variables other
CDOC C                          than the stress tensor components.
CDOC C                          For the linear elastic model, this array
CDOC C                          stores the (engineering) strain
CDOC C                          components.
CDOC DOUBLE_PRECISION STRAN  >  Array of current total (engineering)
CDOC C                          strain components.
CDOC DOUBLE_PRECISION STRES  <  Array of updated stress tensor
CDOC C                          components.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, September 1996: Initial coding
CHST D. de Bortoli  , March     2015: 3-D case added (NTYPE=4)
CHST
      SUBROUTINE SUEL
     1(   NTYPE      ,RPROPS     ,RSTAVA     ,STRAN      ,STRES      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(MSTRE=6)
      DIMENSION
     1    RPROPS(*)          ,RSTAVA(MSTRE)      ,STRAN(*)           ,
     2    STRES(*)    
      DIMENSION
     1    EED(MSTRE)
      DATA
     1    RP5  ,R2   ,R3   ,R4   / 
     2    0.5D0,2.0D0,3.0D0,4.0D0/
C***********************************************************************
C STATE UPDATE PROCEDURE FOR LINEAR ELASTIC MATERIAL MODEL
C
C REFERENCE: Expression (4.43)
C***********************************************************************
C
C Set shear and bulk modulus
C
      GMODU=RPROPS(2)
      BULK=RPROPS(3)
C
C Decompose strain into deviatoric and volumetric components
C ----------------------------------------------------------
C
      R2G=R2*GMODU
      IF(NTYPE.EQ.1)THEN
C for plane stress
        R4G=R4*GMODU
        R4GD3=R4G/R3
        FACTOR=R2G/(BULK+R4GD3)
        EEV=(STRAN(1)+STRAN(2))*FACTOR
        EEVD3=EEV/R3
        EED(1)=STRAN(1)-EEVD3
        EED(2)=STRAN(2)-EEVD3
      ELSEIF(NTYPE.EQ.2.OR.NTYPE.EQ.3)THEN
C for plane strain and axisymmetric cases
        EEV=STRAN(1)+STRAN(2)+STRAN(4)
        EEVD3=EEV/R3
        EED(1)=STRAN(1)-EEVD3
        EED(2)=STRAN(2)-EEVD3
        EED(4)=STRAN(4)-EEVD3
      ELSEIF(NTYPE.EQ.4)THEN
C three-dimensional case
        EEV=STRAN(1)+STRAN(2)+STRAN(3)
        EEVD3=EEV/R3
        EED(1)=STRAN(1)-EEVD3
        EED(2)=STRAN(2)-EEVD3
        EED(3)=STRAN(3)-EEVD3
      ELSE
        CALL ERRPRT('EI0018')
      ENDIF
C Convert engineering shear component into physical component
      IF((NTYPE.EQ.1).OR.(NTYPE.EQ.2).OR.(NTYPE.EQ.3))THEN
        EED(3)=STRAN(3)*RP5
      ELSEIF(NTYPE.EQ.4)THEN
        EED(4)=STRAN(4)*RP5
        EED(5)=STRAN(5)*RP5
        EED(6)=STRAN(6)*RP5
      ENDIF
C
C Update stress using linear elastic law
C ---------------------------------------
C
C hydrostatic stress
      P=BULK*EEV
C stress tensor components
      IF((NTYPE.EQ.1).OR.(NTYPE.EQ.2).OR.(NTYPE.EQ.3))THEN
        STRES(1)=R2G*EED(1)+P
        STRES(2)=R2G*EED(2)+P
        STRES(3)=R2G*EED(3)
        IF(NTYPE.EQ.2.OR.NTYPE.EQ.3)STRES(4)=R2G*EED(4)+P
      ELSEIF(NTYPE.EQ.4)THEN
        STRES(1)=R2G*EED(1)+P
        STRES(2)=R2G*EED(2)+P
        STRES(3)=R2G*EED(3)+P
        STRES(4)=R2G*EED(4)
        STRES(5)=R2G*EED(5)
        STRES(6)=R2G*EED(6)
      ENDIF
C
C Store elastic engineering strain in RSTAVA
C ------------------------------------------
C
      RSTAVA(1)=STRAN(1)
      RSTAVA(2)=STRAN(2)
      RSTAVA(3)=STRAN(3)
      IF(NTYPE.EQ.1)THEN
        R3BULK=R3*BULK
        RSTAVA(4)=-(STRAN(1)+STRAN(2))*(R3BULK-R2G)/(R3BULK+R4G)
      ELSEIF(NTYPE.EQ.2.OR.NTYPE.EQ.3)THEN
        RSTAVA(4)=STRAN(4)
      ELSEIF(NTYPE.EQ.4)THEN
        RSTAVA(4)=STRAN(4)
        RSTAVA(5)=STRAN(5)
        RSTAVA(6)=STRAN(6)
      ENDIF
C
      RETURN
      END
CDOC END_SUBROUTINE SUEL
