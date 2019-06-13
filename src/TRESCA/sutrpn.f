CDOC BEGIN_SUBROUTINE SUTRPN
CDOC State update procedure for the Tresca model in plane stress.
CDOC
CDOC This routine uses the fully implicit elastic predictor/return
CDOC mapping algorithm as the state update procedure for the
CDOC Tresca elasto-plastic material model with piece-wise linear
CDOC isotropic hardening under plane stress condition.
CDOC The algorithm used here is based on the nested iteration approach
CDOC for enforcement of the plane stress constraint at the Gauss point
CDOC level.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DGAM   <  Array of incremental plastic
CDOC C                          multipliers.
CDOC INTEGER          IPROPS >  Array of integer material properties.
CDOC C                          The number of points on the piece-wise
CDOC C                          linear hardening curve is the only
CDOC C                          element of this array used here.
CDOC C                          This array is set in routines
CDOC C                          INDATA and RDTR.
CDOC LOGICAL          LALGVA <  Array of logical algorithmic flags.
CDOC C                          For the present material model, this
CDOC C                          array contains the plastic yielding
CDOC C                          flag, IFPLAS; the return
CDOC C                          algorithm failure flag, SUFAIL;
CDOC C                          the two-vector return flag,
CDOC C                          TWOVEC and the right corner return flag,
CDOC C                          RIGHT.
CDOC C                          The plastic yielding flag is set
CDOC C                          to .TRUE. if plastic yielding
CDOC C                          has occurred and to .FALSE. if
CDOC C                          the step is elastic. The algorithm
CDOC C                          failure flag is set to .FALSE. if
CDOC C                          the state update algorithm has been
CDOC C                          successful and to .TRUE. if the
CDOC C                          return mapping algorithm has failed
CDOC C                          to converge.
CDOC C                          TWOVEC is set to .TRUE.
CDOC C                          if the selected return mapping is to
CDOC C                          a corner (right or left) and is set to
CDOC C                          .FALSE. otherwise.
CDOC C                          RIGHT is set to .TRUE.
CDOC C                          if the selected return is to the right
CDOC C                          corner and is set to .FALSE. otherwise.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC DOUBLE_PRECISION RPROPS >  Array of real material properties.
CDOC C                          This array contains the density
CDOC C                          (not used in this routine), the elastic
CDOC C                          properties: Young's modulus and
CDOC C                          Poisson's ratio, and the plastic
CDOC C                          properties: the pairs
CDOC C                          ``accumulated plastic strain-cohesion''
CDOC C                          defining the (user
CDOC C                          supplied) piece-wise linear hardening
CDOC C                          curve. This array is set in routine
CDOC C                          RDTR.
CDOC DOUBLE_PRECISION RSTAVA <> Array of real state variables other
CDOC C                          than the stress tensor components.
CDOC C                          Previous converged values on entry,
CDOC C                          updated values on exit.
CDOC C                          The state variables stored in
CDOC C                          this array are the (engineering)
CDOC C                          elastic strain components and the
CDOC C                          accumulated plastic strain.
CDOC DOUBLE_PRECISION STRAT  >  Array of elastic trial (engineering)
CDOC C                          strain components.
CDOC DOUBLE_PRECISION STRES  <  Array of updated stress tensor
CDOC C                          components.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, Sept 1998: Initial coding
CHST
      SUBROUTINE SUTRPN
     1(   DGAM       ,IPROPS     ,LALGVA     ,NTYPE      ,RPROPS     ,
     2    RSTAVA     ,STRAT      ,STRES      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(IPHARD=4  ,MSTRE=4)
C Arguments
      LOGICAL
     1    LALGVA(4)
      DIMENSION
     1    DGAM(2)            ,IPROPS(*)          ,RPROPS(*)          ,
     2    RSTAVA(MSTRE+1)    ,STRAT(*)           ,STRES(*)
C Local arrays and variables
      LOGICAL EPFLAG, IFPLAS, SUFAIL
      DIMENSION
     1    DMATX(MSTRE,MSTRE) ,RSTAUX(MSTRE+1)
      DATA
     1    R0    ,TOL   / 
     2    0.D0  ,1.D-08/
      DATA MXITER / 20 /
C***********************************************************************
C STRESS UPDATE PROCEDURE FOR TRESCA TYPE ELASTO-PLASTIC MATERIAL WITH
C PIECE-WISE LINEAR ISOTROPIC HARDENING IN PLANE STRESS. NESTED
C ITERATION APPROACH.
C
C REFERENCE: Section 9.2.2
C***********************************************************************
C Stop program if not plane stress
      IF(NTYPE.NE.1)CALL ERRPRT('EI0060')
C Initialise state update failure flag
      SUFAIL=.FALSE.
C Set some material properties
      NHARD=IPROPS(3)
C
C Newton-Raphson iteration loop for plane stress enforcement
C
C Set initial guess for elastic trial thickness strain. Use previously
C converged elastic trial thickness strain.
      E33TRL=RSTAVA(4)
C Set axisymmetric state flag
      NTYPAX=3
C Start N-R loop
      DO 20 ITER=1,MXITER
C Set state variables to values at beginning of increment
        DO 10 I=1,MSTRE+1
          RSTAUX(I)=RSTAVA(I)
   10   CONTINUE
C Use axisymmetric integration algorithm to compute stresses
        STRAT(4)=E33TRL
        CALL SUTR
     1(   DGAM       ,IPROPS     ,LALGVA     ,NTYPAX     ,RPROPS     ,
     2    RSTAUX     ,STRAT      ,STRES      )
        SUFAIL=LALGVA(2)
C ...emergency exit in case of failure of state update procedure
	IF(SUFAIL)THEN
	  GOTO 999
	ENDIF
        IFPLAS=LALGVA(1)
C Check plane stress convergence
        EPBAR=RSTAUX(MSTRE+1)
	SIGMAY=PLFUN(EPBAR,NHARD,RPROPS(IPHARD))
	RES=ABS(STRES(4))
C ...use normalised out-of-plane stress
        IF(SIGMAY.NE.R0)RES=RES/SIGMAY
        IF(RES.LE.TOL)THEN
C ...and break N-R loop in case of convergence
          GOTO 30
        ENDIF
C Compute axisymmetric consistent tangent components
        EPFLAG=IFPLAS
        CALL CTTR
     1(   DMATX      ,EPFLAG     ,IPROPS     ,LALGVA     ,NTYPAX     ,
     2    RPROPS     ,RSTAUX     ,STRAT      ,STRES      )
C Apply Newton-Raphson correction to normal elastic trial strain
        D22=DMATX(4,4)
        E33TRL=E33TRL-STRES(4)/D22
   20 CONTINUE
C Emergency exit in case of failure of plane stress enforcement loop
      SUFAIL=.TRUE.
      CALL ERRPRT('WE0026')
      GOTO 999
   30 CONTINUE
C Set state variables to current updated values
      DO 40 I=1,MSTRE+1
        RSTAVA(I)=RSTAUX(I)
   40 CONTINUE
      RSTAVA(4)=E33TRL
C
  999 CONTINUE
      RETURN
      END
CDOC END_SUBROUTINE SUTRPN
