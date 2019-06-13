CDOC BEGIN_SUBROUTINE SUTR
CDOC State update procedure for the Tresca elasto-plastic model.
CDOC
CDOC This routine uses the fully implicit elastic predictor/return
CDOC mapping algorithm as the state update procedure for the Tresca
CDOC elasto-plastic material model with piece-wise linear isotropic
CDOC hardening. It contains the plane strain and the axisymmetric
CDOC implementations of the model.
CDOC The essential return mapping is carried out here in the principal
CDOC stress space.
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
CDOC C                          corner and is set to .FALSE.
CDOC C                          otherwise.
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
CHST E.de Souza Neto, May 1996: Initial coding
CHST
CHST E.de Souza Neto, 2005 (??): Associative hardening replaces the
CHST                             originally implemented non-associative
CHST                             hardening.
CHST
      SUBROUTINE SUTR
     1(   DGAM       ,IPROPS     ,LALGVA     ,NTYPE      ,RPROPS     ,
     2    RSTAVA     ,STRAT      ,STRES      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(IPHARD=4  ,MSTRE=4)
C Arguments
      LOGICAL
     1    LALGVA(4)
      DIMENSION
     1    DGAM(2)            ,IPROPS(*)          ,RPROPS(*)          ,
     2    RSTAVA(MSTRE+1)    ,STRAT(MSTRE)       ,STRES(MSTRE)
C Local arrays and variables
      LOGICAL
     1    DUMMY, IFPLAS, RIGHT, SUFAIL, TWOVEC
      DIMENSION
     1    EIGPRJ(MSTRE,2)    ,PSTRS(3)           ,STREST(3)

      DATA
     1    R0   ,R1   ,R2   ,R3   ,R4   ,SMALL ,TOL   / 
     2    0.0D0,1.0D0,2.0D0,3.0D0,4.0D0,1.D-10,1.D-10/
      DATA MXITER / 50 /
C***********************************************************************
C STRESS UPDATE PROCEDURE FOR TRESCA TYPE ELASTO-PLASTIC MATERIAL WITH
C PIECE-WISE LINEAR ISOTROPIC HARDENING:
C IMPLICIT ELASTIC PREDICTOR/RETURN MAPPING ALGORITHM (Boxes 8.1-3).
C PLANE STRAIN AND AXISYMMETRIC IMPLEMENTATIONS.
C
C REFERENCE: Boxes 8.1-3
C            Section 8.1.2
C***********************************************************************
C Stops program if neither plane strain nor axisymmetric state
      IF(NTYPE.NE.2.AND.NTYPE.NE.3)CALL ERRPRT('EI0029')
C Initialize some algorithmic and internal variables
      DGAMA=R0
      DGAMB=R0
      IFPLAS=.FALSE.
      SUFAIL=.FALSE.
      EPBARN=RSTAVA(MSTRE+1)
      EPBAR=EPBARN
C Set some material properties
      YOUNG=RPROPS(2)
      POISS=RPROPS(3)
      NHARD=IPROPS(3)
C Set some constants
      GMODU=YOUNG/(R2*(R1+POISS))
      BULK=YOUNG/(R3*(R1-R2*POISS))
      R2G=R2*GMODU
      R4G=R4*GMODU
      R1D3=R1/R3
C Compute elastic trial state
C ---------------------------
C Volumetric strain and pressure stress
      EEV=STRAT(1)+STRAT(2)+STRAT(4)
      P=BULK*EEV
C Spectral decomposition of the elastic trial deviatoric stress
      EEVD3=EEV*R1D3
      STREST(1)=R2G*(STRAT(1)-EEVD3)
      STREST(2)=R2G*(STRAT(2)-EEVD3)
      STREST(3)=GMODU*STRAT(3)
      CALL SPDEC2(EIGPRJ,PSTRS,DUMMY,STREST)
      PSTRS(3)=R2G*(STRAT(4)-EEVD3)
C Identify maximum (PSTRS1) and minimum (PSTRS3) principal stresses
      II=1
      JJ=1
      PSTRS1=PSTRS(II)
      PSTRS3=PSTRS(JJ)
      DO 10 I=2,3
        IF(PSTRS(I).GE.PSTRS1)THEN
          II=I
          PSTRS1=PSTRS(II)
        ENDIF
        IF(PSTRS(I).LT.PSTRS3)THEN
          JJ=I
          PSTRS3=PSTRS(JJ)
        ENDIF
   10 CONTINUE
      IF(II.NE.1.AND.JJ.NE.1)MM=1
      IF(II.NE.2.AND.JJ.NE.2)MM=2
      IF(II.NE.3.AND.JJ.NE.3)MM=3
      PSTRS2=PSTRS(MM)
C Compute trial yield function and check for plastic consistency
C --------------------------------------------------------------
      SHMAXT=PSTRS1-PSTRS3
      SIGMAY=PLFUN(EPBARN,NHARD,RPROPS(IPHARD))
      PHIA=SHMAXT-SIGMAY
      IF(PHIA/SIGMAY.GT.TOL)THEN
C Plastic step: Apply return mapping
C ==================================
        IFPLAS=.TRUE.
C identify possible two-vector return: right or left of main plane
        SCAPRD=PSTRS1+PSTRS3-PSTRS2*R2
        IF(SCAPRD.GE.R0)THEN
          RIGHT=.TRUE.
        ELSE
          RIGHT=.FALSE.
        ENDIF
C Apply one-vector return mapping first (return to main plane)
C ------------------------------------------------------------
        TWOVEC=.FALSE.
C Start Newton-Raphson iterations
        DO 20 NRITER=1,MXITER
C Compute residual derivative
          DENOM=-R4G-DPLFUN(EPBAR,NHARD,RPROPS(IPHARD))
C Compute Newton-Raphson increment and update variable DGAMA
          DDGAMA=-PHIA/DENOM
          DGAMA=DGAMA+DDGAMA
C Compute new residual
          EPBAR=EPBARN+DGAMA
          SIGMAY=PLFUN(EPBAR,NHARD,RPROPS(IPHARD))
          SHMAX=SHMAXT-R4G*DGAMA
          PHIA=SHMAX-SIGMAY
C Check convergence
          RESNOR=ABS(PHIA/SIGMAY)
          IF(RESNOR.LE.TOL)THEN
C Check validity of one-vector return
            S1=PSTRS1-R2G*DGAMA
            S2=PSTRS2
            S3=PSTRS3+R2G*DGAMA
            DELTA=DMAX1(ABS(S1),ABS(S2),ABS(S3))*SMALL
            IF(S1+DELTA.GE.S2.AND.S2+DELTA.GE.S3)THEN
C converged stress is in the same sextant as trial stress -> 1-vector
C return is valid. Update EPBAR and principal deviatoric stresses
              RSTAVA(MSTRE+1)=EPBAR
              PSTRS1=S1
              PSTRS3=S3
              GOTO 50
            ELSE
C 1-vector return is not valid - go to two-vector procedure
              GOTO 30
            ENDIF
          ENDIF
   20   CONTINUE
C failure of stress update procedure
        SUFAIL=.TRUE.
        CALL ERRPRT('WE0001')
        GOTO 999
   30   CONTINUE
C Apply two-vector return mapping (return to corner - right or left)
C ------------------------------------------------------------------
        TWOVEC=.TRUE.
        DGAMA=R0
        DGABAR=R1
        EPBAR=EPBARN
        SIGMAY=PLFUN(EPBARN,NHARD,RPROPS(IPHARD))
        SHMXTA=PSTRS1-PSTRS3
        IF(RIGHT)THEN
          SHMXTB=PSTRS1-PSTRS2
        ELSE
          SHMXTB=PSTRS2-PSTRS3
        ENDIF
        PHIA=SHMXTA-SIGMAY
        PHIB=SHMXTB-SIGMAY
C Start Newton-Raphson iterations
        DO 40 NRITER=1,MXITER
C Compute residual derivative
          HSLOPE=DPLFUN(EPBAR,NHARD,RPROPS(IPHARD))
          DRVAA=-R4G-HSLOPE
          DRVAB=-R2G-HSLOPE
          DRVBA=-R2G-HSLOPE
          DRVBB=-R4G-HSLOPE
C Compute Newton-Raphson increment and update variables DGAMA and DGAMB
          R1DDET=R1/(DRVAA*DRVBB-DRVAB*DRVBA)
          DDGAMA=(-DRVBB*PHIA+DRVAB*PHIB)*R1DDET
          DDGAMB=(DRVBA*PHIA-DRVAA*PHIB)*R1DDET
          DGAMA=DGAMA+DDGAMA
          DGAMB=DGAMB+DDGAMB
C Compute new residual
          DGABAR=DGAMA+DGAMB
          EPBAR=EPBARN+DGABAR
          SIGMAY=PLFUN(EPBAR,NHARD,RPROPS(IPHARD))
          PHIA=SHMXTA-R2G*(R2*DGAMA+DGAMB)-SIGMAY
          PHIB=SHMXTB-R2G*(DGAMA+R2*DGAMB)-SIGMAY
C Check convergence
          RESNOR=(ABS(PHIA)+ABS(PHIB))/SIGMAY
          IF(RESNOR.LE.TOL)THEN
C Update EPBAR and principal deviatoric stresses
            RSTAVA(MSTRE+1)=EPBAR
            IF(RIGHT)THEN
              PSTRS1=PSTRS1-R2G*(DGAMA+DGAMB)
              PSTRS3=PSTRS3+R2G*DGAMA
              PSTRS2=PSTRS2+R2G*DGAMB
            ELSE
              PSTRS1=PSTRS1-R2G*DGAMA
              PSTRS3=PSTRS3+R2G*(DGAMA+DGAMB)
              PSTRS2=PSTRS2-R2G*DGAMB
            ENDIF
            GOTO 50
          ENDIF
   40   CONTINUE
C failure of stress update procedure
        SUFAIL=.TRUE.
        CALL ERRPRT('WE0001')
        GOTO 999
   50   CONTINUE
C update stress components
C ------------------------
        PSTRS(II)=PSTRS1
        PSTRS(JJ)=PSTRS3
        PSTRS(MM)=PSTRS2
        STRES(1)=PSTRS(1)*EIGPRJ(1,1)+PSTRS(2)*EIGPRJ(1,2)+P
        STRES(2)=PSTRS(1)*EIGPRJ(2,1)+PSTRS(2)*EIGPRJ(2,2)+P
        STRES(3)=PSTRS(1)*EIGPRJ(3,1)+PSTRS(2)*EIGPRJ(3,2)
        STRES(4)=PSTRS(3)+P
C and elastic engineering strain
        RSTAVA(1)=(STRES(1)-P)/R2G+EEVD3
        RSTAVA(2)=(STRES(2)-P)/R2G+EEVD3
        RSTAVA(3)=STRES(3)/GMODU
        RSTAVA(4)=PSTRS(3)/R2G+EEVD3
      ELSE
C Elastic step: update stress using linear elastic law
C ====================================================
        STRES(1)=STREST(1)+P
        STRES(2)=STREST(2)+P
        STRES(3)=STREST(3)
        STRES(4)=PSTRS(3)+P
C elastic engineering strain
        RSTAVA(1)=STRAT(1)
        RSTAVA(2)=STRAT(2)
        RSTAVA(3)=STRAT(3)
        RSTAVA(4)=STRAT(4)
      ENDIF
  999 CONTINUE
C Update algorithmic variables before exit
C ========================================
      DGAM(1)=DGAMA
      DGAM(2)=DGAMB
      LALGVA(1)=IFPLAS
      LALGVA(2)=SUFAIL
      LALGVA(3)=TWOVEC
      LALGVA(4)=RIGHT
      RETURN
      END
CDOC END_SUBROUTINE SUTR
