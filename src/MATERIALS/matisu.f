CDOC BEGIN_SUBROUTINE MATISU
CDOC Material interface for state update routine calls
CDOC
CDOC This routine calls the state update routine according to the
CDOC material type. Given the material-independent kinematic quantities
CDOC computed at the element level (and passed into the present routine
CDOC through its list of arguments) this routine identifies the material
CDOC type in question and calls the corresponding material-specific
CDOC state update routine which updates the stress and other state
CDOC variables at one Gauss point.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DETF   >  Determinant of the total deformation
CDOC C                          gradient at the current Gauss point.
CDOC C                          Used for large strain analysis only.
CDOC INTEGER          NLARGE >  Large strain flag. Large strain
CDOC C                          analysis if NLARGE=1 and
CDOC C                          infinitesimal strain analysis otherwise.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC LOGICAL          SUFAIL <  State update failure flag. Return value
CDOC C                          set to .FALSE. if the state update
CDOC C                          procedure was performed successfully
CDOC C                          and set to .TRUE. otherwise.
CDOC DOUBLE_PRECISION THKGP  <> Thickness of the current Gauss point.
CDOC C                          Used only in plane stress analysis.
CDOC C                          Updated here only in large strain
CDOC C                          analysis. In this case, it is the
CDOC C                          initial (reference config.) thickness on
CDOC C                          entry and returns as the current
CDOC C                          thickness.
CDOC DOUBLE_PRECISION DTIME  >  Time increment.
CDOC DOUBLE_PRECISION EINCR  >  Array of incremental engineering strain
CDOC C                          components. Used in infinitesimal strain
CDOC C                          analysis only.
CDOC DOUBLE_PRECISION FINCR  >  Incremental deformation gradient. Used
CDOC C                          in large strain analysis only.
CDOC INTEGER          IPROPS >  Array of integer material properties.
CDOC LOGICAL          LALGVA <> Array of current logical algorithmic
CDOC C                          variables at the current Gauss point.
CDOC DOUBLE_PRECISION RALGVA <> Array of current real algorithmic
CDOC C                          variables at the current Gauss point.
CDOC C                          Previous converged (equilibrium) value
CDOC C                          on entry. Returns as updated value.
CDOC DOUBLE_PRECISION RPROPS >  Array of real material properties.
CDOC DOUBLE_PRECISION RSTAVA <> Array of current real state variables
CDOC C                          at the current Gauss point. Previous
CDOC C                          converged (equilibrium) value on entry.
CDOC C                          Returns as updated value.
CDOC DOUBLE_PRECISION STRES  <> Array of current (Cauchy) stress
CDOC C                          components at the current Gauss point.
CDOC C                          Previous converged (equilibrium) value 
CDOC C                          on entry. Returns as updated value.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, July      1999: Initial coding
CHST
CHST E.de Souza Neto & F.Adziman,
CHST                  September 2012: Visco-plastic material added.
CHST
CHST F. Adziman,      September 2013: Passing DVOLU into MATISU
CHST                                  for calc. of volume fraction
CHST                                  in SUMEPC 
CHST
CHST F. Adziman,      October 2013  : Passing IELEM and IGAUSP into
CHST                                  MATISU for multivariant calc.
CHST                                  in SUMEPC 
CHST
      SUBROUTINE MATISU
     1(   DETF       ,NLARGE     ,NTYPE      ,SUFAIL     ,THKGP      ,
     3    DTIME      ,EINCR      ,FINCR      ,IPROPS     ,LALGVA     ,
     4    RALGVA     ,RPROPS     ,RSTAVA     ,STRES      ,TTIME      ,
     5    FTOTA      ,DVOLU      ,IELEM      ,IGAUSP     ,STREP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE '../MATERIAL.INC'
C
      PARAMETER( MSTRA=6 )
C Arguments
      LOGICAL
     1    SUFAIL             ,LALGVA
      DIMENSION
     1    EINCR(*)           ,FINCR(3,3)         ,IPROPS(*)          ,
     2    LALGVA(*)          ,RALGVA(*)          ,RPROPS(*)          ,
     3    RSTAVA(*)          ,STRES(*)           ,FTOTA(3,3)
C Local arrays
      DIMENSION
     1    B(MSTRA)           ,BETRL(MSTRA)       ,STRAT(MSTRA)
C Local numerical constants
      DATA
     1    R1   /
     2    1.0D0/
C***********************************************************************
C MATERIAL INTERFACE FOR STATE UPDATE ROUTINE CALLS:
C ACCORDING TO THE MATERIAL TYPE, CALLS MATERIAL-SPECIFIC STATE UPDATE
C ROUTINE TO UPDATE STRESS AND OTHER STATE VARIABLES
C
C REFERENCE: Figure 5.4
C            Sections 5.7.2, 5.7.6
C***********************************************************************
C Set up number of stress components
      IF(NTYPE.EQ.1)THEN
        NSTRE=3
      ELSEIF(NTYPE.EQ.2)THEN
        NSTRE=4
      ELSEIF(NTYPE.EQ.3)THEN
        NSTRE=4
      ELSEIF(NTYPE.EQ.4)THEN
        NSTRE=6
      ELSE
        CALL ERRPRT('EI0040')
      ENDIF
C Identify material type and class
      MATTYP=IPROPS(1)
      MATCLS=IPROPS(2)
C
C Then call material class/type-specific routines
C
      IF(MATCLS.EQ.HYPEPL)THEN
C
C Isotropic elastic/elasto-plastic materials with logarithmic finite
C strain extension
C ==================================================================
C
C Compute elastic trial strains. Note that for the purely elastic models
C the elastic trial strain equals the total strain
C ----------------------------------------------------------------------
        IF(NLARGE.EQ.0)THEN
C Small strains: compute elastic trial INFINITESIMAL strain
          DO 10 ISTRE=1,NSTRE
            STRAT(ISTRE)=RSTAVA(ISTRE)+EINCR(ISTRE)
   10     CONTINUE
        ELSEIF(NLARGE.EQ.1)THEN
C Large strains: compute elastic trial LOGARITHMIC strain
C... elastic trial left Cauchy-Green tensor
          CALL BETRIA
     1(   BETRL      ,RSTAVA     ,FINCR       ,NTYPE      )
C... elastic trial eulerian logarithmic strain
          CALL LOGSTR
     1(   BETRL      ,STRAT      ,NTYPE       )
        ENDIF
C Apply small strain material type-specific state updating procedure
C ------------------------------------------------------------------
        IF(MATTYP.EQ.ELASTC)THEN
C Linear elastic (Hencky material in large strains)
          CALL SUEL
     1(   NTYPE      ,RPROPS     ,RSTAVA     ,STRAT      ,STRES      )
C...set elasto-plastic flag and state update failure flag
          LALGVA(1)=.FALSE.
          LALGVA(2)=.FALSE.
        ELSEIF(MATTYP.EQ.TRESCA)THEN
C Tresca elasto-plastic
          IF(NTYPE.EQ.1)THEN
            CALL SUTRPN
     1(   RALGVA     ,IPROPS     ,LALGVA     ,NTYPE      ,RPROPS     ,
     2    RSTAVA     ,STRAT      ,STRES      )
          ELSEIF(NTYPE.EQ.2.OR.NTYPE.EQ.3)THEN
            CALL SUTR
     1(   RALGVA     ,IPROPS     ,LALGVA     ,NTYPE      ,RPROPS     ,
     2    RSTAVA     ,STRAT      ,STRES      )
          ENDIF
        ELSEIF(MATTYP.EQ.VMISES)THEN
C von Mises elasto-plastic
          IF(NTYPE.EQ.1)THEN
            CALL SUVMPS
     1(   RALGVA     ,IPROPS     ,LALGVA     ,NTYPE      ,RPROPS     ,
     2    RSTAVA     ,STRAT      ,STRES      )
          ELSEIF(NTYPE.EQ.2.OR.NTYPE.EQ.3)THEN
            CALL SUVM
     1(   RALGVA     ,IPROPS     ,LALGVA     ,NTYPE      ,RPROPS     ,
     2    RSTAVA     ,STRAT      ,STRES      )
          ENDIF
        ELSEIF(MATTYP.EQ.MOHCOU)THEN
C Mohr-Coulomb elasto-plastic
          CALL SUMC
     1(   RALGVA     ,IPROPS     ,LALGVA     ,NTYPE      ,RPROPS     ,
     2    RSTAVA     ,STRAT      ,STRES      )
        ELSEIF(MATTYP.EQ.DRUPRA)THEN
C Drucker-Prager elasto-plastic
          IF(NTYPE.EQ.4)THEN
            CALL SUDP3D
     1(   RALGVA     ,IPROPS     ,LALGVA     ,NTYPE      ,RPROPS     ,
     2    RSTAVA     ,STRAT       ,STRES     ,STREP  )
          ELSEIF(NTYPE.EQ.2.OR.NTYPE.EQ.3)THEN
            CALL SUDP
     1(   RALGVA     ,IPROPS     ,LALGVA     ,NTYPE      ,RPROPS     ,
     2    RSTAVA     ,STRAT      ,STRES      )
          ENDIF
        ELSEIF(MATTYP.EQ.LEMDAM)THEN
C Lemaitre's ductile damage elasto-plastic model
          CALL SUDAMA
     1(   RALGVA     ,IPROPS     ,LALGVA     ,NTYPE      ,RPROPS     ,
     2    RSTAVA     ,STRAT      ,STRES      )
        ELSEIF(MATTYP.EQ.DAMELA)THEN
C Isotropically damaged isotropic elastic material with crack closure
C effects
          CALL SUDMEL
     1(   NTYPE      ,RPROPS     ,RSTAVA     ,STRAT      ,STRES      ,
     2    LALGVA(2)  )
C...set elasto-plastic flag to false
          LALGVA(1)=.FALSE.
        ELSE
C Error: Material type not recognised
          CALL ERRPRT('EI0042')
        ENDIF
C Exit routine in case of failure of the state update procedure
        SUFAIL=LALGVA(2)
        IF(SUFAIL)GOTO 999
        IF(NLARGE.EQ.1)THEN
C Perform extra updating operations required by this class of
C elasto-pastic materials at large strains only
C -----------------------------------------------------------
          DETFT=DETF
          IF(NTYPE.EQ.1)THEN
C Plane stress: update Gauss point thickness according to material
C model. Also update the total deformation gradient (taking the
C thickness strain into account)
            IF(MATTYP.EQ.ELASTC)THEN
C... Elastic (Hencky material in large strains)
              CALL TUEL
     1(   DETFT      ,RSTAVA     ,THKGP      ,1          )
            ELSEIF(MATTYP.EQ.VMISES)THEN
C... von Mises elasto-plastic
              CALL TUVM
     1(   DETFT      ,RSTAVA     ,THKGP      ,1          )
            ELSE
C... Error: Material type not recognised or not implemented for finite
C    strains under plane stress
              CALL ERRPRT('EI0058')
            ENDIF
          ENDIF
C Transform Kirchhoff into Cauchy stress
          DETFIN=R1/DETFT
          CALL RVSCAL(STRES,NSTRE,DETFIN)
        ENDIF
C
      ELSEIF(MATCLS.EQ.SINCRY)THEN
C
C Single crystal anisotropic finite elasto-plastic models
C =======================================================
C
        IF(MATTYP.EQ.PDSCRY)THEN
C Planar double slip single crystal
          CALL SUPDSC
     1(   RALGVA     ,FINCR      ,IPROPS     ,LALGVA     ,NTYPE      ,
     2    RPROPS     ,RSTAVA     ,STRES      )	 
C Single crystal anisotropic finite viscoplastic models
C =======================================================
C
        ELSEIF(MATTYP.EQ.VSCTWD)THEN
C Planar double slip single crystal
          CALL SUVSC2
     1(   DTIME      ,FINCR      ,IPROPS     ,LALGVA     ,NTYPE      ,
     2    RPROPS     ,RSTAVA     ,STRES      )
C Single crystal martensitic transformation visco-plastic models
C ===================================================================
C
        ELSEIF(MATTYP.EQ.MTCTWD)THEN
C Single crystal martensitic transformation
          CALL SUMTSC
     1(   RALGVA     ,FINCR      ,IPROPS     ,LALGVA     ,NTYPE      ,
     2    RPROPS     ,RSTAVA     ,STRES      ,DVOLU      ,IELEM      ,
     3    IGAUSP     ,DTIME      )
C Single crystal martensitic transformation elasto-plastic models
C ===================================================================
C	 
        ELSEIF(MATTYP.EQ.MTEPTD)THEN         
C Single crystal martensitic transformation
          CALL SUMEPC
     1(   RALGVA     ,FINCR      ,IPROPS     ,LALGVA     ,NTYPE      ,
     2    RPROPS     ,RSTAVA     ,STRES      ,DVOLU      ,IELEM      ,
     3    IGAUSP     ,DTIME      )
        ELSE
C Error: Material type not recognised
          CALL ERRPRT('EI0042')
        ENDIF
        SUFAIL=LALGVA(2)
        IF(SUFAIL)GOTO 999
      ELSEIF(MATCLS.EQ.HYPER)THEN
C
C Generic isotropic finite hyperelasticity models
C ===============================================
C
C First compute current Left Cauchy-Green strain tensor, B
        CALL LEFTCG
     1(   RSTAVA     ,B          ,FINCR      ,NTYPE      )
C Then call the material type-specific state update procedure
        IF(MATTYP.EQ.OGDEN)THEN
C Ogden model
          CALL SUOGD
     1(   B          ,IPROPS     ,NTYPE      ,RPROPS     ,RSTAVA     ,
     2    STRES      ,THKGP      )
          SUFAIL=.FALSE.
        ELSE
C Error: Material type not recognised
          CALL ERRPRT('EI0042')
        ENDIF
      ELSEIF(MATCLS.EQ.PLASTC)THEN
C
C Elasto-plastic/visoplastic materials with small strain implementation
C only
C =====================================================================
C
C compute elastic trial INFINITESIMAL strain
        DO 20 ISTRE=1,NSTRE
          STRAT(ISTRE)=RSTAVA(ISTRE)+EINCR(ISTRE)
   20   CONTINUE
        IF(MATTYP.EQ.VMMIXD)THEN
C von Mises with mixed isotropic/kinematic hardening
          CALL SUVMMX
     1(   RALGVA     ,IPROPS     ,LALGVA     ,NTYPE      ,RPROPS     ,
     2    RSTAVA     ,STRAT      ,STRES      )
        ELSEIF(MATTYP.EQ.VVMIXD)THEN
C von Mises viscoplastic with mixed isotropic/kinematic hardening and
C Peric's power law for visoplastic flow
          CALL SUVVMX
     1(   RALGVA     ,DTIME      ,IPROPS     ,LALGVA     ,NTYPE      ,
     2    RPROPS     ,RSTAVA     ,STRAT      ,STRES      )
        ELSE
C Error: Material type not recognised
          CALL ERRPRT('EI0042')
        ENDIF
        SUFAIL=LALGVA(2)
        IF(SUFAIL)GOTO 999
      ELSE
C Error: Material class not recognised
        CALL ERRPRT('EI0041')
      ENDIF
  999 CONTINUE
      RETURN
      END
CDOC END_SUBROUTINE MATISU
