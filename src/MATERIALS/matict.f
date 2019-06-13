CDOC BEGIN_SUBROUTINE MATICT
CDOC Material interface for consistent tangent routine calls
CDOC
CDOC This routine calls the consistent tangent computation routine
CDOC according to the material type. Given the necessary
CDOC material-independent kinematic quantities computed at the element
CDOC level (and passed into the present routine through its list of
CDOC arguments) this routine identifies the material type in question
CDOC and calls the corresponding material-specific tangent computation
CDOC routine.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DETF   >  Current total deformation gradient at
CDOC C                          the Gauss point of interest.
CDOC INTEGER          IITER  >  Number of the current equilibrium
CDOC C                          iteration.
CDOC INTEGER          KUNLD  >  Unloading flag. KUNLD is set to 1 if the
CDOC C                          loading programme is currently
CDOC C                          unloading.
CDOC INTEGER          MBDIM  >  Dimension of DMATX.
CDOC INTEGER          MGDIM  >  Dimension of AMATX.
CDOC INTEGER          MSTRE  >  Maximum number of stress components.
CDOC INTEGER          NLARGE >  Large deformation flag. Large
CDOC C                          deformation analysis if
CDOC C                          NLARGE=1 and infinitesimal deformation
CDOC C                          analysis otherwise.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC DOUBLE_PRECISION AMATX  <  Consistent spatial tangent modulus.
CDOC DOUBLE_PRECISION DMATX  <  Consistent infinitesimal tangent
CDOC C                          modulus.
CDOC DOUBLE_PRECISION DTIME  >  Time increment.
CDOC DOUBLE_PRECISION EINCR  >  Current incremental engineering strain
CDOC C                          components. Used only in small strain
CDOC C                          analysis.
CDOC DOUBLE_PRECISION FINCR  >  Current incremental deformation
CDOC C                          gradient. Used only in large strain
CDOC C                          analysis.
CDOC INTEGER          IPROPS >  Array of integer material properties.
CDOC LOGICAL          LALGVA >  Array of current logical algorithmic
CDOC C                          variables at the Gauss point of
CDOC C                          interest.
CDOC DOUBLE_PRECISION RALGVA >  Array of current real algorithmic
CDOC C                          variables at the Gauss points of
CDOC C                          interest.
CDOC DOUBLE_PRECISION RPROPS >  Array of real material properties.
CDOC DOUBLE_PRECISION RSTAVA >  Array of current real state variables
CDOC C                          at the Gauss points of interest.
CDOC DOUBLE_PRECISION RSTAV2 >  Array of real state variables at the
CDOC C                          Gauss point of interest at the previous
CDOC C                          equilibrium (converged) configuration.
CDOC DOUBLE_PRECISION STRES  >  Array of current (Cauchy) stress
CDOC C                          components at the current Gauss point.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, July      1999: Initial coding (taken from STSTD2)
CHST
CHST E.de Souza Neto, October   2008: Unused dummy argument IITER
CHST                                  removed
CHST
CHST E.de Souza Neto & F.Adziman,
CHST                  September 2012: Time increment added
CHST
CHST D. de Bortoli,        July 2014: Added MSTRE to input parameters
CHST
      SUBROUTINE MATICT
     1(   DETF       ,KUNLD      ,MBDIM      ,MGDIM      ,MSTRE      ,
     2    NLARGE     ,NTYPE      ,
     3    AMATX      ,DMATX      ,DTIME      ,EINCR      ,FINCR      ,
     4    IPROPS     ,LALGVA     ,RALGVA     ,RPROPS     ,RSTAVA     ,
     5    RSTAV2     ,STRES      ,FINV       )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE '../MATERIAL.INC'
C
      PARAMETER( MSTRA=6 )
C Arguments
      LOGICAL  LALGVA
      DIMENSION
     1    AMATX(MGDIM,MGDIM) ,DMATX(MBDIM,MBDIM) ,EINCR(MBDIM)       ,
     2    FINCR(3,3)         ,IPROPS(*)          ,LALGVA(*)          ,
     3    RALGVA(*)          ,RPROPS(*)          ,RSTAVA(*)          ,
     4    RSTAV2(*)          ,STRES(*)           ,FINV(3,3)         
C Local arrays and variables
      LOGICAL EPFLAG ,IFPLAS, CTFAIL
      DIMENSION
     1    BETRL(MSTRA)       ,STRAT(MSTRA)       ,STRESK(MSTRE)
C***********************************************************************
C MATERIAL INTERFACE FOR CONSISTENT TANGENT COMPUTATION ROUTINE CALLS:
C ACCORDING TO THE MATERIAL TYPE, CALLS MATERIAL-SPECIFIC TANGENT
C COMPUTATION ROUTINE
C
C REFERENCE: Figure 5.5
C            Sections 5.7.4, 5.7.6
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
        CALL ERRPRT('EI0076')
      ENDIF
C Start by identifying the material type and class
      MATTYP=IPROPS(1)
      MATCLS=IPROPS(2)
C
C Then call material class/type-specific routines
C
      IF(MATCLS.EQ.HYPEPL)THEN
C
C Elastic/elasto-plastic materials with logarithmic finite strain
C extension
C ===============================================================
C
C Retrieve current stress
        DO 50 ISTRE=1,NSTRE
          STRESK(ISTRE)=STRES(ISTRE)
   50   CONTINUE
        IF(NLARGE.EQ.1)THEN
C Large strains: compute last elastic trial LOGARITHMIC strain
C ------------------------------------------------------------
C elastic trial left Cauchy-Green tensor
          CALL BETRIA
     1(   BETRL      ,RSTAV2     ,FINCR      ,NTYPE      )
C elastic trial eulerian logarithmic strain
          CALL LOGSTR
     1(   BETRL      ,STRAT      ,NTYPE      )
          DETFT=DETF
          IF(NTYPE.EQ.1)THEN
C compute total deformation gradient (including thickness strain
C contribution) for plane stress, according to material model
            IF(MATTYP.EQ.ELASTC)THEN
C... Elastic (Hencky material in large strain)
              CALL TUEL
     1(   DETFT      ,RSTAVA     ,DUMMY      ,0          )
            ELSEIF(MATTYP.EQ.VMISES)THEN
C... von Mises elasto-plastic
              CALL TUVM
     1(   DETFT      ,RSTAVA     ,DUMMY      ,0          )
            ELSE
C... Error: Material type not recognised or not implemented for finite
C    strains under plane stress
              CALL ERRPRT('EI0059')
            ENDIF
          ENDIF
C retrieve current KIRCHHOFF stress in large strains
          CALL RVSCAL(STRESK,NSTRE,DETFT)
        ELSE
C Small strains: compute last elastic trial INFINITESIMAL strain
C --------------------------------------------------------------
          DO 60 ISTRE=1,NSTRE
            STRAT(ISTRE)=RSTAV2(ISTRE)+EINCR(ISTRE)
   60     CONTINUE
        ENDIF
C Set plastic elasto-plastic tangent flag
C ---------------------------------------
        IF(MATTYP.NE.ELASTC)THEN
          IFPLAS=LALGVA(1)
          IF((.NOT.IFPLAS).OR.KUNLD.EQ.1)THEN
            EPFLAG=.FALSE.
          ELSE
            EPFLAG=.TRUE.
          ENDIF
        ENDIF
C Call material type-specific routines
C
C WARNING:
C Materials routines must be called using the appropriate array section
C of matrices DMATX and AMATX, so that the dummy variable dimensions
C (hard-coded inside the respective subroutines) match the actual
C variable dimensions. 
C For instance, CTTRPN requires DMATX to have dimensions (4,4), so it is
C called using the array section DMATX(1:4,1:4)
C ----------------------------------------------------------------------
        IF(MATTYP.EQ.ELASTC)THEN
C Elastic
          CALL CTEL
     1(   DMATX(1:6,1:6),NTYPE      ,RPROPS     )
        ELSEIF(MATTYP.EQ.TRESCA)THEN
C Tresca
          IF(NTYPE.EQ.1)THEN
            CALL CTTRPN
     1(   DMATX(1:4,1:4),EPFLAG     ,IPROPS     ,LALGVA     ,RPROPS    ,
     2    RSTAVA        ,STRAT      ,STRESK     )
          ELSEIF(NTYPE.EQ.2.OR.NTYPE.EQ.3)THEN
            CALL CTTR
     1(   DMATX(1:4,1:4),EPFLAG     ,IPROPS     ,LALGVA     ,NTYPE     ,
     2    RPROPS        ,RSTAVA     ,STRAT      ,STRESK     )
          ENDIF
        ELSEIF(MATTYP.EQ.VMISES)THEN
C von Mises
          IF(NTYPE.EQ.1)THEN
            CALL CTVMPS
     1(   RALGVA     ,DMATX(1:4,1:4),EPFLAG     ,IPROPS     ,NTYPE     ,
     2    RPROPS     ,RSTAVA        ,STRESK     )
          ELSEIF(NTYPE.EQ.2.OR.NTYPE.EQ.3)THEN
            CALL CTVM
     1(   RALGVA     ,DMATX(1:4,1:4),EPFLAG     ,IPROPS     ,NTYPE     ,
     2    RPROPS     ,RSTAVA     ,STRESK     )
          ENDIF
        ELSEIF(MATTYP.EQ.MOHCOU)THEN
C Mohr-Coulomb
          CALL  CTMC
     1(   DMATX(1:4,1:4),EPFLAG     ,IPROPS     ,LALGVA     ,NTYPE     ,
     2    RPROPS        ,RSTAVA     ,STRAT      ,STRESK     )
        ELSEIF(MATTYP.EQ.DRUPRA)THEN
C Drucker-Prager
          IF(NTYPE.EQ.4)THEN
            CALL  CTDP3D
     1(   RALGVA     ,DMATX         ,EPFLAG     ,IPROPS     ,LALGVA    ,
     2    NTYPE      ,RPROPS        ,RSTAVA     ,STRAT      )
          ELSEIF(NTYPE.EQ.2.OR.NTYPE.EQ.3)THEN
            CALL  CTDP
     1(   RALGVA     ,DMATX(1:4,1:4),EPFLAG     ,IPROPS     ,LALGVA    ,
     2    NTYPE      ,RPROPS        ,RSTAVA     ,STRAT      )
          ENDIF
        ELSEIF(MATTYP.EQ.LEMDAM)THEN
C Lemaitre's ductile damage model
          CALL CTDAMA
     1(   RALGVA     ,DMATX(1:4,1:4),EPFLAG     ,IPROPS     ,NTYPE     ,
     2    RPROPS     ,RSTAVA        ,STRESK     )
        ELSEIF(MATTYP.EQ.DAMELA)THEN
C Isotropically damaged isotropic elastic material with crack closure
C effects
          CALL CTDMEL
     1(   DMATX(1:4,1:4),NTYPE      ,RPROPS     ,STRAT      ,STRESK    )
        ELSE
C Error: Material type not recognised
          CALL ERRPRT('EI0044')
        ENDIF
C Perform extra kinematical operations required by materials of this
C class at large strains for computation of the spatial modulus 'a'
C ------------------------------------------------------------------
        IF(NLARGE.EQ.1)THEN
          CALL CSTEP
     1(   AMATX      ,BETRL      ,DMATX      ,STRES     ,DETFT      ,
     2    NTYPE      )
        ENDIF
      ELSEIF(MATCLS.EQ.SINCRY)THEN
C
C Single crystal anisotropic elasto-plastic models
C ================================================
C
C Set plastic loading/unloading flag
        IFPLAS=LALGVA(1)
        IF((.NOT.IFPLAS).OR.KUNLD.EQ.1)THEN
          EPFLAG=.FALSE.
        ELSE
          EPFLAG=.TRUE.
        ENDIF
C Call material type-specific routines
C ------------------------------------
        IF(MATTYP.EQ.PDSCRY)THEN
C Planar double slip single crystal
          CALL CSTPDS
     1(   AMATX(1:5,1:5),RALGVA     ,EPFLAG     ,FINCR      ,IPROPS    ,
     2    LALGVA        ,NTYPE      ,RPROPS     ,RSTAVA     ,RSTAV2    ,
     3    STRES      )
       ELSEIF(MATTYP.EQ.VSCTWD)THEN
C Viscoplastic single crystal model
          CALL CSTVS2
     1(   AMATX         ,DTIME      ,EPFLAG     ,FINCR      ,IPROPS    ,
     2    NTYPE         ,RPROPS     ,RSTAVA     ,RSTAV2     ,STRES     )
         ELSEIF(MATTYP.EQ.MTCTWD)THEN
C Single crystal martensitic transformation visco-plastic model 
          CALL CSTMTC
     1(   AMATX         ,RALGVA     ,EPFLAG     ,FINCR      ,IPROPS    ,
     2    LALGVA        ,NTYPE      ,RPROPS     ,RSTAVA     ,RSTAV2    ,
     3    STRES         ,CTFAIL     ,DTIME      ,FINV       ) 
         ELSEIF(MATTYP.EQ.MTEPTD)THEN
C Single crystal martensitic transformation elasto-plastic model 
          CALL CSTMEP
     1(   AMATX         ,RALGVA     ,EPFLAG     ,FINCR      ,IPROPS    ,
     2    LALGVA        ,NTYPE      ,RPROPS     ,RSTAVA     ,RSTAV2    ,
     3    STRES         ,CTFAIL     ,DTIME      ,FINV       )
        ELSE
C... Error: Material type not recognised
          CALL ERRPRT('EI0044')
        ENDIF
      ELSEIF(MATCLS.EQ.HYPER)THEN
C
C Generic isotropic finite hyperelasticity models
C ===============================================
C
C Call material type-specific routines
C ------------------------------------
        IF(MATTYP.EQ.OGDEN)THEN
C Ogden model
          CALL CSTOGD
     1(   AMATX(1:5,1:5),RSTAVA     ,IPROPS     ,NTYPE      ,RPROPS    ,
     2    STRES      )
        ELSE
C Error: Material type not recognised
          CALL ERRPRT('EI0044')
        ENDIF
      ELSEIF(MATCLS.EQ.PLASTC)THEN
C
C Elasto-plastic materials with small strain implementation only
C ==============================================================
C
        IFPLAS=LALGVA(1)
        IF((.NOT.IFPLAS).OR.KUNLD.EQ.1)THEN
          EPFLAG=.FALSE.
        ELSE
          EPFLAG=.TRUE.
        ENDIF
        IF(MATTYP.EQ.VMMIXD)THEN
C von Mises with mixed isotropic/kinematic hardening
          CALL CTVMMX
     1(   RALGVA     ,DMATX(1:4,1:4),EPFLAG     ,IPROPS     ,NTYPE     ,
     2    RPROPS     ,RSTAVA        ,RSTAV2     ,STRES      )
        ELSEIF(MATTYP.EQ.VVMIXD)THEN
C von Mises viscoplatic with mixed isotropic/kinematic hardening and
C Peric's power law for plastic flow
          CALL CTVVMX
     1(   RALGVA     ,DTIME      ,DMATX(1:4,1:4),EPFLAG     ,IPROPS    ,
     2    NTYPE      ,RPROPS     ,RSTAVA        ,RSTAV2     ,STRES     )
        ELSE
C Error: Material type not recognised
          CALL ERRPRT('EI0044')
        ENDIF
      ELSE
C Error: Material class not recognised
        CALL ERRPRT('EI0043')
      ENDIF
C
      RETURN
      END
CDOC END_SUBROUTINE MATICT
