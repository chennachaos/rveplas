CDOC BEGIN_SUBROUTINE MATISW
CDOC Material interface for initialisation and switching state variables
CDOC
CDOC This routine calls the state/algorithmic variables
CDOC initialising/switching routines according to the material type.
CDOC Each material type has its own routine that initialises and
CDOC switches Gauss point data (between current and last converged
CDOC values).
CDOC The initialised/switched data comprises state and algorithmic
CDOC variables whose initialisation/switching rules depend on the
CDOC particular material type considered.
CDOC
CDOC BEGIN_PARAMETERS
CDOC INTEGER          MODE   >  Initialisation/switching mode flag.
CDOC INTEGER          NLARGE >  Large strain flag. Large strain
CDOC C                          analysis if NLARGE=1 and infinitesimal
CDOC C                          strain analysis otherwise.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC INTEGER          IPROPS >  Array of integer material properties.
CDOC LOGICAL          LALGVC <> Array of current logical algorithmic
CDOC C                          variables at the current Gauss point.
CDOC LOGICAL          LALGVL <> Array of logical algorithmic variables
CDOC C                          at the last equilibrium
CDOC C                          configuration for the Gauss point in
CDOC C                          question.
CDOC DOUBLE_PRECISION RALGVC <> Array of current real algorithmic
CDOC C                          variables at the current Gauss point.
CDOC DOUBLE_PRECISION RALGVL <> Array of real algorithmic variables
CDOC C                          at the last equilibrium
CDOC C                          configuration for the current Gauss
CDOC C                          point.
CDOC DOUBLE_PRECISION RPROPS >  Array of real material properties.
CDOC DOUBLE_PRECISION RSTAVC <> Array of current real state variables
CDOC C                          at the current Gauss point.
CDOC DOUBLE_PRECISION RSTAVL <> Array of real state variables at the
CDOC C                          last equilibrium configuration
CDOC C                          at the current Gauss point.
CDOC DOUBLE_PRECISION STRESC <> Array of current (Cauchy) stress
CDOC C                          components at the current Gauss point.
CDOC DOUBLE_PRECISION STRESL <> Array of (Cauchy) stress components at
CDOC C                          last equilibrium configuration for
CDOC C                          the current Gauss point.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, July      1999: Initial coding
CHST
      SUBROUTINE MATISW
     1(   MODE       ,NLARGE     ,NTYPE      ,IPROPS     ,LALGVC     ,
     2    LALGVL     ,RALGVC     ,RALGVL     ,RPROPS     ,RSTAVC     ,
     3    RSTAVL     ,STRESC     ,STRESL     )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE '../MATERIAL.INC'
C Arguments
      LOGICAL   
     1     LALGVC            ,LALGVL
      DIMENSION
     1    IPROPS(*)          ,LALGVC(*)          ,LALGVL(*)          ,
     2    RALGVC(*)          ,RALGVL(*)          ,RPROPS(*)          ,
     3    RSTAVC(*)          ,RSTAVL(*)          ,STRESC(*)          ,
     4    STRESL(*)
C***********************************************************************
C MATERIAL INTERFACE FOR INITIALISATION/SWITCHING ROUTINE CALLS:
C ACCORDING TO THE MATERIAL TYPE, CALLS MATERIAL-SPECIFIC ROUTINE TO
C INITIALISE/SWITCH GAUSS POINT STATE AND ALGORITHMIC VARIABLES
C
C REFERENCE: Sections 5.7.3, 5.7.6
C***********************************************************************
C First identify material type and class
C --------------------------------------
      MATTYP=IPROPS(1)
C
C Then call material type-specific routines
C -----------------------------------------
      IF(MATTYP.EQ.ELASTC)THEN
C Elastic (Hencky material in large strains)
        CALL SWEL
     1(   MODE       ,NTYPE      ,RSTAVC     ,RSTAVL     ,STRESC     ,
     2    STRESL     )
      ELSEIF(MATTYP.EQ.TRESCA)THEN
C Tresca elasto-plastic
        CALL SWTR
     1(   MODE       ,NTYPE      ,LALGVC     ,LALGVL     ,RALGVC     ,
     2    RSTAVC     ,RSTAVL     ,STRESC     ,STRESL     )
      ELSEIF(MATTYP.EQ.VMISES)THEN
C von Mises elasto-plastic
        CALL SWVM
     1(   MODE       ,NTYPE      ,LALGVC     ,LALGVL     ,RALGVC     ,
     2    RSTAVC     ,RSTAVL     ,STRESC     ,STRESL     )
      ELSEIF(MATTYP.EQ.MOHCOU)THEN
C Mohr-Coulomb elasto-plastic
        CALL SWMC
     1(   MODE       ,NTYPE      ,LALGVC     ,LALGVL     ,RALGVC     ,
     2    RSTAVC     ,RSTAVL     ,STRESC     ,STRESL     )
      ELSEIF(MATTYP.EQ.DRUPRA)THEN
C Drucker-Prager elasto-plastic
        CALL SWDP
     1(   MODE       ,NTYPE      ,LALGVC     ,LALGVL     ,RALGVC     ,
     2    RSTAVC     ,RSTAVL     ,STRESC     ,STRESL     )
      ELSEIF(MATTYP.EQ.LEMDAM)THEN
C Lemaitre's ductile damage elasto-plastic model
        CALL SWDAMA
     1(   MODE       ,NTYPE      ,LALGVC     ,LALGVL     ,RALGVC     ,
     2    RSTAVC     ,RSTAVL     ,STRESC     ,STRESL     )
      ELSEIF(MATTYP.EQ.DAMELA)THEN
C Isotropically damaged isotropic elastic material with crack closure
C effects
        CALL SWDMEL
     1(   MODE       ,NTYPE      ,RSTAVC     ,RSTAVL     ,STRESC     ,
     2    STRESL     )
      ELSEIF(MATTYP.EQ.PDSCRY)THEN
C Planar double-slip single crystal
        CALL SWPDSC
     1(   MODE       ,LALGVC     ,LALGVL     ,RALGVC     ,RSTAVC     ,
     2    RSTAVL     ,STRESC     ,STRESL     )
      ELSEIF(MATTYP.EQ.VSCTWD)THEN
C Viscoplastic single crystal model (2D/plane strain and 3D)
        CALL SWVSC2
     1(   MODE       ,LALGVC     ,LALGVL     ,RSTAVC     ,RSTAVL     ,
     2    STRESC     ,STRESL     )
      ELSEIF(MATTYP.EQ.MTCTWD)THEN
C Single crystal martensitic transformation visco-plastic model 
        CALL SWMTSC
     1(   MODE       ,LALGVC     ,LALGVL     ,RALGVC     ,RSTAVC     ,
     2    RSTAVL     ,STRESC     ,STRESL     )
      ELSEIF(MATTYP.EQ.MTEPTD)THEN
C Single crystal martensitic transformation elasto-plastic model
        CALL SWMEPC
     1(   MODE       ,LALGVC     ,LALGVL     ,RALGVC     ,RSTAVC     ,
     2    RSTAVL     ,STRESC     ,STRESL     )
      ELSEIF(MATTYP.EQ.OGDEN)THEN
C Ogden hyperelasticity model
        CALL SWOGD
     1(   MODE       ,NTYPE      ,RSTAVC     ,RSTAVL     ,STRESC     ,
     2    STRESL     )
      ELSEIF(MATTYP.EQ.VMMIXD)THEN
C von Mises with mixed isotropic/kinematic hardening (infinitesimal
C only)
        CALL SWVMMX
     1(   MODE       ,NLARGE     ,NTYPE      ,LALGVC     ,LALGVL     ,
     2    RALGVC     ,RSTAVC     ,RSTAVL     ,STRESC     ,STRESL     )
      ELSEIF(MATTYP.EQ.VVMIXD)THEN
C Viscoplastic von Mises with mixed isotropic/kinematic hardening and
C Peric's power law for viscoplastic flow
        CALL SWVVMX
     1(   MODE       ,NLARGE     ,NTYPE      ,LALGVC     ,LALGVL     ,
     2    RALGVC     ,RSTAVC     ,RSTAVL     ,STRESC     ,STRESL     )
      ELSE
C Error: Material type not recognised
        CALL ERRPRT('EI0046')
      ENDIF
      RETURN
      END
CDOC END_SUBROUTINE MATISW
