CDOC BEGIN_SUBROUTINE MATIRD
CDOC Material interface for reading material-specific input data
CDOC
CDOC This routine calls the routines that read material-specific input
CDOC data from the data file.
CDOC
CDOC BEGIN_PARAMETERS
CDOC CHARACTER        MATNAM >  Character string containing the material
CDOC C                          name.
CDOC INTEGER          NLARGE >  Large strain analysis flag.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC LOGICAL          UNSAUX <  Logical unsymetric stiffness flag.
CDOC INTEGER          IPROPS <  Array of integer material properties.
CDOC DOUBLE_PRECISION RPROPS <  Array of real material properties.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, July 1999
CHST
      SUBROUTINE MATIRD
     1(   MATNAM     ,NLARGE     ,NTYPE      ,UNSAUX     ,IPROPS     ,
     2    RPROPS     )
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C Hyplas database
      INCLUDE '../MATERIAL.INC'
C Arguments
      LOGICAL      UNSAUX
      CHARACTER*80 MATNAM
      DIMENSION
     1    IPROPS(*)          ,RPROPS(*)
C***********************************************************************
C MATERIAL INTERFACE FOR READING MATERIAL-SPECIFIC INPUT DATA
C
C REFERENCE: Sections 5.7.1, 5.7.6
C***********************************************************************
C According to MATNAM, call the appropriate routine to read the
C the material-specific data from the input file and store it
C -------------------------------------------------------------
      IF(MATNAM.EQ.'ELASTIC')THEN
C Elastic
        MATTYP=ELASTC
        MATCLS=HYPEPL
        CALL RDEL
     1(   MRPROP ,MRSTAV ,RPROPS ,UNSAUX)
      ELSEIF(MATNAM.EQ.'TRESCA')THEN
C Tresca elasto-plastic
        MATTYP=TRESCA
        MATCLS=HYPEPL
        CALL RDTR
     1(   IPROPS ,MIPROP ,MLALGV ,MRALGV ,MRPROP ,MRSTAV ,
     2    NLARGE ,NTYPE  ,RPROPS ,UNSAUX)
      ELSEIF(MATNAM.EQ.'VON_MISES')THEN
C von Mises elasto-plastic
        MATTYP=VMISES
        MATCLS=HYPEPL
        CALL RDVM
     1(   IPROPS ,MIPROP ,MLALGV ,MRPROP ,MRSTAV ,
     2    RPROPS ,UNSAUX)
      ELSEIF(MATNAM.EQ.'MOHR_COULOMB')THEN
C Mohr-Coulomb elasto-plastic
        MATTYP=MOHCOU
        MATCLS=HYPEPL
        CALL RDMC
     1(   IPROPS ,MIPROP ,MLALGV ,MRALGV ,MRPROP ,MRSTAV ,
     2    RPROPS ,UNSAUX)
      ELSEIF(MATNAM.EQ.'DRUCKER_PRAGER')THEN
C Drucker-Prager elasto-plastic
        MATTYP=DRUPRA
        MATCLS=HYPEPL
        CALL RDDP
     1(   IPROPS  ,MIPROP ,MLALGV ,MRALGV ,MRPROP ,MRSTAV ,
     2    RPROPS  ,UNSAUX)
      ELSEIF(MATNAM.EQ.'LEMAITRE_DAMAGE')THEN
C Lemaitre's ductile damage model
        MATTYP=LEMDAM
        MATCLS=HYPEPL
        CALL RDDAMA
     1(   IPROPS ,MIPROP ,MLALGV ,MRPROP ,MRSTAV ,NTYPE  ,
     2    RPROPS ,UNSAUX)
      ELSEIF(MATNAM.EQ.'DAMAGED_ELASTIC')THEN
C Isotropically damaged isotropic elastic material with crack closure
C effects
        MATTYP=DAMELA
        MATCLS=HYPEPL
        CALL RDDMEL
     1(   NTYPE ,MRPROP ,MRSTAV ,RPROPS ,UNSAUX)
      ELSEIF(MATNAM.EQ.'OGDEN')THEN
C Ogden hyperelastic
        MATTYP=OGDEN
        MATCLS=HYPER
        CALL RDOGD
     1(   IPROPS ,MIPROP ,MRPROP ,MRSTAV ,RPROPS ,UNSAUX)
      ELSEIF(MATNAM.EQ.'PLANAR_DOUBLE_SLIP_SINGLE_CRYSTAL')THEN
C Planar double-slip single crystal elasto-plastic model 2D
        MATTYP=PDSCRY
        MATCLS=SINCRY
        CALL RDPDSC
     1(   IPROPS ,MIPROP ,MLALGV ,MRALGV ,MRPROP , MRSTAV ,
     2    NLARGE ,NTYPE  ,RPROPS ,UNSAUX )
      ELSEIF(MATNAM.EQ.'VON_MISES_MIXED')THEN	  
C von Mises elasto-plastic with mixed hardening (small strains only)
        MATTYP=VMMIXD
        MATCLS=PLASTC
        CALL RDVMMX
     1(   IPROPS     ,MIPROP     ,MLALGV     ,MRPROP     ,MRSTAV     ,
     2    NLARGE     ,NTYPE      ,RPROPS     ,UNSAUX     )
      ELSEIF(MATNAM.EQ.'VON_MISES_MIXED_VISCO')THEN
C Viscoplastic von Mises elasto-plastic with mixed hardening and Peric's
C power law for viscoplastic flow
        MATTYP=VVMIXD
        MATCLS=PLASTC
        CALL RDVVMX
     1(   IPROPS     ,MIPROP     ,MLALGV     ,MRPROP     ,MRSTAV     ,
     2    NLARGE     ,NTYPE      ,RPROPS     ,UNSAUX     )
      ELSEIF(MATNAM.EQ.'VISCO_SINGLE_CRYSTAL_2D')THEN
C Single crystal viscoplastic model 2D
        MATTYP=VSCTWD
        MATCLS=SINCRY
        CALL RDVSC2
     1(   IPROPS ,MIPROP ,MLALGV ,MRALGV ,MRPROP ,MRSTAV ,
     2    NLARGE ,NTYPE  ,RPROPS ,UNSAUX )
      ELSEIF(MATNAM.EQ.'MARTENSITIC_TRANSFORM_SINGLE_CRYSTAL_2D')THEN
C Single crystal martensitic transformation visco-plastic model 
        MATTYP=MTCTWD
        MATCLS=SINCRY
        CALL RDMTSC
     1(   IPROPS ,MIPROP ,MLALGV ,MRALGV ,MRPROP ,MRSTAV ,
     2    NLARGE ,NTYPE  ,RPROPS ,UNSAUX )
      ELSEIF(MATNAM.EQ.'MARTENSITIC_TRANSFORM_ELASTOPLAS_CRYSTAL')THEN
C Single crystal martensitic transformation elasto-plastic model 
        MATTYP=MTEPTD
        MATCLS=SINCRY
        CALL RDMEPC
     1(   IPROPS ,MIPROP ,MLALGV ,MRALGV ,MRPROP ,MRSTAV ,
     2    NLARGE ,NTYPE  ,RPROPS ,UNSAUX )	 
      ELSE
        CALL ERRPRT('ED0015')
      ENDIF
C Store material type and class flags in IPROPS
C ---------------------------------------------
      IPROPS(1)=MATTYP
      IPROPS(2)=MATCLS
C
      RETURN
      END
CDOC END_SUBROUTINE MATIRD
