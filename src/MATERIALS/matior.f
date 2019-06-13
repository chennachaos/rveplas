CDOC BEGIN_SUBROUTINE MATIOR
CDOC Material interface for output result routine calls
CDOC
CDOC This routine calls the material-specific routines to output results
CDOC according to the material type.
CDOC
CDOC BEGIN_PARAMETERS
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC INTEGER          IPROPS >  Array of integer material properties.
CDOC DOUBLE_PRECISION RALGVA >  Array of current real algorithmic
CDOC C                          variables.
CDOC DOUBLE_PRECISION RPROPS >  Array of real material properties.
CDOC DOUBLE_PRECISION RSTAVA >  Array of current real state variables.
CDOC DOUBLE_PRECISION STRES  >  Array of current (Cauchy) stress
CDOC C                          components.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, July      1999: Initial coding
CHST
      SUBROUTINE MATIOR
     1(   NTYPE      ,IPROPS     ,RALGVA     ,RPROPS     ,RSTAVA     ,
     2    STRES      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE '../MATERIAL.INC'
C Arguments
      DIMENSION
     1    IPROPS(*)          ,RALGVA(*)          ,RPROPS(*)          ,
     2    RSTAVA(*)          ,STRES(*)
C***********************************************************************
C MATERIAL INTERFACE FOR OUTPUT RESULT ROUTINE CALLS:
C ACCORDING TO THE MATERIAL TYPE, CALLS MATERIAL-SPECIFIC OUTPUT
C ROUTINE
C
C REFERENCE: Sections 5.7.5-6
C***********************************************************************
C First identify material type
C ----------------------------
      MATTYP=IPROPS(1)
C Then call corresponding routine to output material-specific results
C -------------------------------------------------------------------
      IF(MATTYP.EQ.ELASTC)THEN
C Elastic
        CALL OREL(16      ,NTYPE   ,STRES   )
      ELSEIF(MATTYP.EQ.TRESCA)THEN
C Tresca
        CALL ORTR(RALGVA  ,16      ,NTYPE   ,RSTAVA  ,STRES   )
      ELSEIF(MATTYP.EQ.VMISES)THEN
C von Mises
        CALL ORVM(RALGVA  ,16      ,NTYPE   ,RSTAVA  ,STRES   )
      ELSEIF(MATTYP.EQ.MOHCOU)THEN
C Mohr-Coulomb
        CALL ORMC(RALGVA  ,16      ,NTYPE   ,RSTAVA  ,STRES   )
      ELSEIF(MATTYP.EQ.DRUPRA)THEN
C Drucker-Prager
        CALL ORDP(RALGVA  ,16      ,NTYPE   ,RSTAVA  ,STRES   )
      ELSEIF(MATTYP.EQ.LEMDAM)THEN
C Lemaitre's ductile damage model
        CALL ORDAMA(RALGVA  ,16    ,NTYPE   ,RSTAVA  ,STRES   )
      ELSEIF(MATTYP.EQ.DAMELA)THEN
C Isotropically damaged isotropic elastic material with crack closure
C effects
        CALL ORDMEL(16    ,NTYPE   ,STRES   )
      ELSEIF(MATTYP.EQ.OGDEN)THEN
C Ogden hyperelastic
        CALL OROGD(16     ,NTYPE   ,STRES   )
      ELSEIF(MATTYP.EQ.PDSCRY)THEN
C Planar double-slip single crystal plasticity model
        CALL ORPDSC(RALGVA,16      ,NTYPE   ,RPROPS  ,RSTAVA  ,STRES   )
      ELSEIF(MATTYP.EQ.VSCTWD)THEN
C Single crystal viscoplastic model 2D
        CALL ORVSC2(16      ,NTYPE   ,RPROPS  ,RSTAVA  ,STRES ,IPROPS  )
      ELSEIF(MATTYP.EQ.MTCTWD)THEN
C Single crystal martensitic transformation visco-plastic model
        CALL ORMTSC(16      ,NTYPE   ,RPROPS  ,RSTAVA  ,STRES ,IPROPS  )
      ELSEIF(MATTYP.EQ.MTEPTD)THEN
C Single crystal martensitic transformation elasto-plastic model
        CALL ORMEPC( RALGVA ,16 ,NTYPE ,RPROPS ,RSTAVA ,STRES ,IPROPS )
      ELSEIF(MATTYP.EQ.VMMIXD)THEN
C von Mises with mixed isotropic/kinematic hardening (infinitesimal
C only)
        CALL ORVMMX(RALGVA,16      ,NTYPE   ,RSTAVA  ,STRES   )
      ELSEIF(MATTYP.EQ.VVMIXD)THEN
C Viscoplastic von Mises with mixed isotropic/kinematic hardening
C and Peric's power law for viscoplastic flow
        CALL ORVVMX(RALGVA,16      ,NTYPE   ,RSTAVA  ,STRES   )
      ELSE
C Error: Material type not recognised
        CALL ERRPRT('EI0045')
      ENDIF 
      RETURN
      END
CDOC END_SUBROUTINE MATIOR
