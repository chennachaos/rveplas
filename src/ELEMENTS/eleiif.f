CDOC BEGIN_SUBROUTINE ELEIIF
CDOC Element interface for internal force vector calculation
CDOC
CDOC This routine calls the internal force vector calculation routines
CDOC for all element classes available in HYPLAS.
CDOC If the internal force vector calculation fails for some reason
CDOC (such as due to failure of a material-specific state update
CDOC procedure) the return value of the logical argument IFFAIL will be
CDOC .TRUE., which will activate the increment cutting procedure in the
CDOC main program.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DTIME  >  Time increment.
CDOC INTEGER          IELEM  >  Number of the element whose internal
CDOC C                          force vector is to be computed.
CDOC LOGICAL          IFFAIL <  Internal force computation failure flag.
CDOC C                          Return value set to .FALSE. if internal
CDOC C                          force vector was successfully evaluated
CDOC C                          and set to .TRUE. otherwise.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto,  August 1999: Initial coding (split from INTFOR)
CHST
CHST E.de Souza Neto & F.Adziman,
CHST                September 2012: Time increment added
CHST
      SUBROUTINE ELEIIF
     1(   DTIME      ,IELEM      ,IFFAIL   ,TTIME  ,DVOLU  )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE '../MAXDIM.INC'
      INCLUDE '../MATERIAL.INC'
      INCLUDE '../ELEMENTS.INC'
      INCLUDE '../GLBDBASE.INC'
C
      LOGICAL  IFFAIL
C***********************************************************************
C ELEMENT INTERFACE FOR COMPUTATION OF ELEMENT INTERNAL FORCE VECTOR.
C CALL ELEMENT CLASS-SPECIFIC INTERNAL FORCE VECTOR CALCULATION ROUTINES
C
C REFERENCE: Section 5.5.1
C            Figure 5.4
C***********************************************************************
C Initialise internal force calculation failure flag
      IFFAIL=.FALSE.
C Recover element and material type group identification numbers
C --------------------------------------------------------------
      IGRUP =IGRPID(IELEM)
      IELIDN=IELTID(IGRUP)
      MATIDN=MATTID(IGRUP)
C Identify element class
C ----------------------
      IELCLS=IELPRP(2,IELIDN)
C
C Call internal force computation routine according to element class
C ------------------------------------------------------------------
C
      IF(IELCLS.EQ.STDARD)THEN
C Standard 2-D and 3-D displacement-based isoparametric elements
        CALL IFSTD
     1(   IELEM      ,IFFAIL     ,MDIME      ,MELEM      ,MPOIN      ,
     2    MSTRE      ,MTOTV      ,NAXIS      ,NLARGE     ,NTYPE      ,
     3    COORD(1,1,1)       ,DINCR              ,DTIME              ,
     4    ELOAD(1,IELEM)     ,IELPRP(1,IELIDN)   ,IPROPS(1,MATIDN)   ,
     5    LALGVA(1,1,IELEM,1),LNODS              ,RALGVA(1,1,IELEM,1),
     6    RELPRP(1,IELIDN)   ,RPROPS(1,MATIDN)   ,RSTAVA(1,1,IELEM,1),
     7    STRSG(1,1,IELEM,1) ,THKGP(1,IELEM,1)   ,TDISP              ,
     8    TTIME              ,NDIME              ,NDOFN      ,DVOLU  ,
     9    STREPG(1,1,IELEM,1)                                        )
C 2-D F-bar elements (for large strain formulation only)
      ELSEIF(IELCLS.EQ.FBAR)THEN
        CALL IFFBA
     1(   IELEM      ,IFFAIL     ,MDIME      ,MELEM      ,MPOIN      ,
     2    MSTRE      ,MTOTV      ,NAXIS      ,NTYPE      ,
     3    COORD(1,1,1)       ,DINCR              ,DTIME              ,
     4    ELOAD(1,IELEM)     ,IELPRP(1,IELIDN)   ,IPROPS(1,MATIDN)   ,
     5    LALGVA(1,1,IELEM,1),LNODS              ,RALGVA(1,1,IELEM,1),
     6    RELPRP(1,IELIDN)   ,RPROPS(1,MATIDN)   ,RSTAVA(1,1,IELEM,1),
     7    STRSG(1,1,IELEM,1) ,THKGP(1,IELEM,1)   ,TDISP              ,
     8    TTIME              ,NDIME              ,NDOFN              )
      ENDIF
C
      RETURN
      END
CDOC END_SUBROUTINE ELEIIF
