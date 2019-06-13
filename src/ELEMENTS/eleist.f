CDOC BEGIN_SUBROUTINE ELEIST
CDOC Element interface for element tangent stiffness matrix computation
CDOC
CDOC This routine calls the tangent stiffness matrix calculation
CDOC routines for all element classes available in HYPLAS.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DTIME  >  Time increment.
CDOC DOUBLE_PRECISION ESTIF  <  Tangent stiffness matrix of the current
CDOC C                          element.
CDOC INTEGER          IELEM  >  Number of the element whose stiffness
CDOC C                          matrix is to be computed.
CDOC INTEGER          KUNLD  >  Unloading flag. KUNLD is set to 1 if the
CDOC C                          the loading programme is currently
CDOC C                          unloading.
CDOC LOGICAL          UNSYM  >  Tangent stiffness symmetry flag.
CDOC C                          .TRUE. if tangent stiffness is
CDOC C                          unsymmetric, .FALSE. otherwise.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, June   1996:  initial coding
CHST
CHST E.de Souza Neto, Jan    1998:  re-organised by element class
CHST
CHST E.de Souza Neto, May    1998:  common blocks variables passed to
CHST                                lower level routines via list of
CHST                                arguments
CHST
CHST E.de Souza Neto, August 1999:  name changed (originally was STIFF)
CHST
CHST E.de Souza Neto, October 2008: Unused argument IITER removed
CHST
CHST E.de Souza Neto & F.Adziman,
CHST                September 2012: Time increment added
CHST
      SUBROUTINE ELEIST
     1(   DTIME      ,ESTIF      ,IELEM      ,KUNLD      ,UNSYM      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C Hyplas global database
      INCLUDE '../MAXDIM.INC'
      INCLUDE '../MATERIAL.INC'
      INCLUDE '../ELEMENTS.INC'
      INCLUDE '../GLBDBASE.INC'
C
      LOGICAL UNSYM
      DIMENSION
     1    ESTIF(MEVAB,MEVAB)
C***********************************************************************
C ELEMENT INTERFACE ROUTINE FOR COMPUTATION OF ELEMENT TANGENT STIFFNESS
C MATRIX.
C CALL ELEMENT CLASS-SPECIFIC STIFFNESS MATRIX CALCULATION ROUTINES
C
C REFERENCE: Section 5.6.2
C            Figure 5.5
C***********************************************************************
C
C Recover element and material type group identification numbers
C --------------------------------------------------------------
      IGRUP=IGRPID(IELEM)
      IELIDN=IELTID(IGRUP)
      MATIDN=MATTID(IGRUP)
C Identify element class
      IELCLS=IELPRP(2,IELIDN)
C
C Call stiffness computation routine according to the element class
C -----------------------------------------------------------------
      IF(IELCLS.EQ.STDARD)THEN
        CALL STSTD
     1(   IELEM      ,KUNLD      ,MDIME      ,MELEM      ,
     2    MPOIN      ,MSTRE      ,MTOTV      ,NAXIS      ,NLARGE     ,
     3    NTYPE      ,UNSYM      ,
     4    COORD(1,1,1)       ,DINCR         ,DTIME       ,ESTIF      ,
     5    IELPRP(1,IELIDN)   ,IPROPS(1,MATIDN)   ,LALGVA(1,1,IELEM,1),
     6    LNODS              ,RALGVA(1,1,IELEM,1),RELPRP(1,IELIDN)   ,
     7    RPROPS(1,MATIDN)   ,RSTAVA(1,1,IELEM,1),RSTAVA(1,1,IELEM,2),
     8    STRSG(1,1,IELEM,1) ,THKGP(1,IELEM,1)   ,TDISP              ,
     9    NDIME      ,NDOFN)
      ELSEIF(IELCLS.EQ.FBAR)THEN
        CALL STFBA
     1(   IELEM      ,KUNLD      ,MDIME      ,MELEM      ,
     2    MPOIN      ,MSTRE      ,MTOTV      ,NAXIS      ,
     3    NTYPE      ,UNSYM      ,
     4    COORD(1,1,1)       ,DINCR         ,DTIME       ,ESTIF      ,
     5    IELPRP(1,IELIDN)   ,IPROPS(1,MATIDN)   ,LALGVA(1,1,IELEM,1),
     6    LNODS              ,RALGVA(1,1,IELEM,1),RELPRP(1,IELIDN)   ,
     7    RPROPS(1,MATIDN)   ,RSTAVA(1,1,IELEM,1),RSTAVA(1,1,IELEM,2),
     8    STRSG(1,1,IELEM,1) ,TDISP         ,NDIME       ,NDOFN        )
      ENDIF
C
      RETURN
      END
CDOC END_SUBROUTINE ELEIST
