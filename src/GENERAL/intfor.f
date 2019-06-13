CDOC BEGIN_SUBROUTINE INTFOR
CDOC Calls internal force vector calculation routines
CDOC
CDOC This routine calls the internal force vector calculation routines
CDOC for all element classes available in HYPLAS. It loops over all
CDOC elements of the mesh. If the internal force vector calculation
CDOC fails for some reason (such as due to failure of an state update
CDOC procedure) the return value of the logical argument INCCUT will be
CDOC .TRUE., which will activate the increment cutting procedure in the
CDOC main program.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DTIME  >  Time increment.
CDOC LOGICAL          INCCUT <  Increment cutting flag. Return value set
CDOC C                          to .FALSE. if internal force vector was
CDOC C                          successfully evaluated and set to .TRUE.
CDOC C                          otherwise.
CDOC C                          When the return value is set to .TRUE.,
CDOC C                          the main program will activate increment
CDOC C                          cutting and the current increment will
CDOC C                          be divided into sub-increments.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto,  June   1996: Initial coding
CHST
CHST E.de Souza Neto,  Jan    1998: Re-organised by element class
CHST
CHST E.de Souza Neto,  August 1999: Element interface call introduced
CHST
CHST E.de Souza Neto & F.Adziman,
CHST                September 2012: Time increment added
CHST
      SUBROUTINE INTFOR( DTIME, INCCUT, TTIME )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE '../MAXDIM.INC'
      INCLUDE '../MATERIAL.INC'
      INCLUDE '../ELEMENTS.INC'
      INCLUDE '../GLBDBASE.INC'
C Arguments
      LOGICAL  INCCUT
C Local variables
      LOGICAL  IFFAIL
C***********************************************************************
C LOOPS OVER ALL ELEMENTS OF THE STRUCTURE TO COMPUTE ELEMENT INTERNAL
C FORCE VECTORS
C
C REFERENCE: Figures 5.2-3
C***********************************************************************
C Initialise increment cutting flag
      INCCUT=.FALSE.
C
C Begin loop over elements
C ========================
      DO 50 IELEM=1,NELEM
C
C Call element interface for internal force vector computation
C ------------------------------------------------------------
        CALL ELEIIF
     1(   DTIME      ,IELEM      ,IFFAIL     ,TTIME  ,DVOLU  )
C
        IF(IFFAIL)THEN
C Internal force calculation failed for current element: Break loop
C over elements and return to main program with increment cutting
C flag activated
          INCCUT=.TRUE.
          GOTO 999
        ENDIF
C
   50 CONTINUE
C Emergency exit
  999 CONTINUE
      RETURN
      END
CDOC END_SUBROUTINE INTFOR
