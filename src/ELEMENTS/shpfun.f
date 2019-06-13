CDOC BEGIN_SUBROUTINE SHPFUN
CDOC Call shape function/derivative computation routines
CDOC
CDOC This routine calls the shape function/shape function derivative
CDOC computation routines for all element types. This routine needs to
CDOC be modified for inclusion of new elements in HYPLAS.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DERIV  <  Array of shape function derivatives.
CDOC DOUBLE_PRECISION EISCRD >  Isoparametric coordinates of the
CDOC C                          point where the shape functions or their
CDOC C                          derivatives are to be evaluated, in the
CDOC C                          following order:
CDOC C                            EISCRD(1) -> XI
CDOC C                            EISCRD(2) -> ETA
CDOC C                            EISCRD(3) -> ZETA
CDOC INTEGER          IBOUND >  Boundary interpolation flag. Entry
CDOC C                          must be 0 for domain interpolation.
CDOC C                          Boundary interpolation is assumed
CDOC C                          otherwise.
CDOC INTEGER          IELTYP >  Element type flag. Follows enumeration
CDOC C                          in file ELEMENTS.INC.
CDOC INTEGER          MDIME  >  Maximum permissible number of spatial
CDOC C                          dimensions. Used here only for array
CDOC C                          dimensioning.
CDOC DOUBLE_PRECISION SHAPE  <  Array of shape function values.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, August 1996: Initial coding
CHST   D. de Bortoli,  March 2015: Isoparametric coordinates now grouped
CHST                               in input vector EISCRD
CHST                               Added 3-D elements (linear hexahedra)
CHST
      SUBROUTINE SHPFUN
     1(   DERIV      ,EISCRD      ,IBOUND     ,IELTYP     ,
     2    MDIME      ,SHAPE      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE '../ELEMENTS.INC'
      DIMENSION
     1    DERIV(MDIME,*)     ,SHAPE(*)     ,EISCRD(MDIME)
C***********************************************************************
C CALL SPECIFIC ROUTINES FOR COMPUTATION OF SHAPE FUNCTIONS AND
C SHAPE FUNCTION DERIVATIVES FOR EACH TYPE OF ELEMENT
C
C REFERENCE: Section 5.6.3
C***********************************************************************
C 2-D elements
C ============
      IF(IELTYP.EQ.TRI3)THEN
        ETASP=EISCRD(2)
        EXISP=EISCRD(1)
        CALL SFT3
     1(   DERIV      ,ETASP      ,EXISP      ,IBOUND     ,MDIME      ,
     2    SHAPE      )
      ELSEIF(IELTYP.EQ.TRI6)THEN
        ETASP=EISCRD(2)
        EXISP=EISCRD(1)
        CALL SFT6
     1(   DERIV      ,ETASP      ,EXISP      ,IBOUND     ,MDIME      ,
     2    SHAPE      )
      ELSEIF(IELTYP.EQ.QUAD4)THEN
        ETASP=EISCRD(2)
        EXISP=EISCRD(1)
        CALL SFQ4
     1(   DERIV      ,ETASP      ,EXISP      ,IBOUND     ,MDIME      ,
     2    SHAPE      )
      ELSEIF(IELTYP.EQ.QUAD8)THEN
        ETASP=EISCRD(2)
        EXISP=EISCRD(1)
        CALL SFQ8
     1(   DERIV      ,ETASP      ,EXISP      ,IBOUND     ,MDIME      ,
     2    SHAPE      )
      ELSEIF(IELTYP.EQ.QUA4FB)THEN
        ETASP=EISCRD(2)
        EXISP=EISCRD(1)
        CALL SFQ4FB
     1(   DERIV      ,ETASP      ,EXISP      ,IBOUND     ,MDIME      ,
     2    SHAPE      )
C 3-D elements
C ============
      ELSEIF(IELTYP.EQ.HEXA8)THEN
        ZETASP=EISCRD(3)
        ETASP=EISCRD(2)
        EXISP=EISCRD(1)
        CALL SFH8
     1(   DERIV      ,ZETASP     ,ETASP      ,EXISP      ,IBOUND     ,
     2    MDIME      ,SHAPE      )
      ELSEIF (IELTYP.EQ.TETA7)THEN
       ZETASP=EISCRD(3)
       ETASP=EISCRD(2)
       EXISP=EISCRD(1)
       CALL SFT7E
     1(   DERIV      ,ZETASP     ,ETASP      ,EXISP      ,IBOUND     ,
     2    MDIME      ,SHAPE      )
      ELSEIF(IELTYP.EQ.HEX8FB)THEN
        ZETASP=EISCRD(3)
        ETASP=EISCRD(2)
        EXISP=EISCRD(1)
        CALL SFH8FB
     1(   DERIV      ,ZETASP     ,ETASP      ,EXISP      ,IBOUND     ,
     2    MDIME      ,SHAPE      )
      ELSE
        CALL ERRPRT('EI0005')
      ENDIF
C
      RETURN
      END
CDOC END_SUBROUTINE SHPFUN
