CDOC BEGIN_SUBROUTINE TUEL
CDOC Thickness update for Hencky elastic model in plane stress
CDOC
CDOC This routine updates the thickness for the (Hencky) elastic model
CDOC under plane stress and large strains. It also computes the total
CDOC deformation gradient (including the thickness strain contribution).
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DETF   <  Determinant of the current total
CDOC C                          deformation gradient (including the
CDOC C                          thickness strain contribution).
CDOC DOUBLE_PRECISION RSTAVA >  Array of current (updated) engineering
CDOC C                          logarithmic strain components.
CDOC DOUBLE_PRECISION THICK  <> Initial (reference) thickness on entry.
CDOC C                          Returns as the current (updated)
CDOC C                          thickness.
CDOC INTEGER MODE            >  Flag. If MODE.NE.1, then only the total
CDOC C                          deformation gradient is computed.
CDOC C                          If MODE = 1, then the thickness is
CDOC C                          updated in addition.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, June 2003: Initial coding
CHST
      SUBROUTINE TUEL
     1(   DETF       ,RSTAVA     ,THICK      ,MODE    )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER( MSTRE=4 )
      DIMENSION
     1    RSTAVA(MSTRE)
C***********************************************************************
C THICKNESS UPDATE FOR THE HENCKY ELASTIC MODEL MODEL UNDER LARGE
C STRAINS AND PLANE STRESS
C
C REFERENCE: Section 13.3.2
C***********************************************************************
C Compute determinant of total deformation gradient (including
C out-of-plane contribution).
C ...start by retrieving the diagonal components of the logarithmic
C    strain tensor
      E11=RSTAVA(1)
      E22=RSTAVA(2)
      E33=RSTAVA(4)
C ...then compute determinant of total deformation gradient
      DETF=EXP(E11+E22+E33)
      IF(MODE.EQ.1)THEN
C Compute thickness stretch
        STRTC3=EXP(E33)
C Update thickness
        THICK=THICK*STRTC3
      ENDIF
C
      RETURN
      END
CDOC END_SUBROUTINE TUEL
