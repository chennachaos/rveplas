CDOC BEGIN_SUBROUTINE TUVM
CDOC Thickness update for large strain von Mises model in plane stress
CDOC
CDOC This routine updates the thickness for the von Mises model under
CDOC plane stress and large strains. It also computes the total
CDOC deformation gradient (including the thickness strain contribution).
CDOC The corresponding state update procedure for the von Mises model
CDOC is carried out in subroutine SUVMPS.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DETF   <> Determinant of the current in-plane
CDOC C                          deformation gradient on entry. Returns
CDOC C                          as the determinant of the current total
CDOC C                          deformation gradient (including the
CDOC C                          thickness strain contribution).
CDOC DOUBLE_PRECISION RSTAVA >  Array of current (updated) real state
CDOC C                          variables other than the stress tensor
CDOC C                          components.
CDOC DOUBLE_PRECISION THICK  <> Initial (reference) thickness on entry.
CDOC C                          Returns as the current (updated)
CDOC C                          thickness.
CDOC INTEGER          MODE   >  Flag. If MODE.NE.1, then only the total
CDOC C                          deformation gradient is computed.
CDOC C                          If MODE = 1, then the thickness is
CDOC C                          updated in addition.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, June 2003: Initial coding
CHST
      SUBROUTINE TUVM
     1(   DETF       ,RSTAVA     ,THICK      ,MODE    )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER( MSTRE=4 )
      DIMENSION
     1    RSTAVA(MSTRE+1)
C***********************************************************************
C THICKNESS UPDATE FOR THE VON MISES ELASTO-PLASTIC MODEL UNDER LARGE
C STRAINS AND PLANE STRESS
C
C REFERENCE: Expressions (14.113-115)
C***********************************************************************
C Compute determinant of total deformation gradient (including
C out-of-plane contribution). Note that, for this model, the determinant
C of the total and elastic deformation gradient coincide due to plastic
C incompressibility.
C... start by retrieving the diagonal components of the elastic
C    logarithmic strain tensor
      EE11=RSTAVA(1)
      EE22=RSTAVA(2)
      EE33=RSTAVA(4)
C... then compute determinant of total deformation gradient
      DETFT=EXP(EE11+EE22+EE33)
      IF(MODE.EQ.1)THEN
C Compute thickness stretch
        STRTC3=DETFT/DETF
C Update thickness
        THICK=THICK*STRTC3
      ENDIF
C return total deformation gradient determinant in DETF
      DETF=DETFT
C
      RETURN
      END
CDOC END_SUBROUTINE TUVM
