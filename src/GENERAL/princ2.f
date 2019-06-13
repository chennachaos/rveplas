CDOC BEGIN_SUBROUTINE PRINC2
CDOC Computes the principal stresses and their angle in 2-D.
CDOC
CDOC This routine computes the eigenvalues of the stress tensor in 2-D
CDOC and the angle of the system of eingenvectors relative to the
CDOC system where the stress components are expressed.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION PSTRS  <  Array of principal stresses and angle.
CDOC DOUBLE_PRECISION STRSG  >  Array of stress components.
CDOC END_PARAMETERS
CHST
CHST M.E.Honnor, September 1986
CHST
CHST E.A.de Souza Neto, October 2008: Angle computation changed for
CHST                                  nealy spherical 2-D stress tensor
CHST                                  so that numerically spherical 2-D
CHST                                  stress will render zero angle.
CHST
      SUBROUTINE PRINC2(PSTRS ,STRSG )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PSTRS(3),STRSG(3)
      DATA R0   ,RP01  ,RP5  ,R1   ,R4   ,R90   ,SMALL/
     1     0.0D0,0.01D0,0.5D0,1.0D0,4.0D0,90.0D0,1.D-6/
C***********************************************************************
C COMPUTES THE PRINCIPAL STRESSES FOR TWO-DIMENSIONAL STRESSES
C***********************************************************************
      PI=R4*ATAN(R1)
      XGASH=(STRSG(1)+STRSG(2))*RP5
      XGISH=(STRSG(1)-STRSG(2))*RP5
      XGESH=STRSG(3)
      XGOSH=SQRT(XGISH*XGISH+XGESH*XGESH)
      PSTRS(1)=XGASH+XGOSH
      PSTRS(2)=XGASH-XGOSH
      AUX=SQRT(STRSG(1)**2+STRSG(2)**2+STRSG(3)**2)
      IF(AUX.EQ.R0)AUX=R1
      IF(ABS(XGESH/AUX).LT.SMALL.AND.ABS(XGISH/AUX).LT.SMALL)THEN
        PSTRS(3)=R0
      ELSE
        PSTRS(3)=(ATAN2(XGESH,XGISH))*R90/PI
        IF(PSTRS(3).LT.(-R90+RP01))PSTRS(3)=R90
      ENDIF
      RETURN
      END
CDOC END_SUBROUTINE PRINC2
