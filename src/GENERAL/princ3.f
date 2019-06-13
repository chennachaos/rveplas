CDOC BEGIN_SUBROUTINE PRINC3
CDOC Computes the principal stresses in 3-D and outputs them in sorted
CDOC order.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION PSTRS  <  Array of principal stresses in
CDOC                            decreasing order.
CDOC DOUBLE_PRECISION STRSG  >  Array of stress components.
CDOC END_PARAMETERS
CHST
CHST D. de Bortoli, October 2015: Initial coding
CHST
CHST D. de Bortoli,    June 2016: Fixed bugs with repeated principal
CHST                              stresses 
CHST                              
CHST
      SUBROUTINE PRINC3(PSTRS ,STRSG)
      IMPLICIT NONE
C Arguments
      DOUBLE PRECISION, DIMENSION(6)  , INTENT(IN)  :: STRSG
      DOUBLE PRECISION, DIMENSION(3)  , INTENT(OUT) :: PSTRS
C Local variables
      DOUBLE PRECISION, DIMENSION(3,3) :: STRMTX ! Matrix representation of stress tensor
      DOUBLE PRECISION, DIMENSION(3,3) :: STREVC ! Eigenvectors of stress tensor
C
      DOUBLE PRECISION  :: PSTMAX,PSTMIN,PSTMID
      INTEGER           :: I, IMX, IMN, IMD
C Parameters
      INTEGER         , PARAMETER :: NDIME=3
C
C***********************************************************************
C COMPUTES THE PRINCIPAL STRESSES FOR THREE-DIMENSIONAL STRESSES
C***********************************************************************
C Convert array representation of stress tensor to matrix form
      CALL ATSYM(STRSG, STRMTX, NDIME,.FALSE.)
C Find eigenvalues (and eigenvectors) of stress tensor
      CALL JACOB(STRMTX, PSTRS, STREVC, NDIME)
C Sort eigenvalues in decreasing order
C Location of maximum and minimum values in PSTRS
      IMX=MAXLOC(PSTRS,1)
      IMN=MINLOC(PSTRS,1)
C Locate other value
      DO I=1,3
        IF(I==IMX.OR.I==IMN)CYCLE
        IMD=I
        EXIT
      ENDDO
      PSTRS=PSTRS([IMX,IMD,IMN])
C
      RETURN
      END
CDOC END_SUBROUTINE PRINC3
