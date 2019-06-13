CDOC BEGIN_SUBROUTINE TAUMEP
CDOC Given the elastic deformation gradient, this subroutine evaluates
CDOC the Kirchhoff stress for the viscoplastic martensitic transforma-
CDOC tion material model.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION FE     >  Elastic deformation gradient.
CDOC INTEGER          NDIME  >  Number of spatial dimension.
CDOC DOUBLE_PRECISION RPROPS >  Array of real material properties.
CDOC C                          It is set in routine RDPDSC.
CDOC DOUBLE_PRECISION TAU    <  Kirchhoff stress tensor.
CDOC END_PARAMETERS
CHST
CHST
CHST
      SUBROUTINE TAUMEP(FE, NDIME, RPROPS, TAU)
      IMPLICIT NONE
C Arguments
      DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN)  :: FE
      INTEGER                         , INTENT(IN)  :: NDIME
      DOUBLE PRECISION, DIMENSION(*)  , INTENT(IN)  :: RPROPS
      DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: TAU
C Local variables
      DOUBLE PRECISION :: GMODU, BULK  ! Neo-Hookean shear and bulk moduli
      DOUBLE PRECISION, DIMENSION(3,3) :: FEISO, ! Isochoric part of elastic deformation gradient
     1                                    BEISO  ! Isochoric elastic left Cauchy-Green tensor
      DOUBLE PRECISION, DIMENSION(6)   :: BEDEV  ! Deviatoric part of BEISO (array form)
C
      DOUBLE PRECISION :: VOLFAC,! Volumetric factor in isochoric split 
     1                    DETFE, ! Determinant of elastic deformation gradient
     2                    TRACE, ! Trace of BEISO
     3                    P      ! Hydrostatic pressure
      INTEGER :: I  ! Loop index
C Functions called
      DOUBLE PRECISION  :: DETM23 ! Determinant of 3x3 matrix
C Local parameters
      DOUBLE PRECISION, PARAMETER :: R1D3=1.0D0/3.0D0
C
C***********************************************************************
C    EVALUATE OF KIRCHHOFF STRESS FOR THE MARTENSITIC TRANSFORMATION 
C    MATERIAL MODEL
C***********************************************************************
C Neo-Hookean properties
      GMODU=RPROPS(2)
      BULK=RPROPS(3)
C Isochoric part of elastic deformation gradient
      DETFE=DETM23(3, FE, NDIME, .FALSE.)
      VOLFAC=DETFE**(-R1D3)
      FEISO=VOLFAC*FE
C Deviatoric part of the isochoric elastic left Cauchy-Green tensor (in
C array form)
      BEISO=MATMUL(FEISO,TRANSPOSE(FEISO))
      TRACE=BEISO(1,1)+BEISO(2,2)+BEISO(3,3)
      CALL SYMTA(BEISO, BEDEV, NDIME, .TRUE.) ! Write BEISO in array form (symmetric)
      BEDEV(1)=BEDEV(1)-R1D3*TRACE
      BEDEV(2)=BEDEV(2)-R1D3*TRACE
      IF(NDIME==2)THEN
        BEDEV(4)=BEDEV(4)-R1D3*TRACE
      ELSEIF(NDIME==3)THEN
        BEDEV(3)=BEDEV(3)-R1D3*TRACE
      ELSE
        CALL ERRPRT('EI0081')
      ENDIF
C Kirchoff stress tensor
      P=BULK*LOG(DETFE)
      CALL ATSYM(BEDEV, TAU, NDIME, .TRUE.) ! Write BEDEV in matrix form (symmetric)
      TAU=GMODU*TAU
      DO I=1,3
        TAU(I,I)=TAU(I,I)+P
      ENDDO
C
      RETURN
      END
CDOC END_SUBROUTINE TAUMEP
