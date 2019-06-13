CDOC BEGIN_SUBROUTINE DTADFE
CDOC Derivative of Kirchoff stress with respect to the elastic 
CDOC deformation gradient.
CDOC Current implementation based on a regularised neo-Hookean
CDOC hyperelastic material model.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION FE     >  Elastic deformation gradient.
CDOC DOUBLE_PRECISION RPROPS >  Array of real material properties.
CDOC C                          It is set in routine RDPDSC.
CDOC DOUBLE_PRECISION DTDFE  <  Derivative of the Kirchoff stress with 
CDOC C                          respect to the elastic deformation 
CDOC C                          gradient for the regularised neo-Hookean
CDOC C                          model.
CDOC END_PARAMETERS
CHST
CHST M. Fauzan Adziman, June 2013: Initial coding
CHST
      SUBROUTINE DTADFE( FE, RPROPS, DTDFE )
      IMPLICIT NONE
      INTEGER, PARAMETER :: MDIM=3
C Arguments
      DOUBLE PRECISION, INTENT(IN) , DIMENSION(MDIM,MDIM)     :: FE
      DOUBLE PRECISION, INTENT(IN) , DIMENSION(*)             :: RPROPS
      DOUBLE PRECISION, INTENT(OUT), DIMENSION(MDIM,MDIM,MDIM,MDIM) 
     1                                                        :: DTDFE
C Local variables
      DOUBLE PRECISION :: GMODU, BULK  ! Neo-Hookean shear and bulk moduli
      DOUBLE PRECISION, DIMENSION(MDIM,MDIM) :: FEINV, FEISO ! Inverse and isochoric part of elastic deformation gradient
      DOUBLE PRECISION, DIMENSION(MDIM,MDIM,MDIM,MDIM) :: FTENS, HTENS ! Tensors used in intermediate calculations
C
      INTEGER :: I, J, K, L  ! Loop indices
      DOUBLE PRECISION :: VOLFAC, DETFE
C Local parameters
      DOUBLE PRECISION, PARAMETER :: R0=0.0D0     , R1=1.0D0        ,
     1                               R1D3=R1/3.0D0, R2D3=2.0D0/3.0D0
      DOUBLE PRECISION, PARAMETER, DIMENSION(MDIM,MDIM) :: ! delta_ij
     1     DELTA=reshape((/R1,R0,R0,R0,R1,R0,R0,R0,R1/), (/3,3/) )
C
C***********************************************************************
C    PARTIAL DERIVATIVE OF KIRCHOFF STRESS WITH RESPECT TO THE ELASTIC
C    DEFORMATION GRADIENT FOR A REGULARISED NEO-HOOKEAN HYPERELASTIC 
C    MODEL
C***********************************************************************
C
C Neo-Hookean properties
      GMODU=RPROPS(2)
      BULK=RPROPS(3)
C Inverse of elastic deformation gradient
      CALL INVMT3(FE, FEINV, DETFE)
C Isochoric part of elastic deformation gradient
      VOLFAC=DETFE**(-R1D3)
      FEISO=VOLFAC*FE
C Assemble derivative of the Kirchoff stress with respect to the elastic
C deformation gradient in full 4th order tensor form
C ----------------------------------------------------------------------
C Evaluate tensor F
      DO I=1,MDIM
        DO J=1,MDIM
          DO K=1,MDIM
            DO L=1,MDIM
              FTENS(I,J,K,L)=GMODU*(DELTA(I,K)*FEISO(J,L)
     1                       +FEISO(I,L)*DELTA(J,K)
     2                       -R2D3*DELTA(I,J)*FEISO(K,L))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C Evaluate tensor H
      DO I=1,MDIM
        DO J=1,MDIM
          DO K=1,MDIM
            DO L=1,MDIM
              HTENS(I,J,K,L)=VOLFAC*(DELTA(I,K)*DELTA(J,L)
     1                       -R1D3*FE(I,J)*FEINV(L,K))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C DTDFE = (double contraction of tensors F and H) + tensor P
      DO I=1,MDIM
        DO J=1,MDIM
          DO K=1,MDIM
            DO L=1,MDIM
              DTDFE(I,J,K,L)=SUM(FTENS(I,J,:,:)*HTENS(:,:,K,L))
     1                       +BULK*DELTA(I,J)*FEINV(L,K)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C	  
      RETURN
      END
CDOC END_SUBROUTINE DTADFE
