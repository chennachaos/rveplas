CDOC BEGIN_SUBROUTINE RDPDSC
CDOC Read data for planar visco-plastic single crystal 2D model.
CDOC
CDOC This routine reads from the data file and echos to the results file
CDOC the material parameters necessary for the large strain
CDOC planar visco-plastic single crystal 2D model with
CDOC piece-wise linear isotropic Taylor hardening.
CDOC It also sets the array of real properties and some components of
CDOC the array of integer material properties.
CDOC These arrays are used by subroutines SUVSC2 and CSTVS2.
CDOC This routine also checks whether the number of integer and real
CDOC material properties, the number of real state variables and the
CDOC number of logical and real algorithmic variables required by the
CDOC present model is compatible with the dimension of the corresponding
CDOC global arrays.
CDOC It also sets the unsymmetric tangent stiffness flag.
CDOC
CDOC BEGIN_PARAMETERS
CDOC INTEGER          IPROPS <  Array of integer material properties.
CDOC INTEGER          MIPROP >  Dimension of the global array of integer
CDOC C                          material properties.
CDOC INTEGER          MLALGV >  Dimension of the global array of logical
CDOC C                          algorithmic variables.
CDOC INTEGER          MRALGV >  Dimension of the global array of real
CDOC C                          algorithmic variables.
CDOC INTEGER          MRPROP >  Dimension of the global array of real
CDOC C                          material variables.
CDOC INTEGER          MRSTAV >  Dimension of the global array of real
CDOC C                          state variables.
CDOC INTEGER          NLARGE >  Large strain analysis integer flag.
CDOC C                          Used for checking only.
CDOC C                          Present model is compatible only with
CDOC C                          large strain analysis.
CDOC INTEGER          NTYPE  >  Stress state type integer flag.
CDOC C                          Used for checking only.
CDOC C                          Present model is compatible only with
CDOC C                          plane strain analysis (NTYPE=2).
CDOC DOUBLE_PRECISION RPROPS <  Array of real material properties.
CDOC LOGICAL          UNSYM  <  Unsymmetric tangent stiffness flag.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, November 1998: Initial coding
CHST E.de Souza Neto,    April 2011: I/O error message added
CHST F. Adziman,       January 2013: Multi-slip systems adopted 
CHST
      SUBROUTINE RDVSC2
     1(   IPROPS     ,MIPROP     ,MLALGV     ,MRALGV     ,MRPROP     ,
     2    MRSTAV     ,NLARGE     ,NTYPE      ,RPROPS     ,UNSYM      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER
     1(   NIPROP=19  ,NLALGV=6   ,NRALGV=0   ,NRSTAV=19   ,MDIM=3
     2     ) !, MSYST=50, IPSYST=8
C Arguments
      LOGICAL UNSYM
      DIMENSION
     1    IPROPS(MIPROP)     ,RPROPS(MRPROP)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: VECS0, VECM0
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: S0M0
      DOUBLE PRECISION, DIMENSION(3,3) :: ROTMAT
      DOUBLE PRECISION, DIMENSION(3) :: EULANG
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: BETASL
C Conversion factor: degrees to radians
      DOUBLE PRECISION, PARAMETER :: DEGRAD=ACOS(-1.0D0)/180.0D0
C***********************************************************************
C READ AND ECHO MATERIAL PROPERTIES FOR THE LARGE STRAIN PLANAR
C SINGLE CRYSTAL VISCO-PLASTIC MODEL WITH NON-LINEAR
C (PIECEWISE LINEAR) ISOTROPIC TAYLOR HARDENING
C***********************************************************************
 1000 FORMAT(' Large strain visco-plastic SINGLE CRYSTAL'/)
 1100 FORMAT(
     1' Mass density ...................................... =',G15.6/
     2' Shear modulus ..................................... =',G15.6/
     3' Bulk modulus ...................................... =',G15.6/
     4' Rate sensitivity parameter                          =',G15.5/
     5' Reference shear strain rate                         =',G15.5)
 1200 FORMAT(/ 
     1' Number of SLIP SYSTEM ............................. =',I3// 
     2' Initial orientation of SLIP SYSTEM'                     /
     3' relative to X-axis (degrees, counterclockwise-positive)'// 
     4'     Slip system no.   Initial orientation'              /)                        
 1300 FORMAT(7X,I3,14X,G15.6)     
 1400 FORMAT(/
     1' Number of points on Taylor hardening curve ........ =',I3//
     2'       Accum. slip    Resolved Schmid yield stress '/)
 1500 FORMAT(2(5X,G15.6))
C
C Set unsymmetric tangent stiffness flag
      UNSYM=.TRUE.
C Check that analysis type and large strain flags are valid for this
C model
      IF((NTYPE/=2).AND.(NTYPE/=4)) STOP 'NTYPE invalid RDVSC2'!CALL ERRPRT('ED0242'), ERRPRT('ED0122')
      IF(NLARGE/=1) STOP 'NLARGE invalid RDVSC2' !CALL ERRPRT('ED0243'), ERRPRT('ED0123')
C
C Read and echo some of the real properties
      WRITE(16,1000)
      READ(15,*,ERR=100,END=100)DENSE
      READ(15,*,ERR=100,END=100)GMODU,BULK
C --------------------------------------------------
C    1. MODEL=1: Norton's model (without yield surfaces)
C    2. MODEL=2: Peric's model (with yield surfaces)
C    3. MODEL=3: Perzyna's model (with yield surfaces)
      READ(15,*,ERR=100,END=100)MODEL,CONST1,CONST2
C!! check values of contants are valid!
      WRITE(16,1100)DENSE,GMODU,BULK,CONST1,CONST2
C real array properties
      RPROPS(1)=DENSE
      RPROPS(2)=GMODU
      RPROPS(3)=BULK
      RPROPS(8)=CONST1
      RPROPS(9)=CONST2
      IPROPS(19)=MODEL
C number of slip systems
      READ(15,*,ERR=100,END=100)NSLSYS
      WRITE(16,1200)NSLSYS
      
      ALLOCATE(VECS0(3,NSLSYS))
      VECS0=R0
      ALLOCATE(VECM0(3,NSLSYS))
      VECM0=R0
      ALLOCATE(S0M0(3,3,NSLSYS))
      S0M0=R0
      
C Set up initial slip systems vectors        
C -----------------------------------
      IF(NTYPE==2)THEN
C Slip system orientation (degrees, counterclockwise-positive)
        READ(15,*,ERR=100,END=100)THETSL
C Angle of slip system shear directions (degrees, 
C counterclockwise-positive)
        ALLOCATE(BETASL(NSLSYS))
        READ(15,*,ERR=100,END=100)BETASL
        DO ISYST=1,NSLSYS
          WRITE(16,1300)ISYST,BETASL(ISYST)
        ENDDO
        
C Convert angles to radians
        THETSL=DEGRAD*THETSL
        BETASL=DEGRAD*BETASL
C System shear direction vectors
        VECS0(1,:)=COS(THETSL+BETASL)
        VECS0(2,:)=SIN(THETSL+BETASL)
C System habit plane normal vectors (m)
        VECM0(1,:)=-SIN(THETSL+BETASL)
        VECM0(2,:)=COS(THETSL+BETASL)
        
        

        
      ELSEIF(NTYPE==4)THEN
C Read euler angles of crystal lattice orientation (Bunge convention)
C (degrees)
        READ(15,*,ERR=100,END=100)(EULANG(I),I=1,3)
        EULANG=EULANG*DEGRAD
C Convert to rotation matrix
        C1=COS(EULANG(1))
        C2=COS(EULANG(2))
        C3=COS(EULANG(3))
        S1=SIN(EULANG(1))
        S2=SIN(EULANG(2))
        S3=SIN(EULANG(3))
        ROTMAT(1,1)=C1*C3-S1*S3*C2
        ROTMAT(1,2)=-C1*S3-S1*C3*C2
        ROTMAT(1,3)=S1*S2
        ROTMAT(2,1)=S1*C3+C1*S3*C2
        ROTMAT(2,2)=-S1*S3+C1*C3*C2
        ROTMAT(2,3)=-C1*S2
        ROTMAT(3,1)=S3*S2
        ROTMAT(3,2)=C3*S2
        ROTMAT(3,3)=C2
C Read plane normal and slip direction vectors (in reference
C crystal orientation)
        DO ISLSYS=1,NSLSYS
          READ(15,*,ERR=100,END=100)
     1                   (VECM0(I,ISLSYS),I=1,3),(VECS0(I,ISLSYS),I=1,3)
C Normalise vectors and rotate to actual crystal coordinate system
          VECM0(:,ISLSYS)=VECM0(:,ISLSYS)/NORM2(VECM0(:,ISLSYS))
          VECS0(:,ISLSYS)=VECS0(:,ISLSYS)/NORM2(VECS0(:,ISLSYS))
          VECM0(:,ISLSYS)=MATMUL(ROTMAT,VECM0(:,ISLSYS))
          VECS0(:,ISLSYS)=MATMUL(ROTMAT,VECS0(:,ISLSYS))
        ENDDO
      ENDIF
CC Read and set slip system
C      DO 10 ISYST=1,NSLSYS
C        READ(15,*,ERR=100,END=100)RPROPS(IPSYST+ISYST-1)
C        WRITE(16,1300)ISYST,RPROPS(IPSYST+ISYST-1)
C   10 CONTINUE
C Set up the constant matrix of slip directions
C ---------------------------------------------
      DO I=1,MDIM
        DO J=1,MDIM
          S0M0(I,J,:)=VECS0(I,:)*VECM0(J,:)
        ENDDO
      ENDDO
C Store slip system vectors in RPROPS
      IPROPS(12)=10
      IPROPS(13)=IPROPS(12)+3*NSLSYS
      RPROPS(IPROPS(12):IPROPS(13))=RESHAPE(VECM0,[3*NSLSYS])
C
      IPROPS(14)=IPROPS(13)+1
      IPROPS(15)=IPROPS(14)+3*NSLSYS
      RPROPS(IPROPS(14):IPROPS(15))=RESHAPE(VECS0,[3*NSLSYS])
C
      IPROPS(16)=IPROPS(15)+1
      IPROPS(17)=IPROPS(16)+3*3*NSLSYS
      RPROPS(IPROPS(16):IPROPS(17))=RESHAPE(S0M0,[3*3*NSLSYS])
C number of points on hardening curve
      READ(15,*,ERR=100,END=100)NHARD
      WRITE(16,1400)NHARD
      IF(NHARD.LT.2) STOP '<2 hardening points' !CALL ERRPRT('ED0128')
      IPROPS(3)=NHARD
      IPROPS(4)=NSLSYS
C size of IPHARD
C      IPHARD=IPSYST+NSLSYS
      IPHARD=IPROPS(17)+1
      IPROPS(18)=IPHARD
C check dimensions of RPROPS
      NRPROP=IPHARD+NHARD*2-1
C!!! need to change above formula
      IF(NRPROP.GT.MRPROP)CALL ERRPRT('ED0131')
C Read and set hardening curve
      DO 20 IHARD=1,NHARD
        READ(15,*,ERR=100,END=100)RPROPS(IPHARD+IHARD*2-2),
     1                            RPROPS(IPHARD+IHARD*2-1)
        WRITE(16,1500)RPROPS(IPHARD+IHARD*2-2),
     1                RPROPS(IPHARD+IHARD*2-1)
   20 CONTINUE
C Check dimension of NSLSYS, RSTAVA, RALGVA, and LALGVA and 
C      IF(NSLSYS.GT.MSYST)STOP 'error reading RDVSC2'!CALL ERRPRT('ED0118')
      IF(NIPROP.GT.MIPROP)STOP 'error reading RDVSC2'!CALL ERRPRT('ED0130')
      IF(NRSTAV.GT.MRSTAV)STOP 'error reading RDVSC2'!CALL ERRPRT('ED0119')
      IF(NRALGV.GT.MRALGV)STOP 'error reading RDVSC2'!CALL ERRPRT('ED0120')
      IF(NLALGV.GT.MLALGV)STOP 'error reading RDVSC2'!CALL ERRPRT('ED0121')
C
      GOTO 200
C Issue error message and abort program execution in case of I/O error
  100 STOP 'error reading RDVSC2'!CALL ERRPRT('ED0124')
C
  200 CONTINUE
      RETURN
      END
CDOC END_SUBROUTINE RPDSC
