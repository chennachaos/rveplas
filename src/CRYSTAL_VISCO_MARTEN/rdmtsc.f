CDOC BEGIN_SUBROUTINE RDMTSC      
      SUBROUTINE RDMTSC
     1(   IPROPS     ,MIPROP     ,MLALGV     ,MRALGV     ,MRPROP     ,
     2    MRSTAV     ,NLARGE     ,NTYPE      ,RPROPS     ,UNSYM      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER
     1(   NIPROP=19   ,NLALGV=6   ,NRALGV=4   ,NRSTAV=52  )
      INTEGER, PARAMETER :: MDIM=3
C Arguments
      LOGICAL UNSYM
      DIMENSION
     1    IPROPS(MIPROP)     ,RPROPS(MRPROP)
     
     
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: BETATR, BETASL
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: HVECS0, HVECM0,
     1            HVECD0
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: VECS0, VECM0
      DOUBLE PRECISION, DIMENSION(3,3) :: ROTMAT
      DOUBLE PRECISION, DIMENSION(3) :: EULANG
      
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: HVED0M, S0M0
C Conversion factor: degrees to radians
      DOUBLE PRECISION, PARAMETER :: DEGRAD=ACOS(-1.0D0)/180.0D0
      
      
C***********************************************************************
C READ AND ECHO MATERIAL PROPERTIES FOR THE LARGE STRAIN PLANAR
C DOUBLE-SLIP SINGLE CRYSTAL ELASTO-PLASTIC MODEL WITH NON-LINEAR
C (PIECEWISE LINEAR) ISOTROPIC TAYLOR HARDENING
C***********************************************************************
 1000 FORMAT(' Martensitic transformation elasto-visco-plastic',
     1       ' SINGLE CRYSTAL (MODIFIED: all systems can transform)'/)
 1100 FORMAT(
     1' Mass density ...................................... =',G15.6/
     2' Shear modulus ..................................... =',G15.6/
     3' Bulk modulus ...................................... =',G15.6/
     4' Habit plane orientation relative to X-axis'                 /
     5' (degrees, counterclockwise-positive) .............. =',G15.6/
     6' Required amount of work to start a transformation . =',G15.6/
     7' Shear magnitude factor ............................ =',G15.6/
     8' Dilatation magnitude factor ....................... =',G15.6)
C
C Set unsymmetric tangent stiffness flag
      UNSYM=.TRUE.
C Check that analysis type and large strain flags are valid for this
C model
      IF((NTYPE/=2).AND.(NTYPE/=4))CALL ERRPRT('ED0242')
      IF(NLARGE/=1) CALL ERRPRT('ED0243')
      
C TODO
C * see notes on RDMEPC
C * update error codes
C * ERROR.RUN - remove any codes related to old CRYSTAL_VISCO_MARTEN model

      

C
C Read and echo some of the real properties
      WRITE(16,1000)
C Mass density
      READ(15,*,ERR=100,END=100)DENSE
C Neo-Hookean hyperelastic properties (bulk, shear moduli)
      READ(15,*,ERR=100,END=100)GMODU,BULK
C Critical energy (DGCR) and number of transformation systems (NTRSYS)
C and (pseudo-)visco transformation properties (Peric model, MU and EPS)
      READ(15,*,ERR=100,END=100)DGCR,NTRSYS,TRMU,TREPS
      ALLOCATE(HVECS0(3,NTRSYS))
      HVECS0=R0
      ALLOCATE(HVECM0(3,NTRSYS))
      HVECM0=R0
      ALLOCATE(HVECD0(3,NTRSYS))
      HVECD0=R0
      ALLOCATE(HVED0M(3,3,NTRSYS))
      HVED0M=R0
C	
C Set up transformation system vectors
C ------------------------------------
      IF(NTYPE==2)THEN
C Transformation system orientation (degrees, counterclockwise-positive)
        READ(15,*,ERR=100,END=100)THETTR
C Angle between two adjacent transformation system habit planes
C (degrees, counterclockwise-positive)
        ALLOCATE(BETATR(NTRSYS))
        READ(15,*,ERR=100,END=100)BETATR
C Convert angles to radians
        THETTR=DEGRAD*THETTR
        BETATR=DEGRAD*BETATR
C System shear direction vectors
        HVECS0(1,:)=COS(THETTR+BETATR)
        HVECS0(2,:)=SIN(THETTR+BETATR)
C System habit plane normal vectors (m)
        HVECM0(1,:)=-SIN(THETTR+BETATR)
        HVECM0(2,:)=COS(THETTR+BETATR)
C
C!!!!!!!!!!!!!!!????
        IF(NTRSYS==4)THEN
          HVECM0(:,2)=-HVECM0(:,2)
          HVECM0(:,3)=-HVECM0(:,3)
        ENDIF
C!!!!!!!!!!!!!!!????
C
C Shear and dilatation magnitude factors (GTRANS, DTRANS)
        READ(15,*,ERR=100,END=100)GTRANS,DTRANS
C System transformation vectors (d)
        HVECD0=GTRANS*HVECS0+DTRANS*HVECM0
C
      ELSEIF(NTYPE==4)THEN
C Read euler angles of crystal lattice orientation (Bunge convention)
C (degrees)
        READ(15,*,ERR=100,END=100)(EULANG(I),I=1,3)
        EULANG=EULANG*DEGRAD ! angles now in radians
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
C Read transformation system vectors in reference (parent lattice)
C coordinate system and rotate them to crystal orientation
        DO ITRSYS=1,NTRSYS
          READ(15,*,ERR=100,END=100)
     1                 (HVECM0(I,ITRSYS),I=1,3),(HVECD0(I,ITRSYS),I=1,3)
C Normalise habit plane normal vector and rotate all vectors to crystal
C orientation
          HVECM0(:,ITRSYS)=HVECM0(:,ITRSYS)/NORM2(HVECM0(:,ITRSYS))
          HVECM0(:,ITRSYS)=MATMUL(ROTMAT,HVECM0(:,ITRSYS))
          HVECD0(:,ITRSYS)=MATMUL(ROTMAT,HVECD0(:,ITRSYS))
        ENDDO
C Calculate shear and dilatation magnitude factors
        DTRANS=DOT_PRODUCT(HVECM0(:,1),HVECD0(:,1))
        GTRANS=SQRT(DOT_PRODUCT(HVECD0(:,1),HVECD0(:,1))-DTRANS**2)
C !!! they are assumed to be the same for all systems!
      ENDIF
C
      WRITE(16,1100)DENSE,GMODU,BULK,THETTR,DGCR,GTRANS,DTRANS
C Compute the tensor product d (x) m for all systems
      DO I=1,MDIM
        DO J=1,MDIM
          HVED0M(I,J,:)=HVECD0(I,:)*HVECM0(J,:)
        ENDDO
      ENDDO
      
      RPROPS(1)=DENSE
      RPROPS(2)=GMODU
      RPROPS(3)=BULK
      RPROPS(5)=DGCR
C
      IPROPS(5)=NTRSYS
C        
C Store transformation system vectors in RPROPS starting at position 10
      IPROPS(6)=10
      IPROPS(7)=IPROPS(6)+3*NTRSYS
      RPROPS(IPROPS(6):IPROPS(7))=RESHAPE(HVECM0,[3*NTRSYS])
C
      IPROPS(8)=IPROPS(7)+1
      IPROPS(9)=IPROPS(8)+3*NTRSYS
      RPROPS(IPROPS(8):IPROPS(9))=RESHAPE(HVECD0,[3*NTRSYS])
C
      IPROPS(10)=IPROPS(9)+1
      IPROPS(11)=IPROPS(10)+3*3*NTRSYS
      RPROPS(IPROPS(10):IPROPS(11))=RESHAPE(HVED0M,[3*3*NTRSYS])
      
      
C Read austenite crystal visco-plasticity properties
C --------------------------------------------------
C    1. MODEL=1 for the Norton's model (without yield surfaces)  
C    2. MODEL=2 for the Peric's model (with yield surfaces) 
      READ(15,*,ERR=100,END=100)MODEL,CONSTM,G0DOT
C      WRITE(16,1100)CONSTM,G0DOT
C number of slip systems
      READ(15,*,ERR=100,END=100)NSLSYS
C      WRITE(16,1200)NSYST
C
      ALLOCATE(VECS0(3,NSLSYS))
      VECS0=R0
      ALLOCATE(VECM0(3,NSLSYS))
      VECM0=R0
      ALLOCATE(S0M0(3,3,NSLSYS))
      S0M0=R0
C
C Set up initial slip systems vectors        
C -----------------------------------
      IF(NTYPE==2)THEN
C Slip system orientation (degrees, counterclockwise-positive)
        READ(15,*,ERR=100,END=100)THETSL
C Angle of slip system shear directions (degrees, 
C counterclockwise-positive)
        ALLOCATE(BETASL(NSLSYS))
        READ(15,*,ERR=100,END=100)BETASL
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
C Read plane normal and slip direction (in parent austenite lattice
C orientation)
        DO ISLSYS=1,NSLSYS
          READ(15,*,ERR=100,END=100)
     1                   (VECM0(I,ISLSYS),I=1,3),(VECS0(I,ISLSYS),I=1,3)
C Normalise and rotate to crystal system
          VECM0(:,ISLSYS)=VECM0(:,ISLSYS)/NORM2(VECM0(:,ISLSYS))
          VECS0(:,ISLSYS)=VECS0(:,ISLSYS)/NORM2(VECS0(:,ISLSYS))
          VECM0(:,ISLSYS)=MATMUL(ROTMAT,VECM0(:,ISLSYS))
          VECS0(:,ISLSYS)=MATMUL(ROTMAT,VECS0(:,ISLSYS))
        ENDDO
      ENDIF
      
C1 Read and set (visco) plastic slip systems
C1        DO ISLYS=1,NSLSYS
C1          READ(15,*,ERR=100,END=100)RPROPS(IPSYST+ISLYS-1)
C1C          WRITE(16,1300)ISLYS,RPROPS(IPSYST+ISLYS-1)
C1        ENDDO
C
C Set up the constant matrix of slip directions
C ---------------------------------------------
      DO I=1,MDIM
        DO J=1,MDIM
          S0M0(I,J,:)=VECS0(I,:)*VECM0(J,:)
        ENDDO
      ENDDO
C
C Store slip system vectors in RPROPS
      IPROPS(12)=IPROPS(11)+1
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
C      
C number of points on hardening curve
      READ(15,*,ERR=100,END=100)NHARD
C      WRITE(16,1400)NHARD
      IF(NHARD.LT.2) STOP '<2 hardening points' !CALL ERRPRT('ED0128')
C size of IPHARD
      IPHARD=IPROPS(17)+1
      IPROPS(18)=IPHARD
C Read and set hardening curve
      DO 20 IHARD=1,NHARD
        READ(15,*,ERR=100,END=100)RPROPS(IPHARD+IHARD*2-2),
     1                            RPROPS(IPHARD+IHARD*2-1)
C        WRITE(16,1500)RPROPS(IPHARD+IHARD*2-2),
C     1                RPROPS(IPHARD+IHARD*2-1)
   20 CONTINUE
C
      IPROPS(19)=MODEL
      RPROPS(8)=CONSTM
      RPROPS(9)=G0DOT
      IPROPS(3)=NHARD
      IPROPS(4)=NSLSYS
      
      
      RPROPS(6)=TRMU
      RPROPS(7)=TREPS
      
C
C Check dimensions of IPROPS
      IF(NIPROP.GT.MIPROP)CALL ERRPRT('ED0130')
C Check dimensions of RPROPS
      NRPROP=IPHARD+NHARD*2-1
      IF(NRPROP.GT.MRPROP)CALL ERRPRT('ED0131')
C Check dimension of RSTAVA, RALGVA and LALGVA
      IF(NRSTAV.GT.MRSTAV)CALL ERRPRT('ED0151')
      IF(NRALGV.GT.MRALGV)CALL ERRPRT('ED0152')
      IF(NLALGV.GT.MLALGV)CALL ERRPRT('ED0153')
C
      GOTO 200
C Issue error message and abort program execution in case of I/O error
  100 CALL ERRPRT('ED0200')
C
  200 CONTINUE
      RETURN
      END
CDOC END_SUBROUTINE RDMTSC