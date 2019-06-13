CDOC BEGIN_SUBROUTINE ORVSC2
CDOC Output results for the visco-plastic single crystal 2D model.
CDOC
CDOC This routine writes to the results file the internal and
CDOC algorithmic variables of the visco-plastic single crystal
CDOC material model with non-linear isotropic Taylor hardening.
CDOC
CDOC BEGIN_PARAMETERS
CDOC INTEGER          NOUTF  >  Results file unit identifier.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC DOUBLE_PRECISION RPROPS >  Array of real material properties.
CDOC INTEGER          IPROPS >  Array of integer material properties.
CDOC C                          The number of points on the piece-wise
CDOC C                          linear hardening curve is the only
CDOC C                          element of this array used here.
CDOC C                          This array is set in routines RVDATA and
CDOC C                          RDVSC2.
CDOC DOUBLE_PRECISION RSTAVA >  Array of real state variables other than
CDOC C                          the stress tensor components.
CDOC DOUBLE_PRECISION STRES  >  Array of stress tensor components.
CDOC 
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, November 1998: Initial coding
CHST F. Adziman, January       2013: Multi-slip systems adopted
CHST
      SUBROUTINE ORVSC2
     1(   NOUTF      ,NTYPE      ,RPROPS     ,RSTAVA     ,STRES      ,
     2    IPROPS     )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER
     1(   NDIM=2     ,NRSTAV=19   ,NIPROP=19, NSTRE=6) !IPSYST=8   ,
C Arguments
      DIMENSION
     1    RPROPS(*)          ,RSTAVA(NRSTAV)     ,STRES(NSTRE)       ,
     2    IPROPS(NIPROP) 
C Local arrays and variables
C      DIMENSION
C     1    AUX1(NDIM,NDIM)    ,DKIRCH(NDIM,NDIM)  ,FE(NDIM,NDIM)      ,
C     2    RE(NDIM,NDIM)      ,ROTSTR(NDIM,NDIM)  ,SCHMID(2*IPROPS(4)),
C     3    UE(NDIM,NDIM)      ,VECM0(NDIM,2*IPROPS(4))                ,
C     4    VECS0(NDIM,2*IPROPS(4))
      DOUBLE PRECISION
     1     VECM0(3,IPROPS(4)) , VECS0(3,IPROPS(4))  , RE(2,2),
     1     VECM(3,IPROPS(4))  , VECS(3,IPROPS(4))   , UE(2,2),
     1     SCHMID(2*IPROPS(4)) , FE(3,3),FEISO(3,3) !S0M0(3,3,IPROPS(4)), 
      DATA R0   ,R1   ,R3   ,R180   /
     1     0.0D0,1.0D0,3.0D0,180.0D0/
C***********************************************************************
C OUTPUT RESULTS (INTERNAL AND ALGORITHMIC VARIABLES) FOR THE 
C SINGLE CRYSTAL VISCO-PLASTIC MODEL
C***********************************************************************
C Two times of slip systems
C      NSYST=2*IPROPS(4)
 1000 FORMAT(' Clrot = ',G12.4, ' Acslip= ',G12.4)      
 1100 FORMAT(' tau',I2,' = ',G12.4)
CC Stops program if not plane strain
C      IF(NTYPE.NE.2)CALL ERRPRT('EI0065')
CC Retrieve stored state variables
CC -------------------------------
CC... current elastic deformation gradient
C      FE(1,1)=RSTAVA(1)
C      FE(2,1)=RSTAVA(2)
C      FE(1,2)=RSTAVA(3)
C      FE(2,2)=RSTAVA(4)
CC... current Taylor hardening variable (acummulated slip)
C      HRVAR=RSTAVA(NRSTAV)
C
      FE=RESHAPE(RSTAVA(1:9),[3,3])
      HRVAR=RSTAVA(19)
C Retrieve material and system properties
C ---------------------------------------
      R1D3=R1/R3
      GMODU=RPROPS(2)
      BULK=RPROPS(3)
C Input parameters for the transformation
C ---------------------------------------
C Recover systems information from RPROPS
      NSLSYS=IPROPS(4)
C      
      VECM0=RESHAPE(RPROPS(IPROPS(12):IPROPS(13)),[3,NSLSYS])
      VECS0=RESHAPE(RPROPS(IPROPS(14):IPROPS(15)),[3,NSLSYS])
C      S0M0=RESHAPE(RPROPS(IPROPS(16):IPROPS(17)),[3,3,NSLSYS])
      
C Compute lattice rotation
C ------------------------
C Perform polar decomposition of the elastic deformation gradient
C      CALL PODEC2
C     1(   FE         ,RE         ,UE         )
      CALL PODEC2
     1(   FE(1:2,1:2)       ,RE         ,UE         )
C From the elastic rotation tensor, compute crystal lattice rotation
      SINE=RE(2,1)
      IF(SINE.GT.R1)SINE=R1
      IF(SINE.LT.-R1)SINE=-R1
      COSINE=RE(1,1)
      IF(COSINE.GT.R1)COSINE=R1
      IF(COSINE.LT.-R1)COSINE=-R1
      DEGRAD=R180/ACOS(-R1)
      SANGLE=DEGRAD*ASIN(SINE)
      CANGLE=DEGRAD*ACOS(COSINE)
      IF(SINE.GE.R0)THEN
        CLROT=CANGLE
      ELSEIF(SINE.LT.R0.AND.COSINE.LT.R0)THEN
        CLROT=-CANGLE
      ELSE
        CLROT=SANGLE
      ENDIF
C Evaluate resolved Schmid stresses (these are not stored in memory)
C ------------------------------------------------------------------
C      DO 10 ISYST=1,NSYST*1/2
CC... retrieve initial system orientation
CC... and convert angles into radians
C        RADEG=ACOS(-R1)/R180
C        THETA=RADEG*RPROPS(IPSYST+ISYST-1)
CC... multi-slip systems vectors
C        VECS0(1,ISYST)=COS(THETA)
C        VECS0(2,ISYST)=SIN(THETA)
C        VECM0(1,ISYST)=-SIN(THETA)
C        VECM0(2,ISYST)=COS(THETA)
C   10 CONTINUE
CC Set up mirrored slip systems vectors
C      DO 20 ISYST=NSYST*1/2+1,NSYST
C        VECS0(1,ISYST)=-VECS0(1,ISYST-NSYST*1/2)
C        VECS0(2,ISYST)=-VECS0(2,ISYST-NSYST*1/2)
C        VECM0(1,ISYST)=VECM0(1,ISYST-NSYST*1/2)
C        VECM0(2,ISYST)=VECM0(2,ISYST-NSYST*1/2)
C   20 CONTINUE
CC compute rotated stress tensor
C      DETF=FE(1,1)*FE(2,2)-FE(1,2)*FE(2,1)
CC... deviatoric Kirchhoff stress
C      PKIRCH=R1/R3*DETF*(STRES(1)+STRES(2)+STRES(4))
C      DKIRCH(1,1)=DETF*STRES(1)-PKIRCH
C      DKIRCH(2,2)=DETF*STRES(2)-PKIRCH
C      DKIRCH(1,2)=DETF*STRES(3)
C      DKIRCH(2,1)=DKIRCH(1,2)
CC... rotated deviatoric Kirchhoff stress
C      CALL RVZERO(AUX1,NDIM*NDIM)
C      DO 50 I=1,NDIM
C        DO 40 J=1,NDIM
C          DO 30 K=1,NDIM
C            AUX1(I,J)=AUX1(I,J)+RE(K,I)*DKIRCH(K,J)
C   30     CONTINUE
C   40   CONTINUE
C   50 CONTINUE
C      CALL RVZERO(ROTSTR,NDIM*NDIM)
C      DO 80 I=1,NDIM
C        DO 70 J=1,NDIM
C          DO 60 K=1,NDIM
C            ROTSTR(I,J)=ROTSTR(I,J)+AUX1(I,K)*RE(K,J)
C   60     CONTINUE
C   70   CONTINUE
C   80 CONTINUE
CC Current Schmid resolved stresses
C      CALL RVZERO(SCHMID,NSYST)
C      DO 110 ISYST=1,NSYST
C        DO 100 I=1,NDIM
C          DO 90 J=1,NDIM
C            SCHMID(ISYST)=SCHMID(ISYST)+
C     1                    ROTSTR(I,J)*VECS0(I,ISYST)*VECM0(J,ISYST)
C   90     CONTINUE                         
C  100   CONTINUE        
C  110 CONTINUE
      NSYST=IPROPS(4)
C Set up initial slip systems vectors
      DETFE=DETM23(3, FE, NDIM, .FALSE.)
      VOLFAC=DETFE**(-R1D3)
      FEISO=VOLFAC*FE
      DO ISYST=1,NSYST
        VECS(:,ISYST)=MATMUL(FEISO, VECS0(:,ISYST))
        VECM(:,ISYST)=MATMUL(FEISO, VECM0(:,ISYST))
        SCHMID(ISYST)=GMODU*DOT_PRODUCT(VECS(:,ISYST),VECM(:,ISYST))
        SCHMID(ISYST+NSYST)=-SCHMID(ISYST)
      ENDDO
C!!!!! store SCHMID as RSTAVA??
C
C Write results to output file
C ----------------------------
      WRITE(NOUTF,1000)CLROT,HRVAR
      DO ISYST=1,2*NSYST
        WRITE(NOUTF,1100)ISYST,SCHMID(ISYST)
      ENDDO
C
      RETURN
      END
CDOC END_SUBROUTINE ORPDSC
