CDOC BEGIN_SUBROUTINE ORMTSC
CDOC Output results for the  martensitic transformation single crystal 
CDOC 2D model.
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
      SUBROUTINE ORMTSC
     1(   NOUTF      ,NTYPE      ,RPROPS     ,RSTAVA     ,STRES      ,
     2    IPROPS     )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER
     1(   IPSYST=8   ,NDIM=2     ,NRSTAV=5   ,NIPROP=4)
C Arguments
      DIMENSION
     1    RPROPS(*)          ,RSTAVA(NRSTAV)     ,STRES(*)           ,
     2    IPROPS(NIPROP) 
C Local arrays and variables
      DIMENSION
     1    AUX1(NDIM,NDIM)    ,DKIRCH(NDIM,NDIM)  ,FE(NDIM,NDIM)      ,
     2    RE(NDIM,NDIM)      ,ROTSTR(NDIM,NDIM)  ,SCHMID(2*IPROPS(4)),
     3    UE(NDIM,NDIM)      ,VECM0(NDIM,2*IPROPS(4))                ,
     4    VECS0(NDIM,2*IPROPS(4))  
      DATA R0   ,R1   ,R3   ,R180   /
     1     0.0D0,1.0D0,3.0D0,180.0D0/
C***********************************************************************
C OUTPUT RESULTS (INTERNAL AND ALGORITHMIC VARIABLES) FOR THE 
C SINGLE CRYSTAL VISCO-PLASTIC MODEL
C***********************************************************************
C Two times of slip systems
      NSYST=2*IPROPS(4)
 1000 FORMAT(' Clrot = ',G12.4, ' Acslip= ',G12.4)      
 1100 FORMAT(' tau',I2,' = ',G12.4)
C Stops program if not plane strain
CC1      IF(NTYPE.NE.2)CALL ERRPRT('EI0065')
      RETURN
CC1
C Retrieve stored state variables
C -------------------------------
C... current elastic deformation gradient
      FE(1,1)=RSTAVA(1)
      FE(2,1)=RSTAVA(2)
      FE(1,2)=RSTAVA(3)
      FE(2,2)=RSTAVA(4)
C... current Taylor hardening variable (acummulated slip)
      HRVAR=RSTAVA(NRSTAV)
C Compute lattice rotation
C ------------------------
C Perform polar decomposition of the elastic deformation gradient
      CALL PODEC2
     1(   FE         ,RE         ,UE         )
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
C Set up initial slip systems vectors
      DO 10 ISYST=1,NSYST*1/2
C... retrieve initial system orientation
C... and convert angles into radians
        RADEG=ACOS(-R1)/R180
        THETA=RADEG*RPROPS(IPSYST+ISYST-1)
C... multi-slip systems vectors
        VECS0(1,ISYST)=COS(THETA)
        VECS0(2,ISYST)=SIN(THETA)
        VECM0(1,ISYST)=-SIN(THETA)
        VECM0(2,ISYST)=COS(THETA)
   10 CONTINUE
C Set up mirrored slip systems vectors
      DO 20 ISYST=NSYST*1/2+1,NSYST
        VECS0(1,ISYST)=-VECS0(1,ISYST-NSYST*1/2)
        VECS0(2,ISYST)=-VECS0(2,ISYST-NSYST*1/2)
        VECM0(1,ISYST)=VECM0(1,ISYST-NSYST*1/2)
        VECM0(2,ISYST)=VECM0(2,ISYST-NSYST*1/2)
   20 CONTINUE
C compute rotated stress tensor
      DETF=FE(1,1)*FE(2,2)-FE(1,2)*FE(2,1)
C... deviatoric Kirchhoff stress
      PKIRCH=R1/R3*DETF*(STRES(1)+STRES(2)+STRES(4))
      DKIRCH(1,1)=DETF*STRES(1)-PKIRCH
      DKIRCH(2,2)=DETF*STRES(2)-PKIRCH
      DKIRCH(1,2)=DETF*STRES(3)
      DKIRCH(2,1)=DKIRCH(1,2)
C... rotated deviatoric Kirchhoff stress
      CALL RVZERO(AUX1,NDIM*NDIM)
      DO 50 I=1,NDIM
        DO 40 J=1,NDIM
          DO 30 K=1,NDIM
            AUX1(I,J)=AUX1(I,J)+RE(K,I)*DKIRCH(K,J)
   30     CONTINUE
   40   CONTINUE
   50 CONTINUE
      CALL RVZERO(ROTSTR,NDIM*NDIM)
      DO 80 I=1,NDIM
        DO 70 J=1,NDIM
          DO 60 K=1,NDIM
            ROTSTR(I,J)=ROTSTR(I,J)+AUX1(I,K)*RE(K,J)
   60     CONTINUE
   70   CONTINUE
   80 CONTINUE
C Current Schmid resolved stresses
      CALL RVZERO(SCHMID,NSYST)
      DO 110 ISYST=1,NSYST
        DO 100 I=1,NDIM
          DO 90 J=1,NDIM
            SCHMID(ISYST)=SCHMID(ISYST)+
     1                    ROTSTR(I,J)*VECS0(I,ISYST)*VECM0(J,ISYST)
   90     CONTINUE                         
  100   CONTINUE        
  110 CONTINUE
C
C Write results to output file
C ----------------------------
      WRITE(NOUTF,1000)CLROT,HRVAR
      DO 120 ISYST=1,NSYST
      WRITE(NOUTF,1100)ISYST,SCHMID(ISYST)
  120 CONTINUE                                
C
      RETURN
      END
CDOC END_SUBROUTINE ORMTSC
