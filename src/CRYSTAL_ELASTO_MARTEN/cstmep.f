CCDOC BEGIN_SUBROUTINE CSTMEP
CDOC Spatial tangent modulus for anisotropic martensitic transformation 
CDOC elasto-plastic-like single crystal model. HYPLAS
CDOC
CDOC This routine computes the (elastic-viscoplastic-transform) spatial
CDOC tangent modulus, a, for martensitic transformation elastoplastic
CDOC single crystal model. The tangent modulus computed here is fully 
CDOC consistent with the elastic predictor/return mapping integration 
CDOC algorithm implemented in subroutine SUMEPC.
CDOC
CDOC This model is restricted to the plane strain case.
CDOC The elastic behaviour of the single crystal is assumed isotropic
CDOC (regularised) Neo-Hookean.
CDOC
CDOC BEGIN_PARAMETERS	
CDOC DOUBLE_PRECISION AMATX  <  Matrix of components of the spatial
CDOC C                          tangent (fourth order) modulus.
CDOC DOUBLE_PRECISION DGAM   >  Array of incremental transformation
CDOC C                          multipliers last determined in
CDOC C                          SUMEPC. All multipliers are set to zero
CDOC C                          for the first iteration of each load
CDOC C                          increment.
CDOC LOGICAL          EPFLAG >  Elasto-transformation flag.
CDOC C                          If .FALSE., AMATX returns as the elastic
CDOC C                          modulus. If .TRUE., AMATX returns as the
CDOC C                          elasto plastic transform modulus
CDOC C                          (step-dependent).      
CDOC DOUBLE_PRECISION FINCR  >  Last incremental deformation gradient.
CDOC INTEGER          IPROPS >  Array of integer material properties.
CDOC LOGICAL          LALGVA >  Array of logical algorithmic flags.
CDOC INTEGER          NTYPE  >  Stress state type flag. The present
CDOC C                          implementation is compatible only with
CDOC C                          plane strain (NTYPE=2).
CDOC DOUBLE_PRECISION RPROPS >  Array of real material properties.
CDOC DOUBLE_PRECISION RSTAVA >  Array of current real state variables
CDOC C                          other than the stress tensor components.
CDOC DOUBLE_PRECISION RSTAVN >  Array of real state variables other than
CDOC C                          the stress tensor components at the
CDOC C                          beginning of the curent load increment.
CDOC DOUBLE_PRECISION STRES  >  Array of current Cauchy stress tensor
CDOC C                          components.
CDOC END_PARAMETERS
CHST
CHST M. Fauzan Adziman, April 2013: Initial coding 
CHST
      SUBROUTINE CSTMEP
     1(   AMATX      ,DGAM       ,EPFLAG     ,FINCR      ,IPROPS     ,
     2    LALGVA     ,NTYPE      ,RPROPS     ,RSTAVA     ,RSTAVN     ,
     3    STRES      ,CTFAIL     ,DTIME      ,FINV       )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER
     1(   MADIM=9    ,NIPROP=19  ,NLALGV=6   ,NRSTAV=52  ,NSTRE=6    )
C Arguments
      LOGICAL
     1    EPFLAG     ,LALGVA(NLALGV)         ,CTFAIL
      DIMENSION
     1    AMATX(MADIM,MADIM) ,FINCR(3,3) ,
     2    IPROPS(NIPROP)     ,RPROPS(*)          ,RSTAVA(NRSTAV)     ,
     3    RSTAVN(NRSTAV)     ,STRES(NSTRE)       ,FINV(3,3) 
C Local arrays and variables
      LOGICAL
     1    IFTRAN     ,IFPLAS  
      DIMENSION
     1    A4TH(3,3,3,3)      ,A4TH2(3,3,3,3)     ,A4TH3(3,3,3,3)     ,
     6    AUXSM1(3,3,3,3)    ,AUXSM2(3,3,3,3)    ,AUXSM3(3,3,3,3)    ,
     8    AUXSM4(3,3,3,3)    ,AUXSM5(3,3)        ,QMATX(3,3,3,3)     ,
     5    DTAUDF(3,3,3,3)    ,DTDFE(3,3,3,3)     ,DFEDF(3,3,3,3)     ,
     1    DFEDDG(3,3)        ,DEVPRJ(9,9)        ,DEVSTR(9)          ,
     4    DPHIA(3,3)         ,DPHIB(3,3)         ,SALPHA(3,3)        ,
     5    SIGMA(3,3)         ,TAUC(3,3)          ,HVD0M(3,3,IPROPS(5))
      DOUBLE PRECISION, DIMENSION(3,3) :: 
     1          FE, FINEL ,FININV, FTOINV, FPAUS, FTRA, FTRINV, FTOTA
      LOGICAL IS2D
      
      
      
      DOUBLE PRECISION AUXSM8(3,3,3,3), FPAINV(3,3)
      
      
      
      DOUBLE PRECISION DTEMP(9)
      
      
      
      DOUBLE PRECISION VECM0(3,IPROPS(4)),VECS0(3,IPROPS(4))
      DOUBLE PRECISION FTOTAW(3,3), FEN(3,3)
      
      
C Local parameters
      DOUBLE PRECISION, PARAMETER :: R0=0.0D0        , RP5=0.5D0       ,
     1                               R1=1.0D0        , R2=2.0D0        ,
     2                               R1D3=1.0D0/3.0D0, R2D3=2.0D0/3.0D0
      DOUBLE PRECISION, PARAMETER, DIMENSION(3,3) :: ! delta_ij
     1     DELTA=reshape((/R1,R0,R0,R0,R1,R0,R0,R0,R1/), (/3,3/) )
C Second-order identity tensor - vector representation using G matrix
C ordering
      DOUBLE PRECISION, PARAMETER, DIMENSION(4) :: SOID=[R1,R0,R0,R1] ! 2D case
      DOUBLE PRECISION, PARAMETER, DIMENSION(9) :: 
     1                            SOID3D=[R1,R0,R0,R0,R1,R0,R0,R0,R1] ! 3D case
C Fourth-order identity tensor (symmetric subspace) - matrix
C representation using G matrix ordering
      DOUBLE PRECISION, PARAMETER, DIMENSION(4,4) :: 
     1 FOIDS=reshape((/R1 ,R0 ,R0 ,R0 , ! 1
     2                 R0 ,RP5,RP5,R0 , ! 2
     3                 R0 ,RP5,RP5,R0 , ! 3
     4                 R0 ,R0 ,R0 ,R1 /), (/4,4/) ) ! 4
C                      1   2   3   4
C 3D G matrix ordering (11,21,31,12,22,32,13,23,33)
C                        1  2  3  4  5  6  7  8  9
C 21,12 -> 2,4
C 31,13 -> 3,7
C 32,23 -> 6,8
      DOUBLE PRECISION, PARAMETER, DIMENSION(9,9) :: 
     1        FOIDS3D=reshape((/R1 ,R0 ,R0 ,R0 ,R0 ,R0 ,R0 ,R0 ,R0 , ! 1
     2                 R0 ,RP5,R0 ,RP5,R0 ,R0 ,R0 ,R0 ,R0 , ! 2
     3                 R0 ,R0 ,RP5,R0 ,R0 ,R0 ,RP5,R0 ,R0 , ! 3
     4                 R0 ,RP5,R0 ,RP5,R0 ,R0 ,R0 ,R0 ,R0 , ! 4
     5                 R0 ,R0 ,R0 ,R0 ,R1 ,R0 ,R0 ,R0 ,R0 , ! 5
     6                 R0 ,R0 ,R0 ,R0 ,R0 ,RP5,R0 ,RP5,R0 , ! 6
     7                 R0 ,R0 ,RP5,R0 ,R0 ,R0 ,RP5,R0 ,R0 , ! 7
     8                 R0 ,R0 ,R0 ,R0 ,R0 ,RP5,R0 ,RP5,R0 , ! 8
     9                 R0 ,R0 ,R0 ,R0 ,R0 ,R0 ,R0 ,R0 ,R1 /), (/9,9/) ) ! 9
C                      1   2   3   4   5   6   7   8   9
C      DOUBLE PRECISION, PARAMETER, DIMENSION(9,9) :: 
CC FOIDS3D(I,J)-R1D3*SOID3D(I)*SOID3D(J)
C     1        DEVPRJ2=[R2D3 ,R0 ,R0 ,R0 ,-R1D3,R0 ,R0 ,R0 ,-R1D3, ! 1
C     2                 R0   ,RP5,R0 ,RP5,R0   ,R0 ,R0 ,R0 ,R0   , ! 2
C     3                 R0   ,R0 ,RP5,R0 ,R0   ,R0 ,RP5,R0 ,R0   , ! 3
C     4                 R0   ,RP5,R0 ,RP5,R0   ,R0 ,R0 ,R0 ,R0   , ! 4
C     5                 -R1D3,R0 ,R0 ,R0 ,R2D3 ,R0 ,R0 ,R0 ,-R1D3, ! 5
C     6                 R0   ,R0 ,R0 ,R0 ,R0   ,RP5,R0 ,RP5,R0   , ! 6
C     7                 R0   ,R0 ,RP5,R0 ,R0   ,R0 ,RP5,R0 ,R0   , ! 7
C     8                 R0   ,R0 ,R0 ,R0 ,R0   ,RP5,R0 ,RP5,R0   , ! 8
C     9                 -R1D3,R0 ,R0 ,R0 ,-R1D3,R0 ,R0 ,R0 ,R2D3 ] ! 9
CC                      1     2   3   4   5     6   7   8   9
C***********************************************************************
C COMPUTATION OF THE CONSISTENT SPATIAL TANGENT MODULUS 'a' FOR
C THE ANISTROPIC TRANSFORMATION SINGLE CRYSTAL ELASTO PLASTIC-LIKE MODEL
C IDEALLY FOR STRESS INDUCED MARTENSITIC TRANSFORMATION.
C MODEL VALID FOR PLANE STRAIN ONLY.
C
C Simplified algorithm of this spatial consistent tangent: 
C --------------------------------------------------------
C   I.  Retrieve material and system properties
C   II. Check transformation consistency flag
C
C A.If EPFLAG=.FALSE. Then
C   Compute elastic tangent modulus a^e
C V.If EPFLAG=.TRUE. Then
C   Compute elastotransformation spatial tangent modulus 
C   B1. Compute DFEDF  (depends on the return-mapping)
C   B2. Compute other tensorial components then assemble DTAUDF
C   B3. Assemble spatial tangent modulus a^ep
C
C***********************************************************************
      
C!!!!!
C!!!!!
C * remove DGAM from argument list
C * remove GOTOs - use IF(ISTRAN), ELSE(ISPLAS), ELSE or something
C * CTFAIL: never set here
C * NGDIM should be 4 or 5 in 2D?
C * FTRINV or FTOINV on calculations -> see F.PhD pages 82-83
C   need to check derivations still stand for elasto-viscoplastic model
C * AMATX dimensions: 4 or 5
C   take into account stress 3,3 in plane strain calculations?
C * Use matching expressions with index notation from write-up
C * change the way IFTRAN and IFPLAS relate to states of transformation:
C   IFTRAN should only be activated when transformation is happening,
C   then .FALSE. once it ends
C   IFPLAS same
C   elastic tangent should be calculated the simple way after transformation ends   
C!!!!!!
C!!!!!!
C!!!!!
C!!!!!
      
C Check analysis type is valid for this model
      IF(NTYPE==2)THEN
        IS2D=.TRUE.
        NDIM=2
        NGDIM=4
      ELSEIF(NTYPE==4)THEN
        IS2D=.FALSE.
        NDIM=3
        NGDIM=9
      ELSE
        CALL ERRPRT('EI0082')
      ENDIF
C Retrieve some state and algorithmic variables
C =============================================
C Retrieve viscoplastic admissibility
      IFTRAN=LALGVA(1)      
      IFPLAS=LALGVA(3)
C Retrieve deformation gradients and compute their inverses:
      FE=RESHAPE(RSTAVA(1:9),[3,3])      ! Elastic
      FPAUS=RESHAPE(RSTAVA(42:50),[3,3]) ! Austenite plastic
      FTRA=RESHAPE(RSTAVA(10:18),[3,3])  ! Transformation
      FTOTA=RESHAPE(RSTAVA(29:37),[3,3]) ! Total
      FINEL=MATMUL(FTRA, FPAUS)          ! Inelastic
      CALL INVMT3(FTRA, FTRINV, DETFTR)
      CALL INVMT3(FTOTA, FTOINV, DETFTO)
      CALL INVMT3(FINEL, FININV, DETFNE)
      
C8
      CALL INVMT3(FPAUS, FPAINV, DETFPA)
C8
      
C??? how can FTRINV be expressed in terms of 
      
C -> check expressions for where it needs to be FININV/FINEL instead of FTRINV/FTRAN
      
C Current value of transformation multiplier
      GAMMA=RSTAVA(19)
C Retrieve material and system properties
C ---------------------------------------
C Neo-Hookean and other constants
      GMODU=RPROPS(2)
      BULK=RPROPS(3)
C Input parameters for the transformation
C ---------------------------------------
C Required amount of work to start a transformation	  
      DGCR=RPROPS(5)
C Retrieve activated variant of the transformation system
      IACTRS=RSTAVA(39)
C Recover systems information from RPROPS
      NTRSYS=IPROPS(5)
      HVD0M=RESHAPE(RPROPS(IPROPS(10):IPROPS(11)),[3,3,NTRSYS])
      
      
C5      IF(IACTRS==0)THEN
C5        FTRA=DELTA
C5      ELSE
C5        FTRA=DELTA+GAMMA*HVD0M(:,:,IACTRS)
C5      ENDIF
C5      FTOTA=MATMUL(MATMUL(FE, FTRAN), FPAUS)
C5      FINEL=MATMUL(FTRA, FPAUS)
C5      CALL INVMT3(FTRA, FTRINV, DETFTR)
C5      CALL INVMT3(FTOTA, FTOINV, DETFTO)
C5      CALL INVMT3(FINEL, FININV, DETFNE)
      
C Retrieve transformation admissibility
      IF(IFTRAN)THEN
C Step is TRANSFROM: Compute ELASTOTRANSFORM spatial 
C                    tangent modulus a^etrans        
C ==================================================
C Compute Kirchoff stress components  
C ----------------------------------	  
        CALL TAUMEP(FE, NDIM, RPROPS, TAUC)
C
C B1. Compute DFEDF, derivative of elastic deformation gradient with
C     respect to total deformation gradient, that depends on the return-
C     mapping algorithm (or simply the algorithmic modulus) 
C ======================================================================
C Compute the Jacobian matrix of the return-mapping
C -------------------------------------------------
C DFEDDG: derivative of the elastic deformation gradient with respect to
C the increment of gamma (transformation multiplier)
C (F.PhD eq. 5.71 - or, equivalently, transpose of tensor M on F.PhD 
C  eq. 5.42)
C ---------------------------------------------------------------      
        DFEDDG=-MATMUL(MATMUL(FE,HVD0M(:,:,IACTRS)),FTRINV)
C DTDFE: derivative of the Kirchoff stress with respect to the elastic
C deformation gradient (only part of jacobian that elastic constitutive 
C model - in this case, regularised Neo-Hookean)
C ---------------------------------------------------------------------
        CALL DTADFE(FE, RPROPS, DTDFE)
C DPHIA: double contraction of DTDFE and DFEDDG
        DO I=1,3
          DO J=1,3
            DPHIA(I,J)=SUM(DTDFE(I,J,:,:)*DFEDDG)
          ENDDO
        ENDDO
C DPHIB: (DFEDDG^T * Tau + Fe^T * DPHIA) * F^-T
C8
C8        DPHIB=MATMUL(
C8     1        MATMUL(TRANSPOSE(DFEDDG),TAUC)+MATMUL(TRANSPOSE(FE),DPHIA)
C8     2        ,TRANSPOSE(FTOINV))
C8
        DPHIB=MATMUL(MATMUL(
     1        MATMUL(TRANSPOSE(DFEDDG),TAUC)+MATMUL(TRANSPOSE(FE),DPHIA)
     2        ,TRANSPOSE(FTOINV)),TRANSPOSE(FPAUS))
C8
C DPHI: double contraction of DPHIB and HVD0M (of active transformation
C system)
        ALPHA=SUM(DPHIB*HVD0M(:,:,IACTRS))
C 
C Compute derivative of transformation multiplier with respect to
C total deformation gradient, denoted as the 2nd order tensor S_alpha
C -------------------------------------------------------------------      
C Compute partial derivative of yield function w.r.t. total deformation
C gradient (DPHIDF)
C... assemble DPHIDF
C      DO I=1,3
C        DO J=1,3
C          DO K=1,3
C            DO L=1,3
C              SUM1=R0
C              SUM2=R0
C              DO M=1,3
C                DO N=1,3
C                  SUM1=SUM1+
C     1               FE(M,I)*SUM(DTDFE(M,N,K,:)*FININV(L,:))*FTOINV(J,N)
C                ENDDO
CC                SUM1=SUM1*FE(M,I)
C                SUM2=SUM2-FE(M,I)*SUM(TAUC(M,:)*FTOINV(L,:))*FTOINV(J,K)
C              ENDDO
C              AUXSM4(I,J,K,L)=FININV(L,I)*SUM( TAUC(K,:)*FTOINV(J,:) )
C     2                        +SUM1+SUM2
C            ENDDO
C          ENDDO
C        ENDDO
C      ENDDO
C DPHIDF_ij = AUXSM4_mnij * HVD0M_mn
          DO M=1,3
            DO N=1,3
              DO K=1,3
                DO L=1,3
                  AUXSM1(M,N,K,L)=SUM( DTDFE(M,N,K,:)*FININV(L,:) ) !mnkl
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          DO M=1,3
            DO J=1,3
              DO K=1,3
                DO L=1,3
                  AUXSM2(M,J,K,L)=SUM( AUXSM1(M,:,K,L)*FTOINV(J,:) ) !mjkl
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          DO M=1,3
            DO J=1,3
              DO K=1,3
                DO L=1,3
                  AUXSM3(M,J,K,L)=-SUM( TAUC(M,:)*FTOINV(L,:))
     1                            *FTOINV(J,K)+AUXSM2(M,J,K,L) !mjkl
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          DO I=1,3
            DO J=1,3
              DO K=1,3
                DO L=1,3
                  AUXSM4(I,J,K,L)=SUM( FE(:,I)*AUXSM3(:,J,K,L) )
     1                           +FININV(L,I)*SUM(TAUC(K,:)*FTOINV(J,:)) !ijkl
                ENDDO
              ENDDO
            ENDDO
          ENDDO
  
          
          
C8
C8
        AUXSM8=AUXSM4
        DO I=1,3
          DO J=1,3
            DO K=1,3
              DO L=1,3
                AUXSM4(I,J,K,L)=SUM(AUXSM8(I,:,K,L)*FPAUS(:,J))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
C8
C8
          
        
        
C Assemble S_ij (F.PhD. eq 5.99)
        DO I=1,3
          DO J=1,3
            SALPHA(I,J)=-SUM(AUXSM4(:,:,I,J)*HVD0M(:,:,IACTRS))/ALPHA
          ENDDO
        ENDDO
C Compute the fourth order tensor Q matrix, or the inelastic part of 
C derivative of elastic deformation gradient with respect to total 
C deformation gradient, Q_ijkl=A_ij S_kl
C ----------------------------------------------------------------
C... compute the 2nd order tensor A_ij (F.PhD. eq 5.91)
        
        
C8
C8        AUXSM5=MATMUL(FE, MATMUL(HVD0M(:,:,IACTRS), FTRINV))
C8
        AUXSM5=MATMUL(FE, MATMUL(FPAINV, 
     1                           MATMUL(HVD0M(:,:,IACTRS), FTRINV)))
C8 FPAINV
        
        
C... finally, assemble Q matrix. Note that Q matrix is only activated 
C    during the transformation, implying that without the Q matrix the
C    consistent tangent is simply returning the elastic component
C    (F.PhD eq 5.90)
        IF(GAMMA < R1)THEN
          DO I=1,3
            DO J=1,3
              DO K=1,3
                DO L=1,3
                  QMATX(I,J,K,L)=AUXSM5(I,J)*SALPHA(K,L) ! Tensor product of AUXSM5 and SALPHA
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ELSE
          QMATX=R0
        ENDIF
C Assemble the algorithmic modulus, DFEDF (F.PhD eq 5.90)
C -------------------------------------------------------
        DO I=1,3
          DO J=1,3
            DO K=1,3
              DO L=1,3
                DFEDF(I,J,K,L)=DELTA(I,K)*FININV(L,J)-QMATX(I,J,K,L)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
C
C B2. Assemble derivative of Kirchoff stress with respect to total 
C     deformation gradient or the DTAUDF part of the consistent spatial 
C     tangent modulus  
C =====================================================================	  	 	      
        DO I=1,3
          DO J=1,3
            DO K=1,3
              DO L=1,3
                DTAUDF(I,J,K,L)=SUM( DTDFE(I,J,:,:)*DFEDF(:,:,K,L) ) ! Double contraction of DTDFE and DFEDF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
C
C B3. Assemble elastotransform spatial tangent modulus (the 'a' matrix)
C =====================================================================	
C Compute NON sigma_il delta_jk part
C ----------------------------------
        DO I=1,3
          DO J=1,3
            DO K=1,3
              DO L=1,3
                A4TH3(I,J,K,L)=SUM( DTAUDF(I,J,K,:)*FTOTA(L,:) )
              ENDDO
            ENDDO
          ENDDO
        ENDDO
C Compute sigma_il delta_jk part
C ------------------------------
        CALL ATSYM( STRES, SIGMA, NDIM, .TRUE.)
C!!!!!! .TRUE. or .FALSE. ??
        DO I=1,3
          DO J=1,3
            DO K=1,3
              DO L=1,3
                A4TH2(I,J,K,L)=SIGMA(I,L)*DELTA(J,K)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
C Finally, assemble elastotransform "a" matrix
C ------------------------------------------
        A4TH=A4TH3/DETFTO-A4TH2
C... store in matrix form (G-matrix component ordering)
        IF(IS2D)THEN
          CALL ARRGAX( A4TH, AMATX(1:5,1:5)  )
        ELSE
          CALL ARRGO3( A4TH, AMATX  )
        ENDIF
C Note: no plasticity is set to occur after transformation
C
      ELSEIF(IFPLAS)THEN  !IF(IFTRAN)
C Step is PLASTIC: Compute ELASTOPLASTIC tangent modulus a^ep
C ===========================================================
        
C5          CALL CSTVCP
C5     1  ( AMATX      ,CTFAIL     ,DTIME      ,FINV       ,FINCR      ,
C5     2    IPROPS     ,NTYPE      ,RPROPS     ,RSTAVA     ,RSTAVN     ,
C5     3    STRES      )
        EPFLAG=.TRUE.
C reorganise RSTAVA
        DTEMP=RSTAVA(10:18)
        RSTAVA(10:18)=RSTAVA(42:50) ! copy FP to its place
        DTEMP2=RSTAVA(19)
        RSTAVA(19)=RSTAVA(41) ! copy HRVAR
        CALL CSTVS2
     1(   AMATX ,  DTIME ,  EPFLAG,  FINCR,  IPROPS,  NTYPE,  RPROPS,
     2    RSTAVA,  RSTAVN,  STRES   )
        RSTAVA(10:18)=DTEMP
        RSTAVA(19)=DTEMP2
        
        
        
      
C          WRITE(*,*)'AMATX HYPLAS ==========================='
C          DO I=1,9
C            WRITE(*,'(9G20.6)')AMATX(I,:)
C          ENDDO
      
C!!!!
C Use Norton model to compare with ELFEM routine
C        IPROPS(19)=1
C!!!!
C          AMATX=R0
C          FEN=RESHAPE(RSTAVN(1:9),[3,3])
C          VECM0=RESHAPE(RPROPS(IPROPS(12):IPROPS(13)),[3,IPROPS(4)])
C          VECS0=RESHAPE(RPROPS(IPROPS(14):IPROPS(15)),[3,IPROPS(4)])
CC ! Wrong total deformation gradient - should be FTOTA from RSTAVA
C          CALL INVMT3(FINV, FTOTAW, DETFIN)
C          IF(IS2D)THEN
C            CALL CTVSC1
C     1(   FTOTAW     ,FINCR      ,FE         ,FEN       ,AMATX(1:4,1:4),
C     2    STRES      ,DUMMY      ,DUMMY      ,VECS0       ,VECM0      ,
C     3    RSTAVA(41) ,RPROPS(2)  ,RPROPS(3)  ,RPROPS(8)   ,RPROPS(9)  ,
C     4    DTIME      ,3          ,4          ,IPROPS(4)   ,4          ,
C     5    IPROPS(3)  ,IDUMMY     ,3      
C     6    ,IPROPS(18), RPROPS)
C1   1(   F          ,FINCR      ,FE         ,FEN         ,AMATX      ,
C1   2    STRESC     ,HEPCUR     ,HSSCUR     ,VECS0       ,VECM0      ,
C1   3    ACSLIP     ,GMODU      ,BULK       ,CONSTM      ,G0DOT      ,
C1   4    DTIME      ,NDIMN      ,NGDIM      ,NSYST       ,NSTREC     ,
C1   5    NHARDS     ,IGROUP     ,NFDIM      
C1   6    ,IPHARD, RPROPS)
C          ELSE
C            CALL CTVSC2
C     1(   FTOTAW     ,FINCR      ,FE         ,FEN         ,AMATX,
C     2    STRES      ,DUMMY      ,DUMMY      ,VECS0       ,VECM0      ,
C     3    RSTAVA(41) ,RPROPS(2)  ,RPROPS(3)  ,RPROPS(8)   ,RPROPS(9)  ,
C     4    DTIME      ,3          ,9          ,IPROPS(4)   ,6          ,
C     5    IPROPS(3)  ,IDUMMY     ,3      
C     6    ,IPROPS(18), RPROPS)
C1   1(   F          ,FINCR      ,FE         ,FEN         ,AMATX      ,
C1   2    STRESC     ,HEPCUR     ,HSSCUR     ,VECS0       ,VECM0      ,
C1   3    ACSLIP     ,GMODU      ,BULK       ,CONSTM      ,G0DOT      ,
C1   4    DTIME      ,NDIMN      ,NGDIM      ,NSYST       ,NSTREC     ,
C1   5    NHARDS     ,IGROUP     ,NFDIM      
C1   6    ,IPHARD, RPROPS)
C          ENDIF
      
C1          WRITE(*,*)'AMATX ELFEM ==========================='
C1          WRITE(*,*)RSTAVA(41)
C1          DO I=1,9
C1            WRITE(*,'(9G20.6)')AMATX(I,:)
C1          ENDDO
          
          !READ(*,*)
      
      
      
      ELSE
C        
C Step is ELASTIC: Compute ELASTIC tangent modulus a^e
C ====================================================
C Get current Cauchy deviatoric stress and Cauchy hydrostatic pressure
        IF(IS2D)THEN
          P=R1D3*(STRES(1)+STRES(2)+STRES(4))
C... use G matrix component ordering to store in-plane deviatoric
C    Cauchy stress components
          DEVSTR(1)=STRES(1)-P
          DEVSTR(2)=STRES(3)
          DEVSTR(3)=STRES(3)
          DEVSTR(4)=STRES(2)-P
        ELSE
          P=R1D3*(STRES(1)+STRES(2)+STRES(3))
          DEVSTR(1)=STRES(1)-P
          DEVSTR(2)=STRES(4)
          DEVSTR(3)=STRES(6)
          DEVSTR(4)=STRES(4)
          DEVSTR(5)=STRES(2)-P
          DEVSTR(6)=STRES(5)
          DEVSTR(7)=STRES(6)
          DEVSTR(8)=STRES(5)
          DEVSTR(9)=STRES(3)-P
        ENDIF
C Trace of isochoric component of Be
C (TRBISO=BEISO_ii=FEISO_ik FEISO_ik = VOLFAC^2 FE_ik FE_ik)
        DETFE=DETM23(3, FE, NDIM, .FALSE.)
        VOLFAC=DETFE**(-R1D3)
        TRBISO=VOLFAC**2 * SUM(FE*FE)
C
        GFAC=R2D3*GMODU*TRBISO/DETFE
        BULFAC=BULK/DETFE
        R2P=R2*P
C
C Assemble deviatoric projection tensor (use G matrix ordering)
C -------------------------------------------------------------
      IF(IS2D)THEN
        DO I=1,NGDIM
          DO J=1,NGDIM
            DEVPRJ(I,J)=FOIDS(I,J)-R1D3*SOID(I)*SOID(J)
          END DO
        END DO
      ELSE
        DO I=1,NGDIM
          DO J=1,NGDIM
            DEVPRJ(I,J)=FOIDS3D(I,J)-R1D3*SOID3D(I)*SOID3D(J)
          END DO
        END DO
      ENDIF
C
C... assemble tensorially compact part (F.PhD eq 5.76)
        IF(IS2D)THEN
          DO I=1,NGDIM
            DO J=1,NGDIM
              AMATX(I,J)=BULFAC*SOID(I)*SOID(J)
     1                   -R2P*FOIDS(I,J)+GFAC*DEVPRJ(I,J)
     2                   -R2D3*(DEVSTR(I)*SOID(J)+SOID(I)*DEVSTR(J))
            END DO
          END DO
        ELSE
          DO I=1,NGDIM
            DO J=1,NGDIM
              AMATX(I,J)=BULFAC*SOID3D(I)*SOID3D(J)
     1                   -R2P*FOIDS3D(I,J)+GFAC*DEVPRJ(I,J)
     2                   -R2D3*(DEVSTR(I)*SOID3D(J)+SOID3D(I)*DEVSTR(J))
            END DO
          END DO
        ENDIF
C
C... add non-compact part: delta_ik sigma_jl (F.PhD eq 5.75)
C AMATX_ijkl = AMATX_ijkl + delta_ik * sigma_jl
C            = AMATX_ijil + sigma_jl
        IF(IS2D)THEN                     !AMATX_ijil - STRES_jl
          AMATX(1,1)=AMATX(1,1)+STRES(1) !AMATX_1111 - STRES_11
          AMATX(1,3)=AMATX(1,3)+STRES(3) !AMATX_1112 - STRES_12
          AMATX(3,1)=AMATX(3,1)+STRES(3) !AMATX_1211 - STRES_21
          AMATX(3,3)=AMATX(3,3)+STRES(2) !AMATX_1212 - STRES_22
          AMATX(2,2)=AMATX(2,2)+STRES(1) !AMATX_2121 - STRES_11
          AMATX(2,4)=AMATX(2,4)+STRES(3) !AMATX_2122 - STRES_12
          AMATX(4,2)=AMATX(4,2)+STRES(3) !AMATX_2221 - STRES_21
          AMATX(4,4)=AMATX(4,4)+STRES(2) !AMATX_2222 - STRES_22
        ELSE
          AMATX(1,1)=AMATX(1,1)+STRES(1) !AMATX_1111 - STRES_11
          AMATX(1,4)=AMATX(1,4)+STRES(4) !AMATX_1112 - STRES_12
          AMATX(1,7)=AMATX(1,7)+STRES(6) !AMATX_1113 - STRES_13
          AMATX(4,1)=AMATX(4,1)+STRES(4) !AMATX_1211 - STRES_21
          AMATX(4,4)=AMATX(4,4)+STRES(2) !AMATX_1212 - STRES_22
          AMATX(4,7)=AMATX(4,7)+STRES(5) !AMATX_1213 - STRES_23
          AMATX(7,1)=AMATX(7,1)+STRES(6) !AMATX_1311 - STRES_31
          AMATX(7,4)=AMATX(7,4)+STRES(5) !AMATX_1312 - STRES_32
          AMATX(7,7)=AMATX(7,7)+STRES(3) !AMATX_1313 - STRES_33
C
          AMATX(2,2)=AMATX(2,2)+STRES(1) !AMATX_2121 - STRES_11
          AMATX(2,5)=AMATX(2,5)+STRES(4) !AMATX_2122 - STRES_12
          AMATX(2,8)=AMATX(2,8)+STRES(6) !AMATX_2123 - STRES_13
          AMATX(5,2)=AMATX(5,2)+STRES(4) !AMATX_2221 - STRES_21
          AMATX(5,5)=AMATX(5,5)+STRES(2) !AMATX_2222 - STRES_22
          AMATX(5,8)=AMATX(5,8)+STRES(5) !AMATX_2223 - STRES_23
          AMATX(8,2)=AMATX(8,2)+STRES(6) !AMATX_2321 - STRES_31
          AMATX(8,5)=AMATX(8,5)+STRES(5) !AMATX_2322 - STRES_32
          AMATX(8,8)=AMATX(8,8)+STRES(3) !AMATX_2323 - STRES_33
C
          AMATX(3,3)=AMATX(3,3)+STRES(1) !AMATX_3131 - STRES_11
          AMATX(3,6)=AMATX(3,6)+STRES(4) !AMATX_3132 - STRES_12
          AMATX(3,9)=AMATX(3,9)+STRES(6) !AMATX_3133 - STRES_13
          AMATX(6,3)=AMATX(6,3)+STRES(4) !AMATX_3231 - STRES_21
          AMATX(6,6)=AMATX(6,6)+STRES(2) !AMATX_3232 - STRES_22
          AMATX(6,9)=AMATX(6,9)+STRES(5) !AMATX_3233 - STRES_23
          AMATX(9,3)=AMATX(9,3)+STRES(6) !AMATX_3331 - STRES_31
          AMATX(9,6)=AMATX(9,6)+STRES(5) !AMATX_3332 - STRES_32
          AMATX(9,9)=AMATX(9,9)+STRES(3) !AMATX_3333 - STRES_33
        ENDIF
C
      ENDIF
C
      RETURN
      END
CDOC END_SUBROUTINE CSTMEP
