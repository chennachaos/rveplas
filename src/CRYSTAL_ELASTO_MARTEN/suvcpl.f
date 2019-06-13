CDOC BEGIN_SUBROUTINE SUVCPL
CDOC State update for the visco-plastic single crystal 2D model.
CDOC
CDOC This routine uses a fully implicit integration algorithm based 
CDOC on the exponential map approximation of the visco-plastic flow 
CDOC rule as the state update procedure for the anisotropic finite  
CDOC strain visco-plastic single crystal 2D model.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION FINCR  >  Incremental deformation gradient.
CDOC INTEGER          IPROPS >  Array of integer material properties.
CDOC C                          The number of points on the piece-wise
CDOC C                          linear hardening curve is the only
CDOC C                          element of this array used here.
CDOC C                          This array is set in routines RVDATA and
CDOC C                          RDVSC2.
CDOC LOGICAL          LALGVA <  Array of logical algorithmic flags.
CDOC C                          For the present material model, this
CDOC C                          array contains the state update failure
CDOC C                          flag, SUFAIL. The algorithm failure flag 
CDOC C                          is set to .FALSE. if the state update 
CDOC C                          algorithm has been successful and to 
CDOC C                          .TRUE. if the return mapping algorithm 
CDOC C                          has failed for any reason.
CDOC INTEGER          NTYPE  >  Stress state type flag. The present
CDOC C                          implementation is compatible only with
CDOC C                          plane strain (NTYPE=2).
CDOC DOUBLE_PRECISION RPROPS >  Array of real material properties.
CDOC C                          This array contains the density
CDOC C                          (not used in this routine), the elastic
CDOC C                          properties: Shear and bulk modulus
CDOC C                          (neo-Hookean constants) and the plastic
CDOC C                          properties: the pairs "Taylor hardening
CDOC C                          variable-resolved reference stress"
CDOC C                          that define the (user supplied)
CDOC C                          piece-wise linear Taylor hardening
CDOC C                          curve. This array is set in routine
CDOC C                          RDVSC2. In addition, this array 
CDOC C                          incorporated strain-rate sensitivity
CDOC C                          exponent and the reference shear strain
CDOC C                          rate     
CDOC DOUBLE_PRECISION RSTAVA <> Array of real state variables other
CDOC C                          than the stress tensor components.
CDOC C                          Previous converged values on entry,
CDOC C                          updated values on exit.
CDOC C                          The state variables stored in
CDOC C                          this array are the in-plane components
CDOC C                          of the elastic deformation gradient and
CDOC C                          the Taylor hardening internal variable.
CDOC DOUBLE_PRECISION STRES  <  Array of updated stress tensor
CDOC C                          components.
CDOC DOUBLE_PRECISION DTIME  >  Time increment.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, April 1999: Initial coding
CHST G.Q.Liu,         April 1999: Merger to new implicit code in ELFEN     
CHST M.F. Adziman,      Nov 2012: Implemented in HYPLAS v_3.0.0, 
CHST                              a rate-dependent version  
CHST M.F. Adziman,  January 2013: Implemented multi-slip systems and
CHST                              Peric's model (with yield surfaces)
CHST
C***********************************************************************
      SUBROUTINE SUVCPL      
     1(   DTIME      ,FINCR      ,IPROPS     ,LALGVA     ,NTYPE      ,
     2    RPROPS     ,RSTAVA     ,STRES      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER
     1(   MDIM=3     ,MGDIM=9    ,MEQUA=10    ,NIPROP=19   ,IPHVAR=41 ,
     2    NSTRE=6    ,NLALGV=6   ,NRSTAV=52   ,IPSYST=10 )
      LOGICAL
     1    LALGVA(NLALGV), IS2D
      DIMENSION
     1    FINCR(MDIM,MDIM)           ,STRES(NSTRE)                   ,
     2    IPROPS(NIPROP)             ,RPROPS(*)                      ,
     3    RSTAVA(NRSTAV)             
C Local arrays and variables
      LOGICAL
     1    SUFAIL     ,NOCONV     ,SINGUL     ,IFPLAS    
      DIMENSION
     1    ABSSCH(IPROPS(4))          ,SCHMID(IPROPS(4))              ,
     2    OVERST(IPROPS(4))          ,SIGNSC(IPROPS(4))              ,
     3    S0M0(MDIM,MDIM,IPROPS(4))  ,SM0MS0(MDIM,MDIM,IPROPS(4))    ,
     4    VECM(MDIM,IPROPS(4))       ,VECS(MDIM,IPROPS(4))           ,
     5    VECM0(MDIM,IPROPS(4))      ,VECS0(MDIM,IPROPS(4))          ,
     1    BEISO(MDIM,MDIM)           ,DEREXP(MDIM,MDIM,MDIM,MDIM),
     5    DRE12(MDIM,MDIM)           ,DRE21(MDIM,MDIM)               ,
     6    DRE4TH(MDIM,MDIM,MDIM,MDIM),DSUMA(MDIM,MDIM)               ,
     7    DSUMFE(MDIM,MDIM,MDIM,MDIM),A2ND(MDIM,MDIM)                ,
     8    A4TH(MDIM,MDIM,MDIM,MDIM)  ,
     4    DRE11(MGDIM,MGDIM)         ,BEDEV(NSTRE)                   ,
     2    DREMTX(MEQUA,MEQUA)        ,RESVEC(MEQUA)
C
      
      DIMENSION DG(IPROPS(4)), DG2(3,3,IPROPS(4)), DTDFEI(3,3,IPROPS(4))
      LOGICAL ISOVER2(3,3,IPROPS(4))
      DOUBLE PRECISION DGDFEI(3,3,IPROPS(4)), DGDA(IPROPS(4)),
     1      FTIDSM(3,3,IPROPS(4))
      DOUBLE PRECISION DF1DF4(3,3,3,3), DF1DFE(MGDIM,MGDIM), DF1DA(3,3), 
     1                 DF2DFE(3,3), DF2DA
      
      DOUBLE PRECISION, DIMENSION(MDIM, MDIM) ::
     1    FPAUS ,FPAUSN , FTOTA ,FE ,FEINV ,FEISO ,FETISO ,FETFPI ,
     9    FETRL ,FPILOG , FPINCI ,FEN
C
      LOGICAL ISOVER(IPROPS(4))
      DOUBLE PRECISION, PARAMETER :: R0=0.0D0        , R1=1.0D0        ,
     1                               R1D3=1.0D0/3.0D0, SMALL=1.0D-12   ,
     2                               TOL=1.0D-10
      INTEGER, PARAMETER :: MXITER=50
      DOUBLE PRECISION, PARAMETER, DIMENSION(MDIM,MDIM) :: ! delta_ij
     1     DELTA=reshape((/R1,R0,R0,R0,R1,R0,R0,R0,R1/), (/3,3/) )
C***********************************************************************
C STRESS UPDATE PROCEDURE FOR THE ANISOTROPIC PLANAR DOUBLE-SLIP SINGLE
C CRYSTAL VISCO-PLASTIC 2D MODEL:
C FULLY IMPLICIT INTEGRATION ALGORITHM BASED ON THE EXPONENTIAL MAP
C APPROXIMATION OF THE VISCO-PLASTIC FLOW RULE.
C***********************************************************************
C Check analysis type is valid for this model
      IF(NTYPE==2)THEN
        NDIM=2
        NGDIM=4
        NEQUA=5
        IS2D=.TRUE.
      ELSEIF(NTYPE==4)THEN
        NDIM=3
        NGDIM=9
        NEQUA=10
        IS2D=.FALSE.
      ELSE
        STOP 'NTYPE invalid CSTVCP.f' !CALL ERRPRT('EI0006')
      ENDIF
C Initialise state update failure flag
      IFPLAS=.FALSE.
      SUFAIL=.FALSE.
C Set previous converged hardening internal variable
      HRVARN=RSTAVA(IPHVAR)
      HRVAR=HRVARN
C... elastic deformation gradient
      FEN=RESHAPE(RSTAVA(1:9),[MDIM,MDIM])
C... austenite plastic deformation gradient
      FPAUSN=RESHAPE(RSTAVA(42:50),[MDIM,MDIM])
C Retrieve material properties
C... neo-Hookean constants
      GMODU=RPROPS(2)
      BULK=RPROPS(3)  
C... number of sampling points on hardening curve
      NHARD=IPROPS(3)
C... number of slip systems
      NSLSYS=IPROPS(4)
C
      VECM0=RESHAPE(RPROPS(IPROPS(12):IPROPS(13)),[MDIM,NSLSYS])
      VECS0=RESHAPE(RPROPS(IPROPS(14):IPROPS(15)),[MDIM,NSLSYS])
      S0M0=RESHAPE(RPROPS(IPROPS(16):IPROPS(17)),[MDIM,MDIM,NSLSYS])
C
      IPHARD=IPROPS(18)
C
C Stops program if the selected slip-rate law is invalid
C    1. MODEL=1: Norton's creep model (no yield surfaces)
C    2. MODEL=2: Peric model (with yield surfaces)
C    3. MODEL=3: Perzyna model (with yield surfaces)
      MODEL=IPROPS(19)
      IF((MODEL/=1).AND.(MODEL/=2).AND.(MODEL/=3))THEN
        STOP 'ERROR: Slip-rate law chosen is invalid in SUVCPL'
      ENDIF
C Compute deformation gradient components
C ---------------------------------------
C Elastic trial deformation gradient
      FETRL=MATMUL(FINCR, FEN)
C Perform isochoric/volumetric split of elastic trial def. grad.
      DETFET=DETM23(MDIM, FETRL, NDIM, .FALSE.)
      VOLFAC=DETFET**(-R1D3)
      FETISO=VOLFAC*FETRL
C Total deformation gradient
      FTOTA=MATMUL(FETRL, FPAUSN)
C Set up initial guess for visco-plastic solution isochoric elastic
C deformation gradient: use previous (equilibrium) converged as a guess
C ---------------------------------------------------------------------
      DETFEN=DETM23(MDIM, FEN, NDIM, .FALSE.)
      VFACN=DETFEN**(-R1D3)
      FEISO=VFACN*FEN
C!!!      
      IF(IS2D)THEN
        FEISO(3,3)=VOLFAC*FEN(3,3)
      ENDIF
C!!!
      
      
      
      
      
      
      
C      WRITE(*,*)'FEN'
C      DO I=1,3
C        WRITE(*,*)FEN(I,:)
C      ENDDO
C      WRITE(*,*)NHARD,IPHARD
C      WRITE(*,'(2G20.10)')RPROPS(IPHARD:IPHARD+2*NHARD-1)
C      READ(*,*)
      
      
      
      
      
C Check plastic consistency
C -------------------------
C Compute yield functions values
C Compute elastic push forward of all slip-systems vectors
C and current resolved Schmid stresses on all slip systems
      DO ISYST=1,NSLSYS
        VECS(:,ISYST)=MATMUL(FETISO, VECS0(:,ISYST))
        VECM(:,ISYST)=MATMUL(FETISO, VECM0(:,ISYST))
        SCHMID(ISYST)=GMODU*DOT_PRODUCT(VECS(:,ISYST),VECM(:,ISYST))
      ENDDO
C... absolute values
C      ABSSCH=ABS(SCHMID)
      RSHEAR=PLFUN(HRVAR,NHARD,RPROPS(IPHARD))
      OVERST=ABS(SCHMID)/RSHEAR
C9!!! changed for Peric and Perzyna models (yield surface)
C9        ISOVER=OVERST > SMALL
      IF(MODEL==1)THEN ! Norton: no yield surface
        ISOVER=OVERST > SMALL
      ELSE ! Peric, Perzyna
        ISOVER=OVERST-R1 > SMALL
      ENDIF
C9!!! changed for Peric and Perzyna models (yield surface)
C
C Extract viscoplastic properties and set constants
      CONST1=RPROPS(8) !Norton:CONSTM       ; Peric, Perzyna: CONSEP
      R1DC=R1/CONST1
      R1DCM1=R1DC-R1
      CONST2=RPROPS(9) !Norton:G0DOT        ; Peric, Perzyna: CONSMU
      IF(MODEL==1)THEN
        DGAMMA=DTIME*CONST2
      ELSE ! Peric, Perzyna
        DGAMMA=DTIME/CONST2
      ENDIF
C     
C... for the model with yield surfaces, if the current state remains 
C    elastic then return back to SUMEPC for elastic update  
      IF(ANY(ISOVER)) IFPLAS=.TRUE.
C
      IF(.NOT.IFPLAS)GOTO 999
C 
C Start Newton-Raphson iterations for visco-plastic system of equations
C =====================================================================
      DO 360 NRITER=1,MXITER
C Compute residual vector
C -----------------------
C Compute elastic push forward of all slip-systems vectors
C and current resolved Schmid stresses on all slip systems
        DO ISYST=1,NSLSYS
          VECS(:,ISYST)=MATMUL(FEISO, VECS0(:,ISYST))
          VECM(:,ISYST)=MATMUL(FEISO, VECM0(:,ISYST))
          SCHMID(ISYST)=GMODU*DOT_PRODUCT(VECS(:,ISYST),VECM(:,ISYST))
        ENDDO
C... absolute values
C        ABSSCH=ABS(SCHMID)
C Sum up inelastic slip contributions from all systems
C !!!(only systems that yield for Peric's law)
        RSHEAR=PLFUN(HRVAR,NHARD,RPROPS(IPHARD))
        OVERST=ABS(SCHMID)/RSHEAR
C9!!! changed to using SIGN function to avoid possible division by zero
C9        SIGNSC=SCHMID/ABSSCH
        SIGNSC=SIGN(R1, SCHMID)
C9!!! changed to using SIGN function to avoid possible division by zero
        
C9!!! changed for Peric and Perzyna models (yield surface)
C9        ISOVER=OVERST > SMALL
        IF(MODEL==1)THEN ! Norton: no yield surface
          ISOVER=OVERST > SMALL
        ELSE ! Peric, Perzyna
          ISOVER=OVERST-R1 > SMALL
        ENDIF
C9!!! changed for Peric model (yield surface)
        FPILOG=R0
        DO ISYST=1,NSLSYS
          IF(ISOVER(ISYST))THEN
            IF(MODEL==1)THEN
              FACTOR=DGAMMA*SIGNSC(ISYST)*OVERST(ISYST)**R1DC
            ELSEIF(MODEL==2)THEN
              FACTOR=DGAMMA*SIGNSC(ISYST)*(OVERST(ISYST)**R1DC-R1)
            ELSE ! Perzyna
              FACTOR=DGAMMA*SIGNSC(ISYST)*(OVERST(ISYST)-R1)**R1DC
            ENDIF
            FPILOG=FPILOG - FACTOR*S0M0(:,:,ISYST) ! FACTOR = gamma^dot * delta_t = delta_Gamma
          ENDIF !FPILOG = log of inverse (?because of - sign?) of incremental plastic deformation gradient
        ENDDO
        
        
        
        
C        WRITE(*,*)NRITER
C        DO I=1,3
C          WRITE(*,*)FPILOG(I,:)
C        ENDDO
C        READ(*,*)
          
        
        

        
        
C!!!!!
C        HSLOPE=DPLFUN(HRVAR,NHARD,RPROPS(IPHARD))
C        IF(MODEL==1)THEN
C          DTIMEC=DTIME*CONST2
C        ELSE
C          DTIMEC=DTIME/CONST2
C        ENDIF
C        FACTOR=DTIMEC*R1DC*HSLOPE/RSHEAR
C        IF(MODEL==1)THEN
C          DG=DTIMEC*SIGNSC*OVERST**R1DC
C          DGDA=-FACTOR*SIGNSC*OVERST**R1DC
C        ELSEIF(MODEL==2)THEN
C          DG=DTIMEC*SIGNSC*(OVERST**R1DC-R1)
C          DGDA=-FACTOR*SIGNSC*OVERST**R1DC
C        ELSE ! Perzyna
C          DG=DTIMEC*SIGNSC*(OVERST-R1)**R1DC
C          DGDA=-FACTOR*SIGNSC*OVERST*(OVERST-R1)**R1DCM1
C        ENDIF
C        FACTOR=DTIMEC*R1DC/RSHEAR*GMODU
C        DO I=1,3
C          DO J=1,3
C            DGDFEI(I,J,:)=FACTOR*OVERST**R1DCM1*( VECS(I,:)*VECM0(J,:)
C     1                                           +VECM(I,:)*VECS0(J,:) )
C          ENDDO
C        ENDDO
C        
C        DO I=1,MDIM
C          DO J=1,MDIM
C            FPILOG(I,J)=-SUM(DG(:)*S0M0(I,J,:), MASK=ISOVER)
C          ENDDO
C        ENDDO
C        
C        CALL DEXPMP(   DEREXP      ,NOCONV     ,FPILOG     )
C        
CC Fe trial iso * (double contraction of DEREXP and S0M0)
CC SKip calculations for systems where .NOT.ISOVER
C        DO ISYST=1,NSLSYS
C          DO I=1,MDIM
C            DO J=1,MDIM
C              FTIDSM(I,J,ISYST)=SUM( DEREXP(I,J,:,:)*S0M0(:,:,ISYST) )
C            ENDDO
C          ENDDO
C          FTIDSM(:,:,ISYST)=MATMUL(FETISO, FTIDSM(:,:,ISYST))
C        ENDDO
C        
C        DO I=1,MDIM
C          DO J=1,MDIM
C            DO K=1,MDIM
C              DO L=1,MDIM
C                DF1DF4(I,J,K,L)=DELTA(I,K)*DELTA(J,L)
C     1                    +SUM(FTIDSM(I,J,:)*DGDFEI(K,L,:), MASK=ISOVER)
C              ENDDO
C            ENDDO
C          ENDDO
C        ENDDO
C        IF(IS2D)THEN
C          CALL ARRGO2
C     1(   DF1DF4(1:NDIM,1:NDIM,1:NDIM,1:NDIM) ,DF1DFE(1:NGDIM,1:NGDIM) )
C        ELSE
C          CALL ARRGO3
C     1(   DF1DF4 ,DF1DFE )
C        ENDIF
C        
C        DO I=1,MDIM
C          DO J=1,MDIM
C            DF1DA(I,J)=SUM(FTIDSM(I,J,:)*DGDA(:), MASK=ISOVER)
C          ENDDO
C        ENDDO
C        DO I=1,MDIM
C          DO J=1,MDIM
C            DF2DFE(I,J)=-SUM(SIGNSC*DGDFEI(I,J,:), MASK=ISOVER)
C          ENDDO
C        ENDDO
C        
C        DF2DA=R1 - SUM(SIGNSC*DGDA, MASK=ISOVER)
C        
C        DREMTX(1:NGDIM,1:NGDIM)=DF1DFE(1:NGDIM,1:NGDIM)
C        DREMTX(1:NGDIM,NEQUA)=RESHAPE(DF1DA(1:NDIM,1:NDIM), [NGDIM])
C        DREMTX(NEQUA,1:NGDIM)=RESHAPE(DF2DFE(1:NDIM,1:NDIM), [NGDIM])
C        DREMTX(NEQUA,NEQUA)=DF2DA
C!!!!!!!
            
        

        
        
        
        
        
        
C Use exponential map to compute inverse of incremental inelastic
C deformation gradient
        CALL EXPMAP
     1(   FPINCI     ,NOCONV     ,FPILOG     )
        IF(NOCONV)THEN
C Exponential map algorithm failed: Break loop and exit
          SUFAIL=.TRUE.
          WRITE(*,*)'exponential map failed - SUVCPL' !CALL ERRPRT('WE0005')
          GOTO 999
        ENDIF
C Compute residual arranged in vector form (using G matrix component
C ordering)
        FETFPI=MATMUL(FETISO, FPINCI)
C Calculate residual vector
        RESVEC(1:NGDIM)=RESHAPE(
     1            FEISO(1:NDIM,1:NDIM) - FETFPI(1:NDIM,1:NDIM), [NGDIM])
C      
        IF(MODEL==1)THEN
          RSUM=SUM(OVERST**R1DC, MASK=ISOVER)
        ELSEIF(MODEL==2)THEN
          RSUM=SUM(OVERST**R1DC-R1, MASK=ISOVER)
        ELSEIF(MODEL==3)THEN
          RSUM=SUM((OVERST-R1)**R1DC, MASK=ISOVER)
        ENDIF
        RESVEC(NEQUA)=HRVAR-HRVARN-DGAMMA*RSUM
C Check for convergence
C ---------------------
        RESNOR=NORM2(RESVEC(1:NEQUA))
        FACTOR=SQRT( NORM2(FEISO(1:NDIM,1:NDIM))**2+ABS(HRVAR-HRVARN) )
        IF(FACTOR/=R0)RESNOR=RESNOR/FACTOR
        IF(RESNOR<TOL)THEN
C Newton-Raphson iterations converged: Break loop and update stresses
C before exit
          GOTO 370
        ENDIF
C Compute derivatives of residual of isochoric elastic def.grad. update
C ---------------------------------------------------------------------
C Exponential map derivative
        CALL DEXPMP(   DEREXP      ,NOCONV     ,FPILOG     )
C!!!!!!
C!!!!!! CAN PASS SUFAIL directly to DEXPMP and other routines
        IF(NOCONV)THEN
C Derivative of exponential map failed: Break loop and exit
          SUFAIL=.TRUE.
          WRITE(*,*)'Derivative exp map failed - SUVCPL'
          GOTO 999
        ENDIF
C!!!!!!
C Derivative of summation of contributions of each slip plane
        DO I=1,MDIM
          DO J=1,MDIM
            SM0MS0(I,J,:)=VECS(I,:)*VECM0(J,:)+VECM(I,:)*VECS0(J,:)
          ENDDO
        ENDDO
C
        HSLOPE=DPLFUN(HRVAR,NHARD,RPROPS(IPHARD))
        DSUMFE=R0
        DSUMA=R0
        AUX1=DGAMMA*GMODU*R1DC/RSHEAR
        AUX2=DGAMMA*HSLOPE*R1DC/RSHEAR
C
        DO ISYST=1,NSLSYS
          IF(ISOVER(ISYST))THEN
            IF((MODEL==1).OR.(MODEL==2))THEN
              FACTOR=AUX1*OVERST(ISYST)**R1DCM1
              FACTA=AUX2*SIGNSC(ISYST)*OVERST(ISYST)**R1DC
            ELSE
              FACTOR=AUX1*(OVERST(ISYST)-R1)**R1DCM1
              FACTA=AUX2*SIGNSC(ISYST)
     1                  *OVERST(ISYST)*(OVERST(ISYST)-R1)**R1DCM1
            ENDIF
            DSUMA=DSUMA + FACTA*S0M0(:,:,ISYST)
C
            DO I=1,NDIM
              DO J=1,NDIM
                DO K=1,NDIM
                  DO L=1,NDIM
                    DSUMFE(I,J,K,L)=DSUMFE(I,J,K,L)
     1                         +FACTOR*S0M0(I,J,ISYST)*SM0MS0(K,L,ISYST)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDDO
C
        A4TH=R0
        DO I=1,NDIM
          DO J=1,NDIM
            DO K=1,NDIM
              DO L=1,NDIM
                A4TH(I,J,K,L)=SUM( DEREXP(I,J,:,:)*DSUMFE(:,:,K,L) )
              ENDDO
            ENDDO
          ENDDO
        ENDDO
C Assemble (4th order tensor) derivative with respect to isochoric
C elastic def.grad. as a 4-index matrix
        DRE4TH=R0
        DO I=1,NDIM
          DO J=1,NDIM
            DO K=1,NDIM
              DO L=1,NDIM
                DRE4TH(I,J,K,L)=DELTA(I,K)*DELTA(J,L)
     1                          +SUM( FETISO(I,:)*A4TH(:,J,K,L) )
              ENDDO
            ENDDO
          ENDDO
        ENDDO
C Rearrange components as a 2-index matrix (DRE11) using G-matrix
C component ordering
        IF(IS2D)THEN
          CALL ARRGO2
     1(   DRE4TH(1:NDIM,1:NDIM,1:NDIM,1:NDIM) ,DRE11(1:NGDIM,1:NGDIM) )
        ELSE
          CALL ARRGO3
     1(   DRE4TH ,DRE11 )
        ENDIF
C Derivative of residual of isochoric elastic def.grad. update with
C respect to hardening variable
C -----------------------------------------------------------------
        A2ND=R0
        DO I=1,NDIM
          DO J=1,NDIM
            A2ND(I,J)=SUM( DEREXP(I,J,:,:)*DSUMA(:,:) )
          ENDDO
        ENDDO
C
        DRE12=-MATMUL(FETISO, A2ND)
C Derivative of hardening residual with respect to isochoric elastic
C def.grad.
C ------------------------------------------------------------------
        DRE21=R0
        AUX1=DGAMMA*R1DC*GMODU/RSHEAR
        DO ISYST=1,NSLSYS
          IF(ISOVER(ISYST))THEN
            IF((MODEL==1).OR.(MODEL==2))THEN
              FACTOR=-AUX1*SIGNSC(ISYST)*OVERST(ISYST)**R1DCM1
            ELSE
              FACTOR=-AUX1*SIGNSC(ISYST)*(OVERST(ISYST)-R1)**R1DCM1
            ENDIF
            DRE21=DRE21 + FACTOR*SM0MS0(:,:,ISYST)
          ENDIF
        ENDDO
C Derivative of hardening residual with respect to hardening variable
C -------------------------------------------------------------------
        FACTOR=DGAMMA*R1DC*HSLOPE/RSHEAR
        IF((MODEL==1).OR.(MODEL==2))THEN
          DF2DA=R1 + FACTOR*SUM(OVERST**R1DC, MASK=ISOVER)
        ELSE
          DF2DA=R1 + FACTOR*SUM(OVERST*(OVERST-R1)**R1DCM1, MASK=ISOVER)
        ENDIF
C Assemble complete Jacobian matrix
C ---------------------------------
        DREMTX(1:NGDIM,1:NGDIM)=DRE11(1:NGDIM,1:NGDIM)
        DREMTX(1:NGDIM,NEQUA)=RESHAPE(DRE12(1:NDIM,1:NDIM), [NGDIM])
        DREMTX(NEQUA,1:NGDIM)=RESHAPE(DRE21(1:NDIM,1:NDIM), [NGDIM])
        DREMTX(NEQUA,NEQUA)=DF2DA
C Solve linear system and apply Newton-Raphson correction to
C isochoric elastic deformation gradient and hardening variable
C -------------------------------------------------------------
C Solve linear algebraic system
        CALL GAUSEL
     1(   DREMTX(1:NEQUA,1:NEQUA) ,RESVEC(1:NEQUA) ,NEQUA  ,SINGUL  )
        IF(SINGUL)THEN
C Linear system is singular: switch state update failure flag on, break 
C Newton-Raphson loop and exit without updating stresses
          SUFAIL=.TRUE.
          WRITE(*,*)'GAUSEL failed! SUVCPL'!CALL ERRPRT('WE0007')
          GOTO 999
        ENDIF
C Apply Newton correction to isochoric elastic deformation gradient
        FEISO(1:NDIM,1:NDIM)=FEISO(1:NDIM,1:NDIM)
     1                       -RESHAPE(RESVEC(1:NGDIM), [NDIM,NDIM])
C Apply Newton correction to hardening variable
        HRVAR=HRVAR-RESVEC(NEQUA)
  360 CONTINUE
C End of Newton-Raphson iterations for visco-plastic system of equations
C ======================================================================
C N-R loop failed to converge: Switch state update failure flag on
C                              and exit without updating stresses
      SUFAIL=.TRUE.
      WRITE(*,*)'N-R loop failed! SUVCPL'!CALL ERRPRT('WE0006')
      GOTO 999
C N-R loop converged: Use neo-Hookean law to update stresses
C ==========================================================
  370 CONTINUE
C Compute inverse of elastic deformation gradient
C... elastic deformation gradient
      FE=FEISO/VOLFAC
C!!!! VOLFAC or VFACN ????
C!!!!!!!! replace FE by FEN (remove FE)
C... inverse of elastic deformation gradient
      CALL INVMT3(FE, FEINV, DETFE)
C Compute austenite plastic deformation gradient
      FPAUS=MATMUL(FTOTA, FEINV)
C Compute elastic left Cauchy-Green tensor
      BEISO=MATMUL(FEISO, TRANSPOSE(FEISO))
C Hydrostatic pressure
      P=BULK*LOG(DETFET)
C Deviatoric component of isochoric elastic left Cauchy-Green tensor
      TRACE=BEISO(1,1)+BEISO(2,2)+BEISO(3,3)
      CALL SYMTA(BEISO, BEDEV, NDIM, .TRUE.) ! Write BEISO in array form (symmetric)
      IF(IS2D)THEN
        BEDEV([1,2,4])=BEDEV([1,2,4])-R1D3*TRACE
      ELSE
        BEDEV([1,2,3])=BEDEV([1,2,3])-R1D3*TRACE
      ENDIF
C Update Cauchy stress components
      DETINV=R1/DETFET
      STRES=GMODU*BEDEV*DETINV
      IF(IS2D)THEN
        STRES([1,2,4])=STRES([1,2,4])+P*DETINV
      ELSE
        STRES([1,2,3])=STRES([1,2,3])+P*DETINV
      ENDIF
C Update state and algorithmic variables
C ======================================
      RSTAVA(1:9)=RESHAPE(FE,[9])
      RSTAVA(41)=HRVAR
      RSTAVA(42:50)=RESHAPE(FPAUS,[9])
  999 CONTINUE 
      LALGVA(2)=SUFAIL
      LALGVA(3)=IFPLAS
      RETURN
      END     
CDOC END_SUBROUTINE SUVCPL
