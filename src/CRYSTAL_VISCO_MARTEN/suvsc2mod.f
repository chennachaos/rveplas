CDOC BEGIN_SUBROUTINE SUVPCR
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
C8    SUBROUTINE SUVPCR
      SUBROUTINE SUVSC2MOD
     1(   DTIME      ,FINCR      ,IPROPS     ,LALGVA     ,NTYPE      ,
     2    RPROPS     ,RSTAVA     ,STRES      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER
     1(   MDIM=3     ,MGDIM=9    ,MEQUA=10    ,NIPROP=19   ,
CC1     2    NSTRE=6    ,NLALGV=2   ,NRSTAV=19   ) !IPSYST=10,IPHVAR=41 
     2    NSTRE=6    ,NLALGV=2   ,NRSTAV=52   ) !IPSYST=10,IPHVAR=41 
      LOGICAL
     1    LALGVA(NLALGV)
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
      DIMENSION DG(IPROPS(4))
      DOUBLE PRECISION DGDFEI(3,3,IPROPS(4)), DGDA(IPROPS(4)),
     1      FTIDSM(3,3,IPROPS(4))
      DOUBLE PRECISION DF1DF4(3,3,3,3), DF1DFE(MGDIM,MGDIM), DF1DA(3,3), 
     1                 DF2DFE(3,3)
      
      DOUBLE PRECISION, DIMENSION(MDIM, MDIM) ::
     1    FEINV ,FEISO ,FETISO ,FETFPI ,
     9    FETRL ,FPILOG , FPINCI ,FEN, FPN
C
C
      LOGICAL ISOVER(IPROPS(4)), IS2D
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
        STOP 'NTYPE invalid SUVSC2.f' !CALL ERRPRT('EI0006')
      ENDIF
C Initialise state update failure flag
      IFPLAS=.FALSE.
      SUFAIL=.FALSE.
C Set previous converged hardening internal variable
      HRVARN=RSTAVA(19)
      HRVAR=HRVARN
C Retrieve elastic deformation gradient
      FEN=RESHAPE(RSTAVA(1:9),[MDIM,MDIM])
C Retrieve material properties:
C - neo-Hookean constants
      GMODU=RPROPS(2)
      BULK=RPROPS(3)
C ... number of sampling points on hardening curve
      NHARD=IPROPS(3)
C - slip systems information
      NSLSYS=IPROPS(4)
      VECM0=RESHAPE(RPROPS(IPROPS(12):IPROPS(13)),[MDIM,NSLSYS])
      VECS0=RESHAPE(RPROPS(IPROPS(14):IPROPS(15)),[MDIM,NSLSYS])
      S0M0=RESHAPE(RPROPS(IPROPS(16):IPROPS(17)),[MDIM,MDIM,NSLSYS])
C
      IPHARD=IPROPS(18)
      MODEL=IPROPS(19)
C??? add checking if MODEL is valid?
C Compute deformation gradient components
C ---------------------------------------
C Elastic trial deformation gradient
      FETRL=MATMUL(FINCR, FEN)
C Perform isochoric/volumetric split of elastic trial def. grad.
      DETFET=DETM23(MDIM, FETRL, NDIM, .FALSE.)
      VOLFAC=DETFET**(-R1D3)
      FETISO=VOLFAC*FETRL
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
C Check plastic consistency
C -------------------------
C Isochoric elastic push forward of system vectors and resolved Schmid
C stresses for all slip systems
      DO ISYST=1,NSLSYS
        VECS(:,ISYST)=MATMUL(FETISO, VECS0(:,ISYST))
        VECM(:,ISYST)=MATMUL(FETISO, VECM0(:,ISYST))
        SCHMID(ISYST)=GMODU*DOT_PRODUCT(VECS(:,ISYST),VECM(:,ISYST))
      ENDDO
C Reference shear stress (for Norton creep model) / Yield stress 
C (for Peric and Perzyna models)
      RSHEAR=PLFUN(HRVAR,NHARD,RPROPS(IPHARD))
      OVERST=ABS(SCHMID)/RSHEAR
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
C      IF(ANY(ISOVER)) IFPLAS=.TRUE.
      IF(ANY(ISOVER))THEN
        IFPLAS=.TRUE.
C
C      IF(IFPLAS)THEN !PLASTIC STEP
C 
C Start Newton-Raphson iterations for visco-plastic system of equations
C =====================================================================
      DO 360 NRITER=1,MXITER
C Compute residual vector
C -----------------------
C Isochoric elastic push forward of system vectors and resolved Schmid
C stresses for all slip systems
        DO ISYST=1,NSLSYS
          VECS(:,ISYST)=MATMUL(FEISO, VECS0(:,ISYST))
          VECM(:,ISYST)=MATMUL(FEISO, VECM0(:,ISYST))
          SCHMID(ISYST)=GMODU*DOT_PRODUCT(VECS(:,ISYST),VECM(:,ISYST))
        ENDDO
        SIGNSC=SIGN(R1, SCHMID) ! Sign of Schmid stresses
C Sum up inelastic slip contributions from all systems
C !!!(only systems that yield for Peric's law)
        RSHEAR=PLFUN(HRVAR,NHARD,RPROPS(IPHARD))
        OVERST=ABS(SCHMID)/RSHEAR
        IF(MODEL==1)THEN ! Norton: no yield surface
          ISOVER=OVERST > SMALL
        ELSE ! Peric, Perzyna
          ISOVER=OVERST-R1 > SMALL
        ENDIF
        HSLOPE=DPLFUN(HRVAR,NHARD,RPROPS(IPHARD))
C!!!!!
        IF(MODEL==1)THEN
          DTIMEC=DTIME*CONST2
        ELSE
          DTIMEC=DTIME/CONST2
        ENDIF
        FACTOR=DTIMEC*R1DC*HSLOPE/RSHEAR
        IF(MODEL==1)THEN
          DG=DTIMEC*SIGNSC*OVERST**R1DC
          DGDA=-FACTOR*SIGNSC*OVERST**R1DC
        ELSEIF(MODEL==2)THEN
          DG=DTIMEC*SIGNSC*(OVERST**R1DC-R1)
          DGDA=-FACTOR*SIGNSC*OVERST**R1DC
        ELSE ! Perzyna
          DG=DTIMEC*SIGNSC*(OVERST-R1)**R1DC
          DGDA=-FACTOR*SIGNSC*OVERST*(OVERST-R1)**R1DCM1
        ENDIF
        FACTOR=DTIMEC*R1DC/RSHEAR*GMODU
        DO I=1,NDIM
          DO J=1,NDIM
            DGDFEI(I,J,:)=FACTOR*OVERST**R1DCM1*( VECS(I,:)*VECM0(J,:)
     1                                           +VECM(I,:)*VECS0(J,:) )
          ENDDO
        ENDDO
        
        
CC1
C Current austenite volume fraction
CC1        AUSVF=R1 - SUM(RSTAVA(20:43))
        
        DO I=1,NDIM
          DO J=1,NDIM
            FPILOG(I,J)=-SUM(DG(:)*S0M0(I,J,:), MASK=ISOVER)
CC1            FPILOG(I,J)=-SUM(AUSVF*DG(:)*S0M0(I,J,:), MASK=ISOVER)
          ENDDO
        ENDDO
        
        CALL DEXPMP(   DEREXP      ,NOCONV     ,FPILOG     )
        IF(NOCONV)THEN
C Derivative of exponential map failed: Break loop and exit
          SUFAIL=.TRUE.
CC1          WRITE(*,*)'Derivative exp map failed - SUVCPL'
          GOTO 999
        ENDIF
        
C Fe trial iso * (double contraction of DEREXP and S0M0)
C SKip calculations for systems where .NOT.ISOVER
        DO ISYST=1,NSLSYS
          IF(ISOVER(ISYST))THEN
            DO I=1,NDIM
              DO J=1,NDIM
                FTIDSM(I,J,ISYST)=SUM( DEREXP(I,J,:,:)*S0M0(:,:,ISYST) )
              ENDDO
            ENDDO
            FTIDSM(:,:,ISYST)=MATMUL(FETISO, FTIDSM(:,:,ISYST))
          ENDIF
        ENDDO
        
        DO I=1,NDIM
          DO J=1,NDIM
            DO K=1,NDIM
              DO L=1,NDIM
                DF1DF4(I,J,K,L)=DELTA(I,K)*DELTA(J,L)
     1                    +SUM(FTIDSM(I,J,:)*DGDFEI(K,L,:), MASK=ISOVER)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        IF(IS2D)THEN
          CALL ARRGO2
     1(   DF1DF4(1:NDIM,1:NDIM,1:NDIM,1:NDIM) ,DF1DFE(1:NGDIM,1:NGDIM) )
        ELSE
          CALL ARRGO3
     1(   DF1DF4 ,DF1DFE )
        ENDIF
C        
        DO I=1,NDIM
          DO J=1,NDIM
            DF1DA(I,J)=SUM(FTIDSM(I,J,:)*DGDA(:), MASK=ISOVER)
          ENDDO
        ENDDO
        DO I=1,NDIM
          DO J=1,NDIM
            DF2DFE(I,J)=-SUM(SIGNSC*DGDFEI(I,J,:), MASK=ISOVER)
          ENDDO
        ENDDO
        
        DF2DA=R1 - SUM(SIGNSC*DGDA, MASK=ISOVER)
C Assemble complete Jacobian matrix
C ---------------------------------
        DREMTX(1:NGDIM,1:NGDIM)=DF1DFE(1:NGDIM,1:NGDIM)
        DREMTX(1:NGDIM,NEQUA)=RESHAPE(DF1DA(1:NDIM,1:NDIM), [NGDIM])
        DREMTX(NEQUA,1:NGDIM)=RESHAPE(DF2DFE(1:NDIM,1:NDIM), [NGDIM])
        DREMTX(NEQUA,NEQUA)=DF2DA
C Use exponential map to compute inverse of incremental inelastic
C deformation gradient
        CALL EXPMAP
     1(   FPINCI     ,NOCONV     ,FPILOG     )
        IF(NOCONV)THEN
C Exponential map algorithm failed: Break loop and exit
          SUFAIL=.TRUE.
CC1          WRITE(*,*)'exponential map failed - SUVCPL' !CALL ERRPRT('WE0005')
          GOTO 999
        ENDIF
C Compute residual arranged in vector form (using G matrix component
C ordering)
        FETFPI=MATMUL(FETISO, FPINCI)
C Calculate residual vector
        RESVEC(1:NGDIM)=RESHAPE(
     1            FEISO(1:NDIM,1:NDIM) - FETFPI(1:NDIM,1:NDIM), [NGDIM])
        RESVEC(NEQUA)=HRVAR-HRVARN-SUM(ABS(DG), MASK=ISOVER)
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
C!!!!!!!
C9        FPILOG=R0
C9        DO ISYST=1,NSLSYS
C9          IF(ISOVER(ISYST))THEN
C9            IF(MODEL==1)THEN
C9              FACTOR=DGAMMA*SIGNSC(ISYST)*OVERST(ISYST)**R1DC
C9            ELSEIF(MODEL==2)THEN
C9              FACTOR=DGAMMA*SIGNSC(ISYST)*(OVERST(ISYST)**R1DC-R1)
C9            ELSE ! Perzyna
C9              FACTOR=DGAMMA*SIGNSC(ISYST)*(OVERST(ISYST)-R1)**R1DC
C9            ENDIF
C9            FPILOG=FPILOG - FACTOR*S0M0(:,:,ISYST) ! FACTOR = gamma^dot * delta_t = delta_Gamma
C9          ENDIF !FPILOG = log of inverse (?because of - sign?) of incremental plastic deformation gradient
C9        ENDDO
C9C Use exponential map to compute inverse of incremental inelastic
C9C deformation gradient
C9        CALL EXPMAP
C9     1(   FPINCI     ,NOCONV     ,FPILOG     )
C9        IF(NOCONV)THEN
C9C Exponential map algorithm failed: Break loop and exit
C9          SUFAIL=.TRUE.
C9          WRITE(*,*)'exponential map failed - SUVCPL' !CALL ERRPRT('WE0005')
C9          GOTO 999
C9        ENDIF
C9C Compute residual arranged in vector form (using G matrix component
C9C ordering)
C9        FETFPI=MATMUL(FETISO, FPINCI)
C9C Calculate residual vector
C9        RESVEC(1:NGDIM)=RESHAPE(
C9     1            FEISO(1:NDIM,1:NDIM) - FETFPI(1:NDIM,1:NDIM), [NGDIM])
C9C      
C9        IF(MODEL==1)THEN
C9          RSUM=SUM(OVERST**R1DC, MASK=ISOVER)
C9        ELSEIF(MODEL==2)THEN
C9          RSUM=SUM(OVERST**R1DC-R1, MASK=ISOVER)
C9        ELSEIF(MODEL==3)THEN
C9          RSUM=SUM((OVERST-R1)**R1DC, MASK=ISOVER)
C9        ENDIF
C9        RESVEC(NEQUA)=HRVAR-HRVARN-DGAMMA*RSUM
C9C Check for convergence
C9C ---------------------
C9        RESNOR=NORM2(RESVEC(1:NEQUA))
C9        FACTOR=SQRT( NORM2(FEISO(1:NDIM,1:NDIM))**2+ABS(HRVAR-HRVARN) )
C9        IF(FACTOR/=R0)RESNOR=RESNOR/FACTOR
C9        IF(RESNOR<TOL)THEN
C9C Newton-Raphson iterations converged: Break loop and update stresses
C9C before exit
C9          GOTO 370
C9        ENDIF
C9C Compute derivatives of residual of isochoric elastic def.grad. update
C9C ---------------------------------------------------------------------
C9C Exponential map derivative
C9        CALL DEXPMP(   DEREXP      ,NOCONV     ,FPILOG     )
C9C!!!!!!
C9C!!!!!! CAN PASS SUFAIL directly to DEXPMP and other routines
C9        IF(NOCONV)THEN
C9C Derivative of exponential map failed: Break loop and exit
C9          SUFAIL=.TRUE.
C9          WRITE(*,*)'Derivative exp map failed - SUVCPL'
C9          GOTO 999
C9        ENDIF
C9C!!!!!!
C9C Derivative of summation of contributions of each slip plane
C9        DO I=1,MDIM
C9          DO J=1,MDIM
C9            SM0MS0(I,J,:)=VECS(I,:)*VECM0(J,:)+VECM(I,:)*VECS0(J,:)
C9          ENDDO
C9        ENDDO
C9C
C9        HSLOPE=DPLFUN(HRVAR,NHARD,RPROPS(IPHARD))
C9        DSUMFE=R0
C9        DSUMA=R0
C9        AUX1=DGAMMA*GMODU*R1DC/RSHEAR
C9        AUX2=DGAMMA*HSLOPE*R1DC/RSHEAR
C9C
C9        DO ISYST=1,NSLSYS
C9          IF(ISOVER(ISYST))THEN
C9            IF((MODEL==1).OR.(MODEL==2))THEN
C9              FACTOR=AUX1*OVERST(ISYST)**R1DCM1
C9              FACTA=AUX2*SIGNSC(ISYST)*OVERST(ISYST)**R1DC
C9            ELSE
C9              FACTOR=AUX1*(OVERST(ISYST)-R1)**R1DCM1
C9              FACTA=AUX2*SIGNSC(ISYST)
C9     1                  *OVERST(ISYST)*(OVERST(ISYST)-R1)**R1DCM1
C9            ENDIF
C9            DSUMA=DSUMA + FACTA*S0M0(:,:,ISYST)
C9C
C9            DO I=1,NDIM
C9              DO J=1,NDIM
C9                DO K=1,NDIM
C9                  DO L=1,NDIM
C9                    DSUMFE(I,J,K,L)=DSUMFE(I,J,K,L)
C9     1                         +FACTOR*S0M0(I,J,ISYST)*SM0MS0(K,L,ISYST)
C9                  ENDDO
C9                ENDDO
C9              ENDDO
C9            ENDDO
C9          ENDIF
C9        ENDDO
C9C
C9        A4TH=R0
C9        DO I=1,NDIM
C9          DO J=1,NDIM
C9            DO K=1,NDIM
C9              DO L=1,NDIM
C9                A4TH(I,J,K,L)=SUM( DEREXP(I,J,:,:)*DSUMFE(:,:,K,L) )
C9              ENDDO
C9            ENDDO
C9          ENDDO
C9        ENDDO
C9C Assemble (4th order tensor) derivative with respect to isochoric
C9C elastic def.grad. as a 4-index matrix
C9        DRE4TH=R0
C9        DO I=1,NDIM
C9          DO J=1,NDIM
C9            DO K=1,NDIM
C9              DO L=1,NDIM
C9                DRE4TH(I,J,K,L)=DELTA(I,K)*DELTA(J,L)
C9     1                          +SUM( FETISO(I,:)*A4TH(:,J,K,L) )
C9              ENDDO
C9            ENDDO
C9          ENDDO
C9        ENDDO
C9C Rearrange components as a 2-index matrix (DRE11) using G-matrix
C9C component ordering
C9        IF(IS2D)THEN
C9          CALL ARRGO2
C9     1(   DRE4TH(1:NDIM,1:NDIM,1:NDIM,1:NDIM) ,DRE11(1:NGDIM,1:NGDIM) )
C9        ELSE
C9          CALL ARRGO3
C9     1(   DRE4TH ,DRE11 )
C9        ENDIF
C9C Derivative of residual of isochoric elastic def.grad. update with
C9C respect to hardening variable
C9C -----------------------------------------------------------------
C9        A2ND=R0
C9        DO I=1,NDIM
C9          DO J=1,NDIM
C9            A2ND(I,J)=SUM( DEREXP(I,J,:,:)*DSUMA(:,:) )
C9          ENDDO
C9        ENDDO
C9C
C9        DRE12=-MATMUL(FETISO, A2ND)
C9C Derivative of hardening residual with respect to isochoric elastic
C9C def.grad.
C9C ------------------------------------------------------------------
C9        DRE21=R0
C9        AUX1=DGAMMA*R1DC*GMODU/RSHEAR
C9        DO ISYST=1,NSLSYS
C9          IF(ISOVER(ISYST))THEN
C9            IF((MODEL==1).OR.(MODEL==2))THEN
C9              FACTOR=-AUX1*SIGNSC(ISYST)*OVERST(ISYST)**R1DCM1
C9            ELSE
C9              FACTOR=-AUX1*SIGNSC(ISYST)*(OVERST(ISYST)-R1)**R1DCM1
C9            ENDIF
C9            DRE21=DRE21 + FACTOR*SM0MS0(:,:,ISYST)
C9          ENDIF
C9        ENDDO
C9C Derivative of hardening residual with respect to hardening variable
C9C -------------------------------------------------------------------
C9        FACTOR=DGAMMA*R1DC*HSLOPE/RSHEAR
C9        IF((MODEL==1).OR.(MODEL==2))THEN
C9          DF2DA=R1 + FACTOR*SUM(OVERST**R1DC, MASK=ISOVER)
C9        ELSE
C9          DF2DA=R1 + FACTOR*SUM(OVERST*(OVERST-R1)**R1DCM1, MASK=ISOVER)
C9        ENDIF
C9C Assemble complete Jacobian matrix
C9C ---------------------------------
C9        DREMTX(1:NGDIM,1:NGDIM)=DRE11(1:NGDIM,1:NGDIM)
C9        DREMTX(1:NGDIM,NEQUA)=RESHAPE(DRE12(1:NDIM,1:NDIM), [NGDIM])
C9        DREMTX(NEQUA,1:NGDIM)=RESHAPE(DRE21(1:NDIM,1:NDIM), [NGDIM])
C9        DREMTX(NEQUA,NEQUA)=DF2DA

        
        
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
CC1          WRITE(*,*)'GAUSEL failed! SUVCPL'!CALL ERRPRT('WE0007')
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
CC1      WRITE(*,*)'N-R loop failed! SUVCPL'!CALL ERRPRT('WE0006')
      GOTO 999
      
      
C4
      ELSE !ELASTIC STEP
C Trial state is correct
        FEISO=FETISO
      ENDIF
C4
      
      
C N-R loop converged: Use neo-Hookean law to update stresses
C ==========================================================
  370 CONTINUE
C Update elastic deformation gradient
      FEN=FEISO/VOLFAC
C9!!!! VOLFAC or VFACN ????
      CALL INVMT3(FEN, FEINV, DETFE)
C Retrieve and update plastic deformation gradient
      FPN=RESHAPE(RSTAVA(10:18),[MDIM,MDIM])
      FPN=MATMUL(MATMUL(FEINV, FETRL), FPN)
C Update Cauchy stresses
      BEISO=MATMUL(FEISO, TRANSPOSE(FEISO))
      CALL SYMTA(BEISO, BEDEV, NDIM, .TRUE.) ! Write BEISO in array form (symmetric)
C Deviatoric component of isochoric elastic left Cauchy-Green tensor
      TRACE=BEISO(1,1)+BEISO(2,2)+BEISO(3,3)
      IF(IS2D)THEN
        BEDEV([1,2,4])=BEDEV([1,2,4])-R1D3*TRACE
      ELSE
        BEDEV([1,2,3])=BEDEV([1,2,3])-R1D3*TRACE
      ENDIF
C Update Cauchy stress components
      STRES=GMODU*BEDEV/DETFET
      P=BULK*LOG(DETFET)/DETFET !Hydrostatic pressure
      IF(IS2D)THEN
        STRES([1,2,4])=STRES([1,2,4])+P
      ELSE
        STRES([1,2,3])=STRES([1,2,3])+P
      ENDIF
C Update state and algorithmic variables
C ======================================
      RSTAVA(1:9)=RESHAPE(FEN,[9])
      RSTAVA(10:18)=RESHAPE(FPN,[9])
      RSTAVA(19)=HRVAR
C!!! Add euler angles output
      
  999 CONTINUE 
      LALGVA(1)=IFPLAS
      LALGVA(2)=SUFAIL
      RETURN
      END     
CDOC END_SUBROUTINE SUVPCR
