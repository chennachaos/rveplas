CDOC BEGIN_SUBROUTINE CSTVCP
CDOC Consistent spatial tangent modulus for planar visco-plastic single 
CDOC crystal 2D model.
CDOC
CDOC This routine computes the spatial tangent modulus consistent with
CDOC the implicit exponential map-based integration algorithm for the
CDOC planar visco-plastic single crystal 2D model
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION AMATX  <  Consistent spatial viscoplastic tangent
CDOC C                          modulus.
CDOC LOGICAL          CTFAIL <  Consistent tangent evaluation failure
CDOC C                          flag. Set to .TRUE. if tangent 
CDOC C                          evaluation fails. Set to .FAILS.
CDOC C                          otherwise.
CDOC DOUBLE_PRECISION DTIME  >  Time increment.
CDOC DOUBLE_PRECISION FINV   >  Inverse of current total deformation 
CDOC                            gradient.
CDOC DOUBLE_PRECISION FINCR  >  Current Incremental deformation
CDOC C                          gradient.
CDOC INTEGER          IPROPS >  Array of integer material properties.
CDOC C                          The number of points on the piece-wise
CDOC C                          linear hardening curve is the only
CDOC C                          element of this array used here.
CDOC C                          This array is set in routines
CDOC C                          RVDATA and RDVSC2.
CDOC INTEGER          NTYPE  >  Stress state type. The present
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
CDOC C                          RDVSC2.
CDOC DOUBLE_PRECISION RSTAVA >  Array of real state variables other
CDOC C                          than the stress tensor components.
CDOC C                          Current values (as output of integration
CDOC C                          algorithm).
CDOC C                          The state variables stored in
CDOC C                          this array are the in-plane components
CDOC C                          of the elastic deformation gradient and
CDOC C                          the Taylor hardening internal variable.
CDOC DOUBLE_PRECISION RSTAVN >  Array of real state variables other
CDOC C                          than the stress tensor components at
CDOC C                          the last (equilibrium) converged
CDOC C                          configuration.
CDOC DOUBLE_PRECISION STRES  <  Array of updated stress tensor
CDOC C                          components.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, April 1999: Initial coding
CHST G.Q.Liu,         April 1999: Merger to new implicit code
CHST                              in ELFEN
CHST F. Adziman,   November 2012: Implemented in HYPLAS v_3.0.0, 
CHST                              a rate-dependent version
CHST F. Adziman,    January 2013: Multi-slip systems adopted and
CHST                              Peric's model (with yield surfaces)  
CHST
      SUBROUTINE CSTVCP 
     1(   AMATX      ,CTFAIL     ,DTIME      ,FINV       ,FINCR       ,
     2    IPROPS     ,NTYPE      ,RPROPS     ,RSTAVA     ,RSTAVN      ,
     3    STRES      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER
     1(   IPHVAR=41  ,MDIM=3     ,NRSTAV=52    ,
     2    NSTRE=6    ,NIPROP=19  ,MGDIM=9)!MGDIM=5)
C Arguments
      LOGICAL
     1    CTFAIL     
      DIMENSION
C1
     1    AMATX(MGDIM,MGDIM)    ,FINV(MDIM,MDIM) ,FINCR(MDIM,MDIM)  ,
     2    IPROPS(NIPROP)        ,RPROPS(*)       ,RSTAVA(NRSTAV)    ,
     3    RSTAVN(NRSTAV)        ,STRES(NSTRE)
C Local arrays and variables
      LOGICAL
     1    IFGERR     ,NOCONV 
      DIMENSION
     1    ABSSCH(IPROPS(4))           ,OVERST(IPROPS(4))              ,
     2    SCHMID(IPROPS(4))           ,SIGNSC(IPROPS(4))              ,
     2    VECM(MDIM,IPROPS(4))        ,VECS(MDIM,IPROPS(4))           ,
     3    VECM0(MDIM,IPROPS(4))       ,VECS0(MDIM,IPROPS(4))          ,
     1    SM0MS0(MDIM,MDIM,IPROPS(4)) ,S0M0(MDIM,MDIM,IPROPS(4))      ,
C
     2    CMATX(MGDIM,MGDIM)          ,HMATX(MGDIM,MGDIM)             ,
     3    DFEDF(MGDIM,MGDIM)          ,AUXM25(MGDIM,MGDIM)            ,
     4    DF1DF(MGDIM,MGDIM)          ,DF1DFE(MGDIM,MGDIM)            ,
     5    DF1DA(MGDIM)                ,DF2DFE(MGDIM)                  ,
     6    AUXM1(MGDIM,MGDIM)          ,AUXM2(MGDIM,MGDIM)             ,
C
     1    A4TH(MDIM,MDIM,MDIM,MDIM)   ,
     2    C4TH(MDIM,MDIM,MDIM,MDIM)   ,H4TH(MDIM,MDIM,MDIM,MDIM)      ,
     6    DRE12(MDIM,MDIM)            ,DRE21(MDIM,MDIM)               ,
     7    DRE4TH(MDIM,MDIM,MDIM,MDIM) ,DSUMA(MDIM,MDIM)               ,
     8    DSUMFE(MDIM,MDIM,MDIM,MDIM) ,A2ND(MDIM,MDIM)                ,
     O    SIGMA(MDIM,MDIM)            ,DEREXP(MDIM,MDIM,MDIM,MDIM)    ,
C
     1    FE(MDIM,MDIM)               ,FEN(MDIM,MDIM)                 , 
     2    F(MDIM,MDIM)                ,FPILOG(MDIM,MDIM)              ,
     3    FPINV(MDIM,MDIM)            ,FPINCI(MDIM,MDIM)              , 
     4    FEINV(MDIM,MDIM)            ,FEISO(MDIM,MDIM)               ,
     5    FETISO(MDIM,MDIM)           ,FETRL(MDIM,MDIM)               ,
     6    FP(MDIM,MDIM)
      
      LOGICAL IS2D, ISOVER(IPROPS(4))
      
      integer seed(1)
      
C Local parameters
      DOUBLE PRECISION, PARAMETER :: R0=0.0D0        , R1=1.0D0        ,
     1                               R1D3=1.0D0/3.0D0, R2D3=2.0D0/3.0D0,
     2                               R4D3=4.0D0/3.0D0, SMALL=1.0D-12
      DOUBLE PRECISION, PARAMETER, DIMENSION(MDIM,MDIM) :: ! delta_ij
     1     DELTA=reshape((/R1,R0,R0,R0,R1,R0,R0,R0,R1/), (/3,3/) )
C**********************************************************************
C COMPUTE CONSISTENT SPATIAL TANGENT MODULUS FOR THE ANISOTROPIC PLANAR
C VISCO-PLASTIC SINGLE CRYSTAL 2D MODEL
C**********************************************************************
C Check analysis type is valid for this model
      IF(NTYPE==2)THEN
        NDIM=2
        NGDIM=4
        IS2D=.TRUE.
      ELSEIF(NTYPE==4)THEN
        NDIM=3
        NGDIM=9
        IS2D=.FALSE.
      ELSE
        STOP 'NTYPE invalid CSTVCP.f' !CALL ERRPRT('EI0007')
      ENDIF
C Initialise state update failure flag
      CTFAIL=.FALSE.
C Set current hardening internal variable
      HRVAR=RSTAVA(IPHVAR)
C... current elastic deformation gradient
      FE=RESHAPE(RSTAVA(1:9),[MDIM,MDIM])
C... recover current total deformation gradient
      CALL INVMT3(FINV, F, DETFIN)
C... recover current elastic deformation gradient
      CALL INVMT3(FE, FEINV, DETFE)
      FP=MATMUL(FEINV, F)
C!!!!!!! F and FP are incorrect!
C  -> F comes from inverting FINV (input argument to this subroutine),
C     that is evaluated incorrectly for FBAR elements (see f-bar routines)
C  -> FP should come from real state variables (FPAUS), or from
C     FTOTA=FE*FTR*FPAUS
C     => FPAUS=FTRINV*FEINV*FTOTA
C!!!!!!!!! not necessary to store all deformation gradients (storing
C          three allows the 4th to be derived directly
C... recover current plastic deformation gradient
      CALL INVMT3(FP, FPINV, DETFP)
C... previous converged elastic deformation gradient (equilibrium)
      FEN=RESHAPE(RSTAVN(1:9),[MDIM,MDIM])
C Retrieve material properties
C... neo-Hookean constants
      GMODU=RPROPS(2)
      BULK=RPROPS(3)
C... number of sampling points on hardening curve
      NHARD=IPROPS(3)
C Stops program if the selected slip-rate law is invalid             
C    1. MODEL=1 for the Norton's model (without yield surfaces)  
C    2. MODEL=2 for the Peric's model (with yield surfaces)
      MODEL=IPROPS(19)
      IF((MODEL/=1).AND.(MODEL/=2))THEN
        STOP 'ERROR: Slip-rate law chosen is invalid in CSTVCP'
      ENDIF
C... number of slip systems
      NSLSYS=IPROPS(4)
C Slip system vectors and matrix of constants
      VECM0=RESHAPE(RPROPS(IPROPS(12):IPROPS(13)),[MDIM,NSLSYS])
      VECS0=RESHAPE(RPROPS(IPROPS(14):IPROPS(15)),[MDIM,NSLSYS])
      S0M0=RESHAPE(RPROPS(IPROPS(16):IPROPS(17)),[MDIM,MDIM,NSLSYS])
C
      IPHARD=IPROPS(18)
C      
C Extract viscoplastic properties and set constants
      CONST1=RPROPS(8) !Norton:CONSTM       ; Peric: CONSEP
      R1DC=R1/CONST1
      R1DCM1=R1DC-R1
      CONST2=RPROPS(9) !Norton:G0DOT        ; Peric: CONSMU
      IF(MODEL==1)THEN
        DGAMMA=DTIME*CONST2
      ELSEIF(MODEL==2)THEN
        DGAMMA=DTIME/CONST2
      ENDIF
C Compute last elastic trial deformation gradient used by the
C integration algorithm
C -----------------------------------------------------------
C Last used elastic trial deformation gradient
      FETRL=MATMUL(FINCR, FEN)
C Isochoric component
      DETFET=DETM23(MDIM, FETRL, NDIM, .FALSE.)
      VOLFAC=DETFET**(-R1D3)
      FETISO=VOLFAC*FETRL
C Set up current isochoric elastic deformation gradient
C -----------------------------------------------------
      VOLFAC=DETFE**(-R1D3)
      FEISO=VOLFAC*FE
C Compute elastic push forward of all slip-systems vectors
C and current resolved Schmid stresses on all slip systems
      DO ISYST=1,NSLSYS
        VECS(:,ISYST)=MATMUL(FEISO, VECS0(:,ISYST))
        VECM(:,ISYST)=MATMUL(FEISO, VECM0(:,ISYST))
        SCHMID(ISYST)=GMODU*DOT_PRODUCT(VECS(:,ISYST),VECM(:,ISYST))
      ENDDO
C... absolute values
      ABSSCH=ABS(SCHMID)
C Current resolved reference stress / yield shear stress
      RSHEAR=PLFUN(HRVAR,NHARD,RPROPS(IPHARD)) !Norton: resolved reference stress
C                                              !Peric: yield resolved shear stress
C!!!!!! add checking if RSHEAR<R0 (like in ELFEM routines)
      OVERST=ABSSCH/RSHEAR
      SIGNSC=SCHMID/ABSSCH
C
C9!!!        
      ISOVER=OVERST > SMALL
C9I        IF(MODEL==1)THEN
C9I          ISOVER=OVERST > SMALL ! No yield surface!
C9I        ELSEIF(MODEL==2)THEN
C9I          ISOVER=OVERST-R1 > SMALL
C9I        ENDIF
C9
C
C Sum up inelastic slip contributions from all systems
C !!!(only systems that yield for Peric's law)
      FPILOG=R0
      DO ISYST=1,NSLSYS
        IF(ISOVER(ISYST))THEN
          IF(MODEL==1)THEN
            FACTOR=DGAMMA*SIGNSC(ISYST)*OVERST(ISYST)**R1DC
          ELSEIF(MODEL==2)THEN
            FACTOR=DGAMMA*SIGNSC(ISYST)*(OVERST(ISYST)**R1DC-R1)
C!!! Need to use signsc with Peric model???
          ENDIF
          FPILOG=FPILOG - FACTOR*S0M0(:,:,ISYST) ! FACTOR = gamma^dot * delta_T = delta_Gamma
        ENDIF !FPILOG = log of inverse (?because of - sign?) of incremental plastic deformation gradient
      ENDDO
C Use exponential map to compute inverse of incremental inelastic
C deformation gradient
      CALL EXPMAP
     1(   FPINCI     ,NOCONV     ,FPILOG     )
C!!!!!!!!!! necessary? -> FPINCI is never used again...
      IF(NOCONV)THEN
C Exponential map algorithm failed: Break loop and exit
        CTFAIL=.TRUE.
        WRITE(*,*)'exponential map failed - CSTVCP' !CALL ERRPRT('WE0030')
        GOTO 999
      ENDIF
C 
C Compute derivatives of residual of isochoric elastic def.grad. update
C =====================================================================
C Exponential map derivative
      CALL DEXPMP(   DEREXP      ,NOCONV     ,FPILOG     )
C!!!!!!!!!! deal with NOCONV here !!!!!!
C Derivative of summation of contributions of each slip plane
C4      DO I=1,NDIM
C4        DO J=1,NDIM
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
          FACTOR=AUX1*OVERST(ISYST)**R1DCM1
          FACTA=AUX2*SIGNSC(ISYST)*OVERST(ISYST)**R1DC
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
C4
      A4TH=R0
      DRE4TH=R0
C4
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
      DO I=1,NDIM
        DO J=1,NDIM
          DO K=1,NDIM
            DO L=1,NDIM
              DRE4TH(I,J,K,L)=DELTA(I,K)*DELTA(J,L)
     1                        +SUM( FETISO(I,:)*A4TH(:,J,K,L) )
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C Rearrange components as a 2-index matrix (DF1DFE) using G-matrix
C component ordering
      IF(IS2D)THEN
        DF1DFE=R0
        CALL ARRGO2
     1(   DRE4TH(1:NDIM,1:NDIM,1:NDIM,1:NDIM), DF1DFE(1:NGDIM,1:NGDIM))
      ELSE
        CALL ARRGO3
     1(   DRE4TH     ,DF1DFE     )
      ENDIF
C Derivative of residual of isochoric elastic def.grad. update with
C respect to hardening variable
C -----------------------------------------------------------------
C4
      A2ND=R0
C4
      DO I=1,NDIM
        DO J=1,NDIM
          A2ND(I,J)=SUM( DEREXP(I,J,:,:)*DSUMA(:,:) )
        ENDDO
      ENDDO
C
      DRE12=-MATMUL(FETISO, A2ND)
C... store in vector form (G-matrix component ordering)
C4
      DF1DA=R0
      DF2DFE=R0
C4
      DF1DA(1:NGDIM)=RESHAPE(DRE12(1:NDIM,1:NDIM),[NGDIM])
C Derivative of hardening residual with respect to isochoric elastic
C deformation gradient
C ------------------------------------------------------------------
      DRE21=R0
      AUX1=DGAMMA*R1DC*GMODU/RSHEAR
      DO ISYST=1,NSLSYS
        IF(ISOVER(ISYST))THEN
          FACTOR=-AUX1*SIGNSC(ISYST)*OVERST(ISYST)**R1DCM1
          DRE21=DRE21 + FACTOR*SM0MS0(:,:,ISYST)
        ENDIF
      ENDDO
C... store in vector form (G-matrix component ordering)
      DF2DFE(1:NGDIM)=RESHAPE(DRE21(1:NDIM,1:NDIM),[NGDIM])
C Derivative of hardening residual with respect to hardening variable
C -------------------------------------------------------------------
      FACTOR=DGAMMA*R1DC*HSLOPE/RSHEAR
      DF2DA=R1 + FACTOR*SUM(OVERST**R1DC, MASK=ISOVER)
C Derivative of residual of isochoric elastic def.grad. update with
C respect to TOTAL isochoric def.grad.
C -----------------------------------------------------------------
C4
      A4TH=R0
C4
      DO I=1,NDIM
        DO J=1,NDIM
          DO K=1,NDIM
            DO L=1,NDIM
              A4TH(I,J,K,L)=-DELTA(I,K)*FPINV(L,J)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C... store in matrix form (G-matrix component ordering)
      IF(IS2D)THEN
        DF1DF=R0
        CALL ARRGO2
     1(   A4TH(1:NDIM,1:NDIM,1:NDIM,1:NDIM)       ,DF1DF(1:4,1:4) )
      ELSE
        CALL ARRGO3
     1(   A4TH       ,DF1DF      )
      ENDIF
C Assemble dFe/dF (algorithm-consistent deriv. of elastic isochoric
C def.grad. with respect to total isochoric def.grad.)
C =================================================================
C4
C4      AUXM1=R0
C4
      DO I=1,NGDIM
        DO J=1,NGDIM
          AUXM1(I,J)=-(DF1DFE(I,J)-DF1DA(I)*DF2DFE(J)/DF2DA)
        END DO
      END DO
C1
      AUXM2=R0
C1
      CALL RMINVE
     1( AUXM1(1:NGDIM,1:NGDIM) ,AUXM2(1:NGDIM,1:NGDIM) ,NGDIM ,IFGERR )
      IF(IFGERR)THEN
        CTFAIL=.TRUE.
        WRITE(*,*)'RMINVE failed! CSTVCP'!CALL ERRPRT('WE0031')
        GOTO 999
      ENDIF
      DFEDF=R0
C1      DFEDF(1:NGDIM,1:NGDIM)=MATMUL(AUXM2(1:NGDIM,1:NGDIM), 
C1     1                              DF1DF(1:NGDIM,1:NGDIM))
      DFEDF=MATMUL(AUXM2, DF1DF)
      DFEDF(5,5)=R1
      
C!!!!!!!!!!!! USE GAUSEL instead of matrix inversion possible?
      
C Compute C matrix
C ================
C First as 4th order tensor
C!!!! NDIM or MDIM??
      DO I=1,MDIM
        DO J=1,MDIM
          DO K=1,MDIM
            DO L=1,MDIM
              C4TH(I,J,K,L)=DELTA(I,K)*FEISO(J,L)+FEISO(I,L)*DELTA(J,K)
     1                      -R2D3*DELTA(I,J)*FEISO(K,L)
C!!!!!!! 1/3 or 2/3 ??? -> see EADSN p 726
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C!!!!!!! storing things in G-matrix ordering allows double contractions
C        between two 4th-order tensors to be expressed as matrix products
C        -> explore in SUMEPC
C...then store in matrix form (G-matrix component ordering)
      IF(IS2D)THEN
        CMATX=R0
C1        CALL ARRGO2
C1     1(   C4TH(1:NDIM,1:NDIM,1:NDIM,1:NDIM)  ,CMATX(1:NGDIM,1:NGDIM) )
        CALL ARRGAX
     1(   C4TH  ,CMATX(1:5,1:5) )
      ELSE
        CALL ARRGO3
     1(   C4TH      ,CMATX       )
      ENDIF
C Compute H matrix
C ================
C First as 4th order tensor      
      DO I=1,MDIM
        DO J=1,MDIM
          DO K=1,MDIM
            DO L=1,MDIM
              H4TH(I,J,K,L)=DELTA(I,K)*F(L,J)-R1D3*F(I,J)*DELTA(K,L)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C...then store in matrix form (G-matrix component ordering)
      IF(IS2D)THEN
        HMATX=R0
C1        CALL ARRGO2
C1     1(   H4TH(1:NDIM,1:NDIM,1:NDIM,1:NDIM)  ,HMATX(1:NGDIM,1:NGDIM) )
        CALL ARRGAX
     1(   H4TH  ,HMATX(1:5,1:5) )
      ELSE
        CALL ARRGO3
     1(   H4TH      ,HMATX       )
      ENDIF
      
C1      WRITE(*,*)'DFEDF'
C1      DO I=1,9
C1        WRITE(*,'(9G20.6)')DFEDF(I,:)
C1      ENDDO
C1      WRITE(*,*)'HMATX'
C1      DO I=1,9
C1        WRITE(*,'(9G20.6)')HMATX(I,:)
C1      ENDDO
C1      WRITE(*,*)'CMATX'
C1      DO I=1,9
C1        WRITE(*,'(9G20.6)')CMATX(I,:)
C1      ENDDO
C1      READ(*,*)
      
      
      
C Assemble spatial tangent "a"
C ============================
      CALL ATSYM( STRES, SIGMA, NDIM, .FALSE.)
C!!!!! .TRUE. or .FALSE. ???????
      FACTOR=BULK/DETFE
C4
      A4TH=R0
C4
      DO I=1,NDIM
        DO J=1,NDIM
          DO K=1,NDIM
            DO L=1,NDIM
              A4TH(I,J,K,L)=FACTOR*DELTA(I,J)*DELTA(K,L)
     1                      -SIGMA(I,L)*DELTA(J,K)
            END DO
          END DO
        END DO
      END DO
      
      IF(IS2D)THEN
        AUXM1=R0
        CALL ARRGO2
     1(   A4TH(1:NDIM,1:NDIM,1:NDIM,1:NDIM) ,AUXM1(1:NGDIM,1:NGDIM)    )
      ELSE
        CALL ARRGO3
     1(   A4TH      ,AUXM1       )
      ENDIF
C4
      AUXM25=R0
C4      AUXM25(1:NGDIM,1:NGDIM)=MATMUL(MATMUL(CMATX(1:NGDIM,1:NGDIM), 
C4     1                                      DFEDF(1:NGDIM,1:NGDIM)), 
C4     2                               HMATX(1:NGDIM,1:NGDIM))
C5      AUXM25(1:5,1:5)=MATMUL(MATMUL(CMATX(1:5,1:5), 
C5     1                              DFEDF(1:5,1:5)) , 
C5     2                       HMATX(1:5,1:5))
      AUXM25=MATMUL(MATMUL(CMATX, DFEDF), HMATX)
C... finally get "a" matrix
      FACTOR=GMODU/(DETFE**(R4D3))
      AMATX=FACTOR*AUXM25+AUXM1
C... if consistent tangent failed
  999 CONTINUE
      RETURN
      END
CDOC END_SUBROUTINE CSTVCP
