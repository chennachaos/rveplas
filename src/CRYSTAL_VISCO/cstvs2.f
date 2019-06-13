CDOC BEGIN_SUBROUTINE CSTVPC
CDOC Consistent spatial tangent modulus for viscoplastic single crystal
CDOC model. Plane strain and 3D implementations.
CDOC
CDOC This routine computes the spatial tangent modulus consistent with
CDOC the implicit exponential map-based integration algorithm for the
CDOC viscoplastic single crystal model.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION AMATX  <  Consistent spatial viscoplastic tangent
CDOC C                          modulus.
CDOC DOUBLE_PRECISION DTIME  >  Time increment.
CDOC LOGICAL          EPFLAG >  Elasto-plastic flag.
CDOC C                          If .FALSE., AMATX returns as the elastic
CDOC C                          modulus. If .TRUE., AMATX returns as the
CDOC C                          elasto-plastic modulus.
CDOC DOUBLE_PRECISION FINCR  >  Current Incremental deformation
CDOC C                          gradient.
CDOC INTEGER          IPROPS >  Array of integer material properties.
CDOC C                          The number of points on the piece-wise
CDOC C                          linear hardening curve is the only
CDOC C                          element of this array used here.
CDOC C                          This array is set in routine RDVSC2.
CDOC INTEGER          NTYPE  >  Stress state type. The present
CDOC C                          implementation is compatible only with
CDOC C                          plane strain (NTYPE=2).
CDOC DOUBLE_PRECISION RPROPS >  Array of real material properties.
CDOC C                          This array contains the density
CDOC C                          (not used in this routine), the elastic
CDOC C                          properties: shear and bulk modulus
CDOC C                          (neo-Hookean constants) and the plastic
CDOC C                          properties: the hardening curve pairs "Taylor hardening
CDOC C                          variable-resolved reference stress"
CDOC C                          that define the (user supplied)
CDOC C                          piece-wise linear Taylor hardening
CDOC C                          curve. This array is set in routine
CDOC C                          RDVPCR.
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
CHST D. de Bortoli,     May 2016: 3-D implementation;
CHST                              added Perzyna slip rate law;
CHST                              removed unused arguments;
CHST                              added argument EPFLAG and elastic
CHST                              tangent calculation
CHST
C2    SUBROUTINE CSTVPC
      SUBROUTINE CSTVS2
     1(   AMATX ,  DTIME ,  EPFLAG,  FINCR,  IPROPS,  NTYPE,  RPROPS,
     2    RSTAVA,  RSTAVN,  STRES   )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER
     1(   MDIM=3, MGDIM=9, NSTRE=6, NRSTAV=19, NIPROP=19)
C Arguments
      INTEGER, PARAMETER :: IPNSYS=4
C!!!! Keep track of NRPROP?
      DIMENSION
C1
     1    AMATX(MGDIM,MGDIM)    ,FINCR(MDIM,MDIM),
     2    IPROPS(NIPROP)        ,RPROPS(*)       ,RSTAVA(NRSTAV)    ,
     3    RSTAVN(NRSTAV)        ,STRES(NSTRE)
C Local arrays and variables
      LOGICAL EPFLAG
      LOGICAL
     1    IFGERR     ,NOCONV 
      DIMENSION ! IPROPS(4) = number of slip systems
     1    ABSSCH(IPROPS(4))           ,OVERST(IPROPS(4))              ,
     2    SCHMID(IPROPS(4))           ,SIGNSC(IPROPS(4))              ,
     2    VECM(MDIM,IPROPS(4))        ,VECS(MDIM,IPROPS(4))           ,
     3    VECM0(MDIM,IPROPS(4))       ,VECS0(MDIM,IPROPS(4))          ,
     1    SM0MS0(MDIM,MDIM,IPROPS(4)) ,S0M0(MDIM,MDIM,IPROPS(4))      ,
C
     2    CMATX(MGDIM,MGDIM)          ,HMATX(MGDIM,MGDIM)             ,
     3    DFEDF(MGDIM,MGDIM)          ,AUXM25(MGDIM,MGDIM)            ,
     4    DF1DF(MGDIM,MGDIM)          ,!DF1DFE(MGDIM,MGDIM)            ,
C     5    DF1DA(MGDIM)                ,DF2DFE(MGDIM)                  ,
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
      
      DOUBLE PRECISION DEVSTR(MGDIM)
      
      
      DIMENSION DG(IPROPS(4))
      DOUBLE PRECISION DGDFEI(3,3,IPROPS(4)), DGDA(IPROPS(4)),
     1      FTIDSM(3,3,IPROPS(4))
      DOUBLE PRECISION DF1DF4(3,3,3,3), DF1DA(3,3), 
     1                 DF2DFE(3,3)
      DOUBLE PRECISION AUXM14(3,3,3,3)
      LOGICAL IS2D, ISOVER(IPROPS(4))
C Local parameters
      DOUBLE PRECISION, PARAMETER :: R0=0.0D0        , RP5=0.5D0       ,
     1                               R1=1.0D0        , R2=2.0D0        ,
     2                               R1D3=1.0D0/3.0D0, R2D3=2.0D0/3.0D0,
     3                               R4D3=4.0D0/3.0D0, SMALL=1.0D-12
      DOUBLE PRECISION, PARAMETER, DIMENSION(MDIM,MDIM) :: ! delta_ij
     1     DELTA=reshape((/R1,R0,R0,R0,R1,R0,R0,R0,R1/), (/3,3/) )
      
      
C 3D G matrix ordering (11,21,31,12,22,32,13,23,33)
C                        1  2  3  4  5  6  7  8  9
C 21,12 -> 2,4
C 31,13 -> 3,7
C 32,23 -> 6,8
C
C delta_ij - G-matrix ordering representation
      DOUBLE PRECISION, PARAMETER, DIMENSION(4) :: 
     1                                SOID2=[R1,R0,R0,R1]
      DOUBLE PRECISION, PARAMETER, DIMENSION(9) :: 
     1                                SOID3=[R1,R0,R0,R0,R1,R0,R0,R0,R1]
C Fourth-order identity tensor (symmetric subspace) - matrix
C representation using G matrix ordering
C 1/2*(delta_ik*delta_jl + delta_il*delta_jk)
      DOUBLE PRECISION, PARAMETER, DIMENSION(4,4) :: 
     1    FOIDS2=reshape((/R1 ,R0 ,R0 ,R0 , ! 1
     2                          R0 ,RP5,RP5,R0 , ! 2
     3                          R0 ,RP5,RP5,R0 , ! 3
     4                          R0 ,R0 ,R0 ,R1/), (/4,4/)) ! 4
C                               1   2   3   4
      DOUBLE PRECISION, PARAMETER, DIMENSION(9,9) :: 
     1    FOIDS3=reshape((/R1 ,R0 ,R0 ,R0 ,R0 ,R0 ,R0 ,R0 ,R0 , ! 1
     2                R0 ,RP5,R0 ,RP5,R0 ,R0 ,R0 ,R0 ,R0 , ! 2
     3                R0 ,R0 ,RP5,R0 ,R0 ,R0 ,RP5,R0 ,R0 , ! 3
     4                R0 ,RP5,R0 ,RP5,R0 ,R0 ,R0 ,R0 ,R0 , ! 4
     5                R0 ,R0 ,R0 ,R0 ,R1 ,R0 ,R0 ,R0 ,R0 , ! 5
     6                R0 ,R0 ,R0 ,R0 ,R0 ,RP5,R0 ,RP5,R0 , ! 6
     7                R0 ,R0 ,RP5,R0 ,R0 ,R0 ,RP5,R0 ,R0 , ! 7
     8                R0 ,R0 ,R0 ,R0 ,R0 ,RP5,R0 ,RP5,R0 , ! 8
     9                R0 ,R0 ,R0 ,R0 ,R0 ,R0 ,R0 ,R0 ,R1/), (/9,9/)) ! 9
C                      1   2   3   4   5   6   7   8   9
C Assemble deviatoric projection tensor (use G matrix ordering)
CC FOIDS3D(I,J)-R1D3*SOID3D(I)*SOID3D(J)
C  1/2*(delta_ik*delta_jl + delta_il*delta_jk) - 1/3*(delta_ij*delta_kl)
      DOUBLE PRECISION, PARAMETER, DIMENSION(4,4) :: 
     1            DEVPR2=reshape((/R2D3 ,R0 ,R0 ,-R1D3, ! 1
     2                    R0   ,RP5,RP5,R0   , ! 2
     3                    R0   ,RP5,RP5,R0   , ! 3
     4                    -R1D3,R0 ,R0 ,R2D3 /), (/4,4/)) ! 4
C                         1     2   3   4
      DOUBLE PRECISION, PARAMETER, DIMENSION(9,9) :: 
     1   DEVPR3=reshape((/R2D3 ,R0 ,R0 ,R0 ,-R1D3,R0 ,R0 ,R0 ,-R1D3, ! 1
     2                    R0   ,RP5,R0 ,RP5,R0   ,R0 ,R0 ,R0 ,R0   , ! 2
     3                    R0   ,R0 ,RP5,R0 ,R0   ,R0 ,RP5,R0 ,R0   , ! 3
     4                    R0   ,RP5,R0 ,RP5,R0   ,R0 ,R0 ,R0 ,R0   , ! 4
     5                    -R1D3,R0 ,R0 ,R0 ,R2D3 ,R0 ,R0 ,R0 ,-R1D3, ! 5
     6                    R0   ,R0 ,R0 ,R0 ,R0   ,RP5,R0 ,RP5,R0   , ! 6
     7                    R0   ,R0 ,RP5,R0 ,R0   ,R0 ,RP5,R0 ,R0   , ! 7
     8                    R0   ,R0 ,R0 ,R0 ,R0   ,RP5,R0 ,RP5,R0   , ! 8
     9            -R1D3,R0 ,R0 ,R0 ,-R1D3,R0 ,R0 ,R0 ,R2D3 /), (/9,9/)) ! 9
C                 1     2   3   4   5     6   7   8   9
      
C**********************************************************************
C COMPUTE CONSISTENT SPATIAL TANGENT MODULUS FOR THE ANISOTROPIC 
C VISCO-PLASTIC SINGLE CRYSTAL MODEL (PLANE STRAIN AND 3D 
C IMPLEMENTATIONS)
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
C Retrieve elastic deformation gradient
      FE=RESHAPE(RSTAVA(1:9),[MDIM,MDIM])
      DETFE=DETM23(3, FE, NDIM, .FALSE.)
C Retrieve material properties
      GMODU=RPROPS(2) ! Neo-hookean shear and bulk moduli
      BULK=RPROPS(3)
      MODEL=IPROPS(19) ! Viscoplastic slip law choice
C Store Cauchy stresses in matrix form
      CALL ATSYM( STRES, SIGMA, NDIM, .TRUE.)
C
C ELASTIC tangent modulus
C =======================
C (note: Norton viscoplastic slip law (MODEL==1) is always inelastic,
C  even for unloading - it doesn't have yield surfaces)
      IF((MODEL/=1).AND.(.NOT.EPFLAG))THEN
C Deviatoric Cauchy stresses stored in array form (G-matrix ordering)
        DEVSTR(1:NGDIM)=RESHAPE(SIGMA(1:NDIM,1:NDIM), [NGDIM])
        P=BULK*LOG(DETFE)/DETFE ! Hydrostatic pressure
        IF(IS2D)THEN
          DEVSTR([1,4])=DEVSTR([1,4])-P
        ELSE
          DEVSTR([1,5,9])=DEVSTR([1,5,9])-P
        ENDIF
C Trace of isochoric component of elastic left Cauchy-Green strain
C tensor (TRBISO = BEISO_ii = FEISO_ik FEISO_ik = J^(-2/3) FE_ik FE_ik)
        TRBISO=DETFE**(-R2D3) * SUM(FE*FE)
        GFAC=R2D3*GMODU*TRBISO/DETFE
        BULFAC=BULK/DETFE
        R2P=R2*P
C Assemble tensorially compact part of AMATX
        IF(IS2D)THEN
          DO I=1,NGDIM
            DO J=1,NGDIM
              AMATX(I,J)=BULFAC*SOID2(I)*SOID2(J)
     1                   -R2P*FOIDS2(I,J)+GFAC*DEVPR2(I,J)
     2                   -R2D3*(DEVSTR(I)*SOID2(J)+SOID2(I)*DEVSTR(J))
            END DO
          END DO
        ELSE
          DO I=1,NGDIM
            DO J=1,NGDIM
              AMATX(I,J)=BULFAC*SOID3(I)*SOID3(J)
     1                   -R2P*FOIDS3(I,J)+GFAC*DEVPR3(I,J)
     2                   -R2D3*(DEVSTR(I)*SOID3(J)+SOID3(I)*DEVSTR(J))
            END DO
          END DO
        ENDIF
C Add non-compact part: delta_ik sigma_jl
C AMATX_ijkl = AMATX_ijkl + delta_ik * sigma_jl = AMATX_ijil + sigma_jl
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
C ELASTOPLASTIC tangent modulus
C =============================
      ELSE
C... recover current plastic deformation gradient
        FP=RESHAPE(RSTAVA(10:18),[MDIM,MDIM])
        F=MATMUL(FE, FP)
        CALL INVMT3(FP, FPINV, DETFP)
C... previous converged elastic deformation gradient (equilibrium)
        FEN=RESHAPE(RSTAVN(1:9),[MDIM,MDIM])
C Set current hardening internal variable
        HRVAR=RSTAVA(19)
C... number of sampling points on hardening curve
        NHARD=IPROPS(3)
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
      SIGNSC=SIGN(R1, SCHMID) ! Sign of Schmid stresses
C Current resolved reference stress / yield shear stress
      RSHEAR=PLFUN(HRVAR,NHARD,RPROPS(IPHARD)) !Norton: resolved reference stress
C                                              !Peric: yield resolved shear stress
C??? add checking if RSHEAR<R0 (like in ELFEM routines)
      OVERST=ABS(SCHMID)/RSHEAR
C
      IF(MODEL==1)THEN
        ISOVER=OVERST > SMALL ! No yield surface!
      ELSE
        ISOVER=OVERST-R1 > SMALL
      ENDIF
C
      HSLOPE=DPLFUN(HRVAR,NHARD,RPROPS(IPHARD))
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
     1                                         +VECM(I,:)*VECS0(J,:) )
        ENDDO
      ENDDO
        
      DO I=1,NDIM
        DO J=1,NDIM
          FPILOG(I,J)=-SUM(DG(:)*S0M0(I,J,:), MASK=ISOVER)
        ENDDO
      ENDDO
        
      CALL DEXPMP(   DEREXP      ,NOCONV     ,FPILOG     )
      IF(NOCONV)THEN
C Derivative of exponential map failed: Break loop and exit
        STOP 'Derivative exp map failed - CTVCPL'
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
C Derivative of residual of isochoric elastic def.grad. update with
C respect to TOTAL isochoric def.grad.
C -----------------------------------------------------------------
      A4TH=R0
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
      DO I=1,NDIM
        DO J=1,NDIM
          DO K=1,NDIM
            DO L=1,NDIM
              AUXM14(I,J,K,L)=DF1DA(I,J)*DF2DFE(K,L)/DF2DA
     1                        -DF1DF4(I,J,K,L)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      IF(IS2D)THEN
        CALL ARRGO2
     1(   AUXM14(1:NDIM,1:NDIM,1:NDIM,1:NDIM) ,AUXM1(1:NGDIM,1:NGDIM) )
      ELSE
        CALL ARRGO3
     1(   AUXM14 ,AUXM1 )
      ENDIF
      AUXM2=R0
      CALL RMINVE
     1( AUXM1(1:NGDIM,1:NGDIM) ,AUXM2(1:NGDIM,1:NGDIM) ,NGDIM ,IFGERR )
      IF(IFGERR)THEN
        STOP 'RMINVE failed! CSTVPC'!CALL ERRPRT('WE0031')
      ENDIF
      DFEDF=R0
      DFEDF=MATMUL(AUXM2, DF1DF)
      IF(IS2D)THEN
        DFEDF(5,5)=R1
      ENDIF
      
C Compute C matrix
C ================
C First as 4th order tensor
      DO I=1,MDIM
        DO J=1,MDIM
          DO K=1,MDIM
            DO L=1,MDIM
              C4TH(I,J,K,L)=DELTA(I,K)*FEISO(J,L)+FEISO(I,L)*DELTA(J,K)
     1                      -R2D3*DELTA(I,J)*FEISO(K,L)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C...then store in matrix form (G-matrix component ordering)
      IF(IS2D)THEN
        CMATX=R0
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
        CALL ARRGAX
     1(   H4TH  ,HMATX(1:5,1:5) )
      ELSE
        CALL ARRGO3
     1(   H4TH      ,HMATX       )
      ENDIF
C Assemble spatial tangent "a"
C ============================
      FACTOR=BULK/DETFE
      A4TH=R0
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
      AUXM25=R0
      AUXM25=MATMUL(MATMUL(CMATX, DFEDF), HMATX)
C... finally get "a" matrix
      FACTOR=GMODU/(DETFE**(R4D3))
      AMATX=FACTOR*AUXM25+AUXM1
      ENDIF
      
      RETURN
      END
CDOC END_SUBROUTINE CSTVPC
