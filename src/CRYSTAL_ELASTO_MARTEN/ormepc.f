CDOC BEGIN_SUBROUTINE ORMEPC
CDOC Output results for martensitic transformation elasto-plastic
CDOC single crystal model
CDOC
CDOC This routine writes to the results file the internal and
CDOC algorithmic variables of the martensitic transformation elasto-
CDOC plastic single crystal material model.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DGAM   >  Array of incremental plastic
CDOC C                          multipliers.
CDOC C                          Computed in routine SUTR.
CDOC INTEGER          NOUTF  >  Results file unit identifier.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC DOUBLE_PRECISION RPROPS >  Array of real material properties.
CDOC DOUBLE_PRECISION RSTAVA >  Array of real state variables other than
CDOC C                          the stress tensor components.
CDOC DOUBLE_PRECISION STRES  >  Array of stress tensor components.
CDOC END_PARAMETERS
CHST
CHST M. Fauzan Adziman, April 2013: Initial coding
CHST
      SUBROUTINE ORMEPC
     1(   DGAM       ,NOUTF      ,NTYPE      ,RPROPS     ,RSTAVA     ,
     2    STRES      ,IPROPS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER
     1(   NDIM=2     ,NRSTAV=52  ,NFDIM=3    , NIPROP=19)
C Arguments
      DIMENSION
C1     1    DGAM(NTSYST)        ,
     1    RPROPS(*)          ,RSTAVA(NRSTAV)     ,
     2    STRES(*)           ,IPROPS(NIPROP)
C Local arrays and variables
      DIMENSION
C     1    FE(NDIM,NDIM)      ,RE(NDIM,NDIM)      ,
     1    FE(NFDIM,NFDIM)    ,RE(NDIM,NDIM)      ,
     2    UE(NDIM,NDIM)      ,FEISO(NFDIM,NFDIM) ,BEISO(NFDIM,NFDIM) ,
     3    BEDEV(4)           ,CAUCHE(4)          ,TN(NFDIM,NFDIM)    ,
     4    PAUXTR(2,1)        !,TNMD(1)            ,
C     5    VECS0(3,1)         ,VECM0(3,1)         ,VECD0(3,1)         ,
C1     5    GAM(1)             ,D(1)
      
      

      
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C #DEBUG-BEG#
C EXTRA VARIABLES FOR VISCOPLASTIC SYSTEMS
C Number of slip systems: IPROPS(4) -> 2*IPROPS(4) using the double slip
C convention
C      INTEGER, PARAMETER :: NTRSYS=4
      DOUBLE PRECISION
     1     VECM0(3,IPROPS(4)), VECS0(3,IPROPS(4)),
     1     VECM(3,IPROPS(4)), VECS(3,IPROPS(4)),
     1     S0M0(3,3,IPROPS(4)),
     2     DKIRCH(NDIM,NDIM)  ,ROTSTR(NDIM,NDIM)  ,
     3     SCHMID(2*IPROPS(4)),
     4     HVECM0(3,IPROPS(5)),HVECD0(3,IPROPS(5)),HVD0M(3,3,IPROPS(5)),
     5     DT(IPROPS(5))
C #DEBUG-END#
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      
      
      DATA R0   ,R1   ,R3   ,R180   /
     1     0.0D0,1.0D0,3.0D0,180.0D0/
C***********************************************************************
C OUTPUT RESULTS (INTERNAL AND ALGORITHMIC VARIABLES) FOR THE ANISO-
C TROPIC TRANSFORMATION SINGLE CRYSTAL ELASTOPLASTIC MODEL
C***********************************************************************
C~C Full output option
C~C ------------------
C~ 1000 FORMAT(' Clrot = ',G12.4, ' Gamma = ',G12.4, ' D     = ',G12.4/ 
C~     2       ' Sxx-e = ',G12.4, ' Syy-e = ',G12.4, ' Sxy-e = ',G12.4, 
C~     3       ' Szz-e = ',G12.4/
C~     4       ' T-11  = ',G12.4, ' T-21  = ',G12.4, ' T-12  = ',G12.4,
C~     5       ' T-22  = ',G12.4, ' T-33  = ',G12.4/
C~     6       ' Dvolu = ',G12.4)
C~C Simplified output option
C~C ------------------------  
C~ 1000 FORMAT(' Clrot = ',G12.4, ' Gamma = ',G12.4, ' D     = ',G12.4,    
C~     1       ' Dvolu = ',G12.4)   
      
C Stops program if not plane strain
C!      IF(NTYPE.NE.2)CALL ERRPRT('EI0037')
      
C Retrieve material and system properties
C ---------------------------------------
      R1D3=R1/R3
      GMODU=RPROPS(2)
      BULK=RPROPS(3)
C Input parameters for the transformation
C ---------------------------------------
C Recover systems information from RPROPS
      NSLSYS=IPROPS(4)
      NTRSYS=IPROPS(5)
      
      HVECM0=RESHAPE(RPROPS(IPROPS(6):IPROPS(7)),[3,NTRSYS])
      HVECD0=RESHAPE(RPROPS(IPROPS(8):IPROPS(9)),[3,NTRSYS])
      HVD0M=RESHAPE(RPROPS(IPROPS(10):IPROPS(11)),[3,3,NTRSYS])
      VECM0=RESHAPE(RPROPS(IPROPS(12):IPROPS(13)),[3,NSLSYS])
      VECS0=RESHAPE(RPROPS(IPROPS(14):IPROPS(15)),[3,NSLSYS])
      S0M0=RESHAPE(RPROPS(IPROPS(16):IPROPS(17)),[3,3,NSLSYS])
C      
C Retrieve stored state variables
C -------------------------------
C... current elastic deformation gradient
      FE=RESHAPE(RSTAVA(1:9),[3,3])
C Compute lattice rotation
C ------------------------
C Perform polar decomposition of the elastic deformation gradient
      CALL PODEC2
C1     1(   FE         ,RE         ,UE         )
     1(   FE(1:2,1:2)  ,RE         ,UE         )
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
C Compute elastic Cauchy stresses (not stored in memory)
C ------------------------------------------------------
C Perform isochoric/volumetric split of elastic deformation gradient
      DETFE=FE(1,1)*FE(2,2)-FE(1,2)*FE(2,1)
      VOLFAC=DETFE**(-R1D3)
      CALL RVZERO(FEISO,NFDIM*NFDIM)
      FEISO(1,1)=VOLFAC*FE(1,1)
      FEISO(2,1)=VOLFAC*FE(2,1)
      FEISO(1,2)=VOLFAC*FE(1,2)
      FEISO(2,2)=VOLFAC*FE(2,2)
      FEISO(3,3)=VOLFAC*R1 
C Compute elastic left Cauchy-Green tensor
      CALL RVZERO(BEISO,9)
      DO 30 I=1,NDIM
        DO 20 J=1,NDIM
          DO 10 K=1,NDIM
            BEISO(I,J)=BEISO(I,J)+FEISO(I,K)*FEISO(J,K)
   10     CONTINUE
   20   CONTINUE
   30 CONTINUE
      BEISO(3,3)=VOLFAC*VOLFAC
C Hydrostatic pressure
      P=BULK*LOG(DETFE)
C Deviatoric component of isochoric elastic left Cauchy-Green tensor
      TRACE=BEISO(1,1)+BEISO(2,2)+BEISO(3,3)
      BEDEV(1)=BEISO(1,1)-R1D3*TRACE
      BEDEV(2)=BEISO(2,2)-R1D3*TRACE
      BEDEV(3)=BEISO(1,2)
      BEDEV(4)=BEISO(3,3)-R1D3*TRACE
C Update Cauchy stress components as a function of elastic components
      DETINV=R1/DETFE
      CAUCHE(1)=(GMODU*BEDEV(1)+P)*DETINV
      CAUCHE(2)=(GMODU*BEDEV(2)+P)*DETINV
      CAUCHE(3)=(GMODU*BEDEV(3))*DETINV
      CAUCHE(4)=(GMODU*BEDEV(4)+P)*DETINV 
C!!!! ??? -> see if anything equivalent exists in any crystal models output
C	  
C Reset gamma if exceeds 1 implying transformation has completed 
C --------------------------------------------------------------
C... retrive current value of transformation multiplier
      GAMMA=RSTAVA(19)
      IF(GAMMA.GE.R1)THEN
        GAMMA=R1 
      ENDIF
C
C Calculate dissipation
C ---------------------
C... retrieve current work conjugate of plastic deformation
      TN=RESHAPE(RSTAVA(20:28),[3,3])
C... calculate work conjugate of plastic deformation times vector m
C    and internal product with vector of the habit plan system (d)
      CALL RVZERO(PAUXTR,2)
      DO I=1,NDIM
        DO J=1,NDIM   
          PAUXTR(I,1)=PAUXTR(I,1)+TN(I,J)*HVECM0(J,1)
        END DO
      END DO
      TNMD=PAUXTR(1,1)*HVECD0(1,1)+PAUXTR(2,1)*HVECD0(2,1)
C... calculate instantaneous dissipated mechanical work 
      DISS=GAMMA*TNMD

      

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C #DEBUG-BEG#

C Full output option
C ------------------
 1000 FORMAT(' Sxx-e = ',G12.4, ' Syy-e = ',G12.4, ' Sxy-e = ',G12.4, 
     1       ' Szz-e = ',G12.4/
     2       ' T-11  = ',G12.4, ' T-21  = ',G12.4, ' T-12  = ',G12.4,
     3       ' T-22  = ',G12.4, ' T-33  = ',G12.4/
     4       ' Clrot = ',G12.4, ' Gamma = ',G12.4, ' D     = ',G12.4,
     5       ' Dvolu = ',G12.4/ 
     6       ' hrvar = ',G12.4, ' Gam*DGC=',G12.4, ' iacttr= ',I2)
 1100 FORMAT(' tau_',I0.2,'= ',G12.4)
 1200 FORMAT('  D_',I0.2,' = ',G12.4)
      
      
C EXTRA VARIABLES RELATED TO VISCO MODEL
      HRVAR=RSTAVA(41)
C Evaluate resolved Schmid stresses (these are not stored in memory)
C ------------------------------------------------------------------
      NSYST=IPROPS(4)
C1 compute rotated stress tensor
C1... deviatoric Kirchhoff stress
C1      PKIRCH=R1D3*DETFE*(STRES(1)+STRES(2)+STRES(4))
C1      DKIRCH(1,1)=DETFE*STRES(1)-PKIRCH
C1      DKIRCH(2,2)=DETFE*STRES(2)-PKIRCH
C1      DKIRCH(1,2)=DETFE*STRES(3)
C1      DKIRCH(2,1)=DKIRCH(1,2)
C1C... rotated deviatoric Kirchhoff stress
C1      ROTSTR=MATMUL(MATMUL(TRANSPOSE(RE),DKIRCH),RE) ! Re^T Tau_d Re
C1      DO ISYST=1,NSYST
C1        SCHMID(ISYST)=SUM(ROTSTR(:,:)*S0M0(1:NDIM,1:NDIM,ISYST)) !C systems 1...NSYST
C1        SCHMID(ISYST+NSYST)=-SCHMID(ISYST) !C systems NSYST+1...2*NSYST (same normal, reverse shear direction -> S0M0 changes sign)
C1      ENDDO
      
      
C using formula specific for Neo-Hookean
      DO ISYST=1,NSYST
        VECS(:,ISYST)=MATMUL(FEISO, VECS0(:,ISYST))
        VECM(:,ISYST)=MATMUL(FEISO, VECM0(:,ISYST))
        SCHMID(ISYST)=GMODU*DOT_PRODUCT(VECS(:,ISYST),VECM(:,ISYST))
        SCHMID(ISYST+NSYST)=-SCHMID(ISYST)
      ENDDO
      
      
C!!! add output of CURRENT AND ORIGINAL orientation of slip vectors / martensitic transf. vectors
      
      
      
C EXTRA VARIABLES RELATED TO MARTENSITIC TRANSFORMATION
C Calculate dissipation
C ---------------------
C... calculate work conjugate of plastic deformation times vector m
C    and internal product with vector of the habit plan system (d)
      DO ISYST=1,NTRSYS
        DT(ISYST)=GAMMA*
     1                 SUM(TN(1:NDIM,1:NDIM)*HVD0M(1:NDIM,1:NDIM,ISYST))
      ENDDO
C Active transformation system
      IACTRS=INT(RSTAVA(39))
      DGCR=RPROPS(5)
C
C Write results to output file
C ----------------------------
C Full output option
C ------------------
      WRITE(NOUTF,1000)CAUCHE(1:4),RSTAVA([20,21,23,24,28]), !(1,1),(2,1),skip (3,1),(1,2),(2,2),skip (3,2), skip (1,3), skip (2,3), (3,3)
     2            CLROT,GAMMA,DISS,RSTAVA(38),HRVAR,GAMMA*DGCR,IACTRS
CC Schmid stress in viscoplastic crystal slip systems
C      DO ISYST=1,2*NSYST
C        WRITE(NOUTF,1100)ISYST,SCHMID(ISYST)
C      ENDDO
CC Dissipation in martensitic transformation systems
C      DO ISYST=1,NTRSYS
C        WRITE(NOUTF,1200)ISYST,DT(ISYST)
C      ENDDO
      
      
      
C #DEBUG-END#
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      
C      
C~C Full output option
C~C ------------------
C~      WRITE(NOUTF,1000)CLROT,GAM(1),D(1),CAUCHE(1),CAUCHE(2),
C~     1                 CAUCHE(3),CAUCHE(4),RSTAVA(10),RSTAVA(11),
C~     2                 RSTAVA(12),RSTAVA(13),RSTAVA(14),RSTAVA(19)
C~C Simplified output option
C~C ------------------------
C~      WRITE(NOUTF,1000)CLROT,GAM(1),D(1),RSTAVA(19)
C 
      RETURN
      END
CDOC END_SUBROUTINE ORMEPC
