CDOC BEGIN_SUBROUTINE SUMTSC
      SUBROUTINE SUMTSC
     1(   DGAMMA     ,FINCR      ,IPROPS     ,LALGVA     ,NTYPE      ,
     2    RPROPS     ,RSTAVA     ,STRES      ,DVOLU      ,IELEM      ,
     3    IGAUSP     ,DTIME      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)	
      PARAMETER
     1(   NIPROP=19  ,NLALGV=6   ,NRSTAV=52  ,NSTRE=6    )
C Arguments
      LOGICAL
     1    LALGVA(NLALGV)
      DIMENSION
     1    FINCR(3,3)         ,IPROPS(NIPROP)     ,
     2    STRES(NSTRE)       ,RPROPS(*)          ,RSTAVA(NRSTAV)
C Local arrays and variables for martensitic transformation
      LOGICAL
     1    IFTRAN    ,SUFAIL  ,IFPLAS    ,IFTRN2, IFPLUP
      DIMENSION
     1    DPHIA(3,3)         ,DPHIB(3,3)         ,DTDFE(3,3,3,3)
      DIMENSION
     1    FE(3,3)            ,FTOTA(3,3)         ,FTRA(3,3)  ,
     2    FETRL(3,3)         ,FTOINV(3,3)        ,DFEDDG(3,3),
     3    FEN(3,3)           ,FTRAN(3,3)         ,FTRINV(3,3)
      DIMENSION
     2    HVD0M(3,3,IPROPS(5))
      DIMENSION
     2    PHIVAR(IPROPS(5))  ,PHIAUX(IPROPS(5))  ,
     3    PHIASC(IPROPS(5))  ,TAUC(3,3)          ,TNTRL(3,3)         ,
     4    TN(3,3)            ,TAUTRC(3,3)
      DIMENSION
     1    FPAUSN(3,3)        ,FINEL(3,3)         ,FININV(3,3)
      
      DIMENSION FTOTA2(3,3), FPAUS2(3,3)
      
      
      DOUBLE PRECISION FTRANTEMP(3,3), PIOLA(3,3), FTOTN(3,3), 
     1                 SIGMA(3,3), FTOTINV(3,3)
      
      DOUBLE PRECISION DTEMP(9)
      
      
      DOUBLE PRECISION, SAVE :: DISSIP=0.0D0
      
      
      DOUBLE PRECISION GAMMA(IPROPS(5)), GAMMAD(IPROPS(5)),
     1    RESVEC(IPROPS(5)), DREMTX(IPROPS(5),IPROPS(5)), 
     2    TRSTRE(IPROPS(5)), GAMMAP(IPROPS(5)), GAMMADP(IPROPS(5)),
     3    FTRAP(3,3), FINELP(3,3), FININVP(3,3), FEP(3,3), TAUP(3,3),
     4    TP(3,3), TRSTREP(IPROPS(5)), RESVECP(IPROPS(5))
      
      LOGICAL ISTRAN(IPROPS(5)), ISTRANP(IPROPS(5)), SINGUL
      
      
      
C Local parameters
      DOUBLE PRECISION, PARAMETER :: R0=0.0D0     , R1=1.0D0        ,
     1                               TOL=1.0D-8
      INTEGER, PARAMETER          :: MXITER=50
      DOUBLE PRECISION, PARAMETER, DIMENSION(3,3) :: ! delta_ij
     1     DELTA=reshape((/R1,R0,R0,R0,R1,R0,R0,R0,R1/), (/3,3/) )
C***********************************************************************

C Check analysis type is valid for this model
      IF(NTYPE==2)THEN
        NDIM=2
      ELSEIF(NTYPE==4)THEN
        NDIM=3
      ELSE
        CALL ERRPRT('EI0080')
      ENDIF
C Initialise some algorithmic and internal variables
      DGAMMA=R0
      IFTRAN=.FALSE.
      IFPLAS=.FALSE.
      SUFAIL=.FALSE.
      
C7
C Retrieve previous converged transformation state
      IFTRN2=LALGVA(1)
      IFPLUP=.FALSE.
      
C
CC1      GAMMAN=RSTAVA(19)
CC1      FTOTN=RESHAPE(RSTAVA(29:37),[3,3])
C
      
C7
      
C First, perform plastic return-mapping (only if transformation has not already started)
C      IF(.NOT.IFTRN2)THEN
C      IF(.NOT.IFTRN2.AND.GAMMA==R0)THEN
CC1      GAMMA=RSTAVA(19) ! Previous state of transformation multiplier
      GAMMA=RSTAVA(20:43)
      GAMMAT=SUM(GAMMA)
CC2
      IF(GAMMAT==R0)THEN
CC2
C      IF(GAMMAT<R1)THEN
CC2
C Check plastic admissibility 
C ===========================
C... and if PLASTIC then compute the viscoplastic contribution.         
C5        IF(CKPLAS.AND.IFPRTR.EQ.R0)THEN
C5        CALL SUVCPL
C5     1    ( DTIME      ,FINCR      ,IPROPS     ,LALGVA     ,NTYPE     ,
C5     2      RPROPS     ,RSTAVA     ,STRES      )
        IFPLUP=.TRUE.
CC1        DTEMP=RSTAVA(10:18)
CC1        RSTAVA(10:18)=RSTAVA(42:50) ! copy FP to its place
CC1        DTEMP2=RSTAVA(19)
CC1        RSTAVA(19)=RSTAVA(41) ! copy HRVAR
        
        
CC2        CALL SUVSC2
CC2     1(   DTIME      ,FINCR      ,IPROPS     ,LALGVA     ,NTYPE      ,
CC2     2    RPROPS     ,RSTAVA     ,STRES      )
        CALL SUVSC2MOD
     1(   DTIME      ,FINCR      ,IPROPS     ,LALGVA     ,NTYPE      ,
     2    RPROPS     ,RSTAVA     ,STRES      )
CC2
      
      
      
CC1        RSTAVA(42:50)=RSTAVA(10:18) ! new FP
CC1        RSTAVA(41)=RSTAVA(19) ! new HRVAR
CC1        RSTAVA(10:18)=DTEMP
CC1        RSTAVA(19)=DTEMP2
        IFPLAS=LALGVA(1)
        SUFAIL=LALGVA(2)
C
C5        SUFAIL=LALGVA(2)
C5        IFPLAS=LALGVA(3)
        IF(SUFAIL)GOTO 999
      ENDIF
      
      
C
C Retrieve material properties:
C -----------------------------
C Neo-Hookean shear and bulk moduli
      GMODU=RPROPS(2)
      BULK=RPROPS(3)
C Martensitic transformation properties
      DGCR=RPROPS(5) ! Energy barrier (required work to start transformation)
      NTRSYS=IPROPS(5) ! Number of transformation systems
      HVD0M=RESHAPE(RPROPS(IPROPS(10):IPROPS(11)),[3,3,NTRSYS]) ! System constant matrix: d (x) m
C Martensitic transformation state variables
CC1      GAMMA=RSTAVA(19) ! Previous state of transformation multiplier
CC1      IACTRS=INT(RSTAVA(39)) ! Activated transformation system
C7      IFVARS=INT(RSTAVA(40)) ! Identify whether transformation system selection is required
C7      IFPRTR=INT(RSTAVA(51)) ! State of transformation in the previous converged solution
C
C Trial deformation gradients:
C ----------------------------
C Retrieve previous elastic, austenite plastic and transformation
C deformation gradients
      FEN=RESHAPE(RSTAVA(1:9),[3,3])
CC1      FPAUSN=RESHAPE(RSTAVA(42:50),[3,3])
      FPAUSN=RESHAPE(RSTAVA(10:18),[3,3])
CC1      FTRAN=RESHAPE(RSTAVA(10:18),[3,3])
      FTRAN=DELTA
      DO I=1,NTRSYS
        FTRAN=FTRAN + GAMMA(I)*HVD0M(:,:,I)
      ENDDO
      
CC1
C Extra transformation properties
      TRMU=RPROPS(6)
      TREPS=RPROPS(7)
CC1
      
C6      IF(IACTRS==0)THEN
C6        FTRAN=DELTA
C6      ELSE
C6        FTRAN=DELTA+GAMMA*HVD0M(:,:,IACTRS)
C6      ENDIF
C6C  -> FTRAN can be recovered from gamma and IACTRS
C6C      IF(IACTRS==0)THEN
C6C        FTRANTEMP=DELTA
C6C      ELSE
C6C        FTRANTEMP=DELTA+GAMMA*HVD0M(:,:,IACTRS)
C6C      ENDIF
C6C      IF(GAMMA/=R0)THEN
C6C        WRITE(*,*)'FTRAN, FTRANTEMP'
C6C        DO I=1,3
C6C          WRITE(*,'(6g20.6)')FTRAN(I,:), FTRANTEMP(I,:)
C6C        ENDDO
C6C      ENDIF
      
C7C Current elastic trial deformation gradient
C7      FETRL=MATMUL(FINCR,FEN)
      IF(.NOT.IFPLUP)THEN
        FETRL=MATMUL(FINCR,FEN)
      ELSE ! Plastic update already happened, trial is whatever is left of elastic def grad?
        FETRL=FEN
      ENDIF
C Current total deformation gradient and its inverse
      FTOTA=MATMUL(FETRL,MATMUL(FTRAN,FPAUSN))
C!! is this different than FINCR*previous FTOTA ( from RSTAVA) ????
      CALL INVMT3(FTOTA, FTOINV, DETFTO)
C
C Trial work conjugate of transformation deformation gradient (TNTRL):
C ------------------------------------------------------------
      CALL TAUMEP(FETRL, NDIM, RPROPS, TAUTRC) ! Trial Kirchhoff stress
C8      
C8      TNTRL=MATMUL(TRANSPOSE(FETRL),MATMUL(TAUTRC,TRANSPOSE(FTOINV)))
C8      
      TNTRL=MATMUL(MATMUL(TRANSPOSE(FETRL),
     1             MATMUL(TAUTRC,TRANSPOSE(FTOINV))),TRANSPOSE(FPAUSN))
C8
      
C!!!!!!!! if transformation finished, can skip to elastic update
        IF(GAMMAT>=R1)THEN
c          WRITE(*,*)'WARNING: gamma >= 1 in SUMTSC.Skipping transf. R-M'
c          WRITE(*,*)'Gamma total:', GAMMAT
c          WRITE(*,*)GAMMA
c          READ(*,*)
          GOTO 666
        ENDIF
C!!!!!!!!
C!! create logical variable for transformation finished (ISMART ?)
      
      
      
C
C Check martensitic transformation admissibility
C ==============================================
C Compute transformation function value at trial state
C ----------------------------------------------------
      
      
C!!!
C If transformation system has not already been chosen, select the most
C favourable one (highest transformation function value)
C !!!! can remove variable IFVARS
CC1        IF(IACTRS==0)THEN !variant selection still needed
CC1        DO ISYST=1,NTRSYS
CC1          PHIVAR(ISYST)=SUM(TNTRL*HVD0M(:,:,ISYST))-DGCR
CC1        ENDDO
CC1C5        IACTRS=MAXLOC(PHIVAR,1)
CC1C5        PHI=PHIVAR(IACTRS)
CC1        PHI=MAXVAL(PHIVAR,1)
CC1      ELSE
CC1        PHI=SUM(TNTRL*HVD0M(:,:,IACTRS))-DGCR
CC1      ENDIF
CC1        DO ISYST=1,NTRSYS
CC1          PHIVAR(ISYST)=SUM(TNTRL*HVD0M(:,:,ISYST))-DGCR
CC1        ENDDO
        DO ISYST=1,NTRSYS
          TRSTRE(ISYST)=SUM(TNTRL*HVD0M(:,:,ISYST))
        ENDDO
CC1
C!!!
      
C5      DO ISYST=1,NTRSYS
C5        PHIVAR(ISYST)=SUM(TNTRL*HVD0M(:,:,ISYST))-DGCR
C5      ENDDO
C5      PHIAUX=PHIVAR/DGCR-TOL
C5      PHIASC=PHIAUX
C5      
C5      
C5C... select the most favourable variant for transformation             
C5      IF(IFVARS.EQ.0)THEN
C5C... sort variants (in an ascending order) according to their potential   
C5C    to initiate a transformation 
C5        
C5C2        WRITE(*,*)'====== START SELECTION====='
C5C2        WRITE(*,'(*(F15.6))')PHIASC
C5        
C5          DO I=1,NTRSYS-1
C5            DO J=I+1,NTRSYS
C5              IF(PHIAUX(I).GT.PHIAUX(J))THEN
C5                 AUXLIST=PHIASC(I)
C5                 PHIASC(I)=PHIASC(J)
C5                 PHIASC(J)=AUXLIST
C5                 
C5C2        WRITE(*,*)'--------swap------------'
C5C2        WRITE(*,'(*(F15.6))')PHIASC
C5        
C5                 
C5              ENDIF
C5            END DO
C5          END DO
C5C... The bottom data in the sorted array is the most positive value, 
C5C    hence the most favoured variant to transform
C5          PHIDEL=PHIASC(NTRSYS)      
C5        DO ISYST=1,NTRSYS
C5          IF(PHIDEL.EQ.PHIAUX(ISYST))THEN
C5            IACTRS=ISYST
C5            PHI=PHIVAR(IACTRS)
C5          ENDIF	    
C5        END DO
C5        
C5        
C5C2        WRITE(*,*)'--------sys sel------------'
C5C2        WRITE(*,'(*(F15.6))')PHIVAR
C5C2        WRITE(*,'(*(F15.6))')PHIAUX
C5C2        WRITE(*,'(*(F15.6))')PHIASC
C5C2        WRITE(*,'(*(F15.6))')PHIDEL
C5C2        WRITE(*,'(I6,F15.6)')IACTRS, PHI
C5C2        WRITE(*,'(I6,F15.6)')MAXLOC(PHIVAR), PHIVAR(MAXLOC(PHIVAR))
C5C2        WRITE(*,'(I6,F15.6)')MAXLOC(PHIAUX), PHIVAR(MAXLOC(PHIAUX))
C5        
C5        
C5C... Once a transformation has been initiated at a gauss point, the 
C5C    selected variant at the gauss point to be preserved throughout
C5C    the transformation process (see also SWMEPC sub-routine)
C5      ELSE
C5        PHI=PHIVAR(IACTRS)
C5      ENDIF
      
      
C... check transformation consistency
CC1      IF(PHI/DGCR.GT.TOL)THEN
CC1        IFTRAN=.TRUE.
CC1C If transforming variant not yet defined, select it
CC1        IF(IACTRS==0)THEN
CC1          IACTRS=MAXLOC(PHIVAR,1)
CC1        ENDIF
CC1      ENDIF
      ISTRAN=TRSTRE > DGCR
      IF(ANY(ISTRAN))IFTRAN=.TRUE.
                   
C Transformation step: Apply return-mapping
C =========================================        
C Setting initial states at the beginning of the return-mapping loop 
C ------------------------------------------------------------------
      IF(IFTRAN)THEN
        FE=FETRL
        FTRA=FTRAN
        TAUC=TAUTRC
        
CC1 Initial guess for delta_gamma: 0
        GAMMAD=R0
        
CC1
        
C
C Start Newton-Raphson iterations to find solution of the return-
C mapping equation as a function of transformation multiplier 
        nr: DO ITER=1,MXITER
C Compute Jacobian matrix of the return-mapping (DPHI)
C ---------------------------------------------
    
        
C Calculate residual vector
C==========================
C Update all variables: 
          FTRA=DELTA
          DO I=1,NTRSYS
            FTRA=FTRA + GAMMA(I)*HVD0M(:,:,I)
          ENDDO
          FINEL=MATMUL(FTRA,FPAUSN)
          CALL INVMT3(FINEL, FININV, DETFIN)
          FE=MATMUL(FTOTA,FININV)
          CALL TAUMEP(FE, NDIM, RPROPS, TAUC)
          TN=MATMUL(MATMUL(TRANSPOSE(FE),
     1             MATMUL(TAUC,TRANSPOSE(FTOINV))),TRANSPOSE(FPAUSN))
          DO ISYST=1,NTRSYS
            TRSTRE(ISYST)=SUM(TN*HVD0M(:,:,ISYST))
          ENDDO
          
          ISTRAN=TRSTRE > DGCR
        
          DO ISYST=1,NTRSYS
            IF(ISTRAN(ISYST))THEN
              RESVEC(ISYST)=(R1 + TRMU*GAMMAD(ISYST)/DTIME)**TREPS
     1                     -TRSTRE(ISYST)/DGCR
            ELSE
              RESVEC(ISYST)=GAMMAD(ISYST)
            ENDIF
          ENDDO
          
C Check for convergence
C======================
          RESNOR=NORM2(RESVEC(1:NTRSYS))
          FACTOR=NORM2(TRSTRE(1:NTRSYS)/DGCR)
          IF(FACTOR/=R0) RESNOR=RESNOR/FACTOR
          
          !WRITE(*,*)'nr marten: ', ITER, RESNOR
          !WRITE(*,*)GAMMAD
          !READ(*,*)
          
          
          IF(RESNOR<TOL)THEN
C Newton-Raphson iterations converged: Break loop and update stresses
C before exit
CC1            IF(ANY(TRSTRE < DGCR))THEN
CC1              WRITE(*,*)'PHI: N-R loop converged but not valid'
CC1              WRITE(*,*)TRSTRE < DGCR
CC1              WRITE(*,*)TRSTRE
CC1              READ(*,*)
CC1            ENDIF
            IF(ANY(GAMMAD < R0))THEN
              WRITE(*,*)'GAMMAD: N-R loop converged but not valid'
              WRITE(*,*)GAMMAD < R0
              WRITE(*,*)GAMMAD
              !READ(*,*)
              SUFAIL=.TRUE.
              GOTO 999
            ENDIF
            GOTO 200
          ENDIF
          
          
        
C Compute jacobian (numerical)
C=============================
C          EPS=1.0D-8 ! Perturbation factor
          EPS=1.0D-7 ! Perturbation factor
          DREMTX=R0
          DO JSYST=1,NTRSYS
            GAMMADP=GAMMAD
            GAMMADP(JSYST)=GAMMADP(JSYST) + EPS
C update FTRA, FE, TAU, T and TRSTRE (to calculate perturbed residual)
            GAMMAP=RSTAVA(20:43)+GAMMADP
            FTRAP=DELTA
            DO I=1,NTRSYS
              FTRAP=FTRAP + GAMMAP(I)*HVD0M(:,:,I)
            ENDDO
            FINELP=MATMUL(FTRAP,FPAUSN)
            CALL INVMT3(FINELP, FININVP, DETFINP)
            FEP=MATMUL(FTOTA,FININVP)
            CALL TAUMEP(FEP, NDIM, RPROPS, TAUP) ! Perturbed Kirchhoff stress
            TP=MATMUL(MATMUL(TRANSPOSE(FEP),
     1             MATMUL(TAUP,TRANSPOSE(FTOINV))),TRANSPOSE(FPAUSN))
            DO ISYST=1,NTRSYS
              TRSTREP(ISYST)=SUM(TP*HVD0M(:,:,ISYST))
            ENDDO
            
            
            ISTRANP=TRSTREP > DGCR
            
            DO ISYST=1,NTRSYS
              IF(ISTRANP(ISYST))THEN
                RESVECP(ISYST)=(R1 + TRMU*GAMMADP(ISYST)/DTIME)**TREPS
     1                        -TRSTREP(ISYST)/DGCR
              ELSE
                RESVECP(ISYST)=GAMMADP(ISYST)
              ENDIF
            ENDDO
            DREMTX(:,JSYST)=(RESVECP-RESVEC)/EPS
          ENDDO
          
          
C Solve linear algebraic system
C==============================
          CALL GAUSEL
     1(   DREMTX(1:NTRSYS,1:NTRSYS) ,RESVEC(1:NTRSYS) ,NTRSYS  ,SINGUL )
          IF(SINGUL)THEN
C Linear system is singular: switch state update failure flag on, break 
C Newton-Raphson loop and exit without updating stresses
            SUFAIL=.TRUE.
            WRITE(*,*)'GAUSEL failed! SUMTSC'!CALL ERRPRT('WE0007')
            GOTO 999
          ENDIF
C Apply Newton correction
          GAMMAD=GAMMAD-RESVEC
          GAMMA=RSTAVA(20:43)+GAMMAD
        
C End of Newton Raphson for transformation multiplier
        ENDDO nr
C IF stress update procedure failed to converge: Break loop and exit
        SUFAIL=.TRUE.
        WRITE(*,*)'ERROR: Stress update procedure failed to converge'
        GOTO 999
C Check whether or not transformation has completed      
  200   CONTINUE
          
        
C Reset transformation multiplier equal to ONE if transformation 
C has completed, then update transformation and elastic deformation 
C gradients along with work conjugate of transformation deformation 
C -----------------------------------------------------------------
        GAMMAT=SUM(GAMMA)
        
        
CC1 Strategy 1: proportionally reduce all gammas to satisfy gamma = 1
        IF(GAMMAT > R1)THEN
C          WRITE(*,*)'Warning: GAMMA > 1 in SUMTSC - readjusting...'
C Simple strategy for update: 
          GAMMA=GAMMA/GAMMAT
          GAMMAT=SUM(GAMMA)
C and reupdate variables
          FTRA=DELTA
          DO I=1,NTRSYS
            FTRA=FTRA + GAMMA(I)*HVD0M(:,:,I)
          ENDDO
          FINEL=MATMUL(FTRA,FPAUSN)
          CALL INVMT3(FINEL, FININV, DETFIN)
          FE=MATMUL(FTOTA,FININV)
          CALL TAUMEP(FE, NDIM, RPROPS, TAUC)
          TN=MATMUL(MATMUL(TRANSPOSE(FE),
     1             MATMUL(TAUC,TRANSPOSE(FTOINV))),TRANSPOSE(FPAUSN))
          DO ISYST=1,NTRSYS
            TRSTRE(ISYST)=SUM(TN*HVD0M(:,:,ISYST))
          ENDDO
        ENDIF
        
CC1 Strategy 2: sub-stepping
C Bisection of incremental deformation gradient
C When transformation is completed, use remaining incremental 
C deformation gradient in elastic part
C --> consequences for consistent tangent??
CC1 Strategy 3: increment cutting (use only for comparing results!)
c        TOLER=1.0D-3
c        IF(GAMMAT > R1+TOLER)THEN
c          SUFAIL=.TRUE.
c          GOTO 999
c        ENDIF
C
      ELSE
C7C        
C7C Check plastic admissibility 
C7C ===========================
C7C... and if PLASTIC then compute the viscoplastic contribution.         
C7C5        IF(CKPLAS.AND.IFPRTR.EQ.R0)THEN
C7C5        CALL SUVCPL
C7C5     1    ( DTIME      ,FINCR      ,IPROPS     ,LALGVA     ,NTYPE     ,
C7C5     2      RPROPS     ,RSTAVA     ,STRES      )
C7        DTEMP=RSTAVA(10:18)
C7        RSTAVA(10:18)=RSTAVA(42:50) ! copy FP to its place
C7        DTEMP2=RSTAVA(19)
C7        RSTAVA(19)=RSTAVA(41) ! copy HRVAR
C7        CALL SUVSC2
C7     1(   DTIME      ,FINCR      ,IPROPS     ,LALGVA     ,NTYPE      ,
C7     2    RPROPS     ,RSTAVA     ,STRES      )
C7        RSTAVA(42:50)=RSTAVA(10:18) ! new FP
C7        RSTAVA(41)=RSTAVA(19) ! new HRVAR
C7        RSTAVA(10:18)=DTEMP
C7        RSTAVA(19)=DTEMP2
C7        IFPLAS=LALGVA(1)
C7        SUFAIL=LALGVA(2)
C7C
C7C5        SUFAIL=LALGVA(2)
C7C5        IFPLAS=LALGVA(3)
C7        
C7C!!!!!!! updates in SUVCPL will be overwriten by updates below???
C7        
C7          IF(SUFAIL)GOTO 999
C7          IF(IFPLAS)THEN
C7C Check transformation admissibility
C7C Retrieve updated elastic and plastic deformation gradients
C7            FE=RESHAPE(RSTAVA(1:9),[3,3])
C7            FPAUSN=RESHAPE(RSTAVA(42:50),[3,3])
C7            FTOTA=MATMUL(FE,MATMUL(FTRAN,FPAUSN))
C7            CALL INVMT3(FTOTA, FTOINV, DETFTO)
C7C !!!! FTOINV should still remain the same (even when FPAUSN was updated)?
C7
C7C Update work conjugate of transformation deformation gradient and
C7C transformation function
C7            CALL TAUMEP(FE, NDIM, RPROPS, TAUC) ! Kirchhoff stress
C7C8
C7C8            TN=MATMUL(TRANSPOSE(FE),MATMUL(TAUC,TRANSPOSE(FTOINV)))
C7C8
C7            TN=MATMUL(MATMUL(TRANSPOSE(FE),
C7     1                MATMUL(TAUC,TRANSPOSE(FTOINV))),TRANSPOSE(FPAUSN))
C7C8
C7            
C7C If transformation system has not already been chosen, select the most
C7C favourable one (highest transformation function value)
C7            DO ISYST=1,NTRSYS
C7              PHIVAR(ISYST)=SUM(TN*HVD0M(:,:,ISYST))-DGCR
C7            ENDDO
C7            IF(IACTRS==0) IACTRS=MAXLOC(PHIVAR,1)
C7            PHI=PHIVAR(IACTRS)
C7C If transformation happened during plastic update, retry step with
C7C smaller increment
C7            IF(PHI/DGCR > TOL)THEN
C7              SUFAIL=.TRUE.
C7              WRITE(*,*) 'WARNING: Transformation started in plastic',
C7     1                   'update - increment cutting...'
C7            ENDIF
C7          ENDIF
C7          GOTO 999
C       
C Elastic step: Trial states are all actual
C =========================================
C!!!!!!!!
  666   CONTINUE
C!!!!!!!!
        FE=FETRL
        FTRA=FTRAN
        TAUC=TAUTRC
        TN=TNTRL
      ENDIF
C
C Update Cauchy stresses and other real state variables
C =====================================================
      CALL SYMTA(TAUC, STRES, NDIM, .TRUE.)
      STRES=STRES/DETFTO
C... update elastic deformation gradient components
      RSTAVA(1:9)=RESHAPE(FE,[9])
C... update transformation deformation gradient components
CC1      RSTAVA(10:18)=RESHAPE(FTRA,[9])
C... update transformation multiplier
CC1      RSTAVA(19)=GAMMA
      RSTAVA(20:43)=GAMMA
      RSTAVA(44)=GAMMAT
CC1
C... update work conjugate of transformation deformation gradient
CC1      RSTAVA(20:28)=RESHAPE(TN,[9])
C!!!!!!!! Not used anywhere, except for ORMEPC (for dissipation calculation)
C... update total deformation gradient components
CC1      RSTAVA(29:37)=RESHAPE(FTOTA,[9])
C!!! need to store it? can be calculated from the others
C... update differential volume (for volume fraction calculation)
CC1      RSTAVA(38)=DVOLU
C... update active variant of the transformation system
CC1      RSTAVA(39)=IACTRS
C 
      
C Calculate transformation dissipation
CC1      RSTAVA(52)=RSTAVA(52)+(GAMMA-GAMMAN)*SUM(TN*HVD0M(:,:,IACTRS))
C Calculate total dissipation
CC1      CALL ATSYM(STRES,SIGMA,NDIM,.TRUE.)
CC1      CALL INVMT3(FTOTA, FTOTINV, DETFTOT)
CC1      PIOLA=DETFTOT*MATMUL(SIGMA, TRANSPOSE(FTOTINV))!Piola-Kirchhoff stress: P = J sigma F^-T
CC1      RSTAVA(51)=RSTAVA(51)+SUM(PIOLA*(FTOTA-FTOTN))
      
      
CC1      WRITE(16,'(A,2E20.8)')'Dissipation:', RSTAVA(51:52)
      
      IF(IFTRAN)THEN
        !DO I=1,NTRSYS
        !  WRITE(*,*)I,SUM(TN(:,:)*HVD0M(:,:,I))/1.0d6,GAMMAD(I)
        !ENDDO
        !READ(*,*)
      ENDIF
      
C Incremental dissipation
c      DDISSIP=R0
c      IF(GAMMAT==R1)THEN
c        DDISSIP=R0
c      ELSE
c      DO I=1,NTRSYS
C          DDISSIP=DDISSIP+DTIME*GAMMAD(I)*SUM(TN(:,:)*HVD0M(:,:,I))
c        IF(ISTRAN(I))THEN
c          DDISSIP=DDISSIP+GAMMAD(I)*SUM(TN(:,:)*HVD0M(:,:,I))
c        ENDIF
c      ENDDO
c      ENDIF
c      DISSIP=DISSIP+DDISSIP
      
C Total dissipation
c      TDISSIP=R0
c      DO I=1,NTRSYS
c        TDISSIP=TDISSIP+GAMMA(I)*SUM(TN(:,:)*HVD0M(:,:,I))
c      ENDDO
c      WRITE(*,*)'Dissipation: ',DISSIP/1.0d6, TDISSIP/1.0d6
c      IF(GAMMAT>=R1)READ(*,*)
      
      
      
      
  999 CONTINUE
C Update some algorithmic variables before exit
C ---------------------------------------------
      LALGVA(1)=IFTRAN
      LALGVA(2)=SUFAIL
      LALGVA(3)=IFPLAS
C 
      
      
C8      RSTAVA(52)=R0
C8      IF(IFTRAN)RSTAVA(52)=R1
      
      
C	  
      RETURN
      END
CDOC END_SUBROUTINE SUMTSC
