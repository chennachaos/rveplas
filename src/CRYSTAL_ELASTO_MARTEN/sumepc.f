CDOC BEGIN_SUBROUTINE SUMEPC
CDOC State update for thermodynamically consistent martensitic elasto-
CDOC plastic-like single crystal model. HYPLAS
CDOC
CDOC This routine uses an implicit elastic predictor/return mapping
CDOC algorithm as the state update procedure. This model is compatible 
CDOC only with the plane strain assumption.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DVOLU  >  Differential volume.      
CDOC DOUBLE_PRECISION FINCR  >  Incremental deformation gradient.
CDOC INTEGER          IPROPS >  Array of integer material properties.
CDOC C                          This array is set in routines RVDATA 
CDOC C                          and RDMEPC.
CDOC LOGICAL          LALGVA <  Array of logical algorithmic flags.
CDOC C                          For the present material model, this
CDOC C                          array contains the plastic-like yielding
CDOC C                          flag, IFTRAN, the state update failure
CDOC C                          flag, SUFAIL.
CDOC C                          The plastic-like yielding flag is set to
CDOC C                          .TRUE. if transformation has occurred
CDOC C                          and to .FALSE. if the step is elastic.
CDOC C                          The algorithm failure flag is set to
CDOC C                          .FALSE. if the state update algorithm
CDOC C                          has been successful and to .TRUE. if the
CDOC C                          return mapping algorithm has failed for
CDOC C                          any reasons.
CDOC INTEGER          NTYPE  >  Stress state type flag. The present
CDOC C                          implementation is compatible only with
CDOC C                          plane strain (NTYPE=2).
CDOC DOUBLE_PRECISION RPROPS >  Array of real material properties.
CDOC C                          This array contains the density
CDOC C                          (not used in this routine), shear and
CDOC C                          bulk modulus (neo-Hookean constants),
CDOC C                          critical free energy for activation,
CDOC C                          shear and dilatation factor. This array
CDOC C                          is set in routine RDPDSC. 
CDOC DOUBLE_PRECISION RSTAVA <> Array of real state variables other
CDOC C                          than the stress tensor components.
CDOC C                          Previous converged values on entry,
CDOC C                          updated values on exit.
CDOC DOUBLE_PRECISION STRES  <  Array of updated stress tensor
CDOC C                          components.
CDOC END_PARAMETERS
CHST
CHST M. Fauzan Adziman, April 2013: Initial coding
CHST
      SUBROUTINE SUMEPC
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
      
      
C Local parameters
      DOUBLE PRECISION, PARAMETER :: R0=0.0D0     , R1=1.0D0        ,
     1                               TOL=1.0D-8
      INTEGER, PARAMETER          :: MXITER=50
      DOUBLE PRECISION, PARAMETER, DIMENSION(3,3) :: ! delta_ij
     1     DELTA=reshape((/R1,R0,R0,R0,R1,R0,R0,R0,R1/), (/3,3/) )
C***********************************************************************
C STRESS UPDATE PROCEDURE FOR THE ANISOTROPIC TRANSFORMATION SINGLE
C CRYSTAL VISCOPLASTIC MARTENSITIC TRANSFORMATION MODEL 
C
C Simplified algorithm:
C ---------------------
C    1. Retrieve material and system properties
C    2. Compute elastic trial state
C    3. Compute transformation trial state
C    4. Compute total deformation gradient at trial state
C    5. Compute work conjugate of transformation deformation gradient
C    6. Compute yield function value at trial state 
C    7. Check transformation admissibility
C
C A.If IFTRAN=.FALSE. then
C   A1. Elastic step: trial state is the actual state
C   A2. Update cauchy stress components
C   A3. Update elastic and transformation deformation gradient
C   A4. Update some algorithmic variables before exit
C
C B.If IFTRAN=.TRUE. then
C   B1. Transformation step: apply return mapping
C   B2. Start Newton-Raphson iteration for the transformation multiplier
C   B3.    Compute jacobian matrix of the return-mapping 
C   B4.    Apply Newton-Raphson correction 
C   B5.    Update transformation deformation gradient, update FEISO, the
C          work conjugate of transformation deformation gradient and the
C          current state of yield function
C   B6.    Check for convergence
C   B7. End of a Newton Raphson step
C   B8. Update gamma (state of transformation)
C       If(GAMMA.GE.R1) then GAMMA=R1 -> update FEISO endif
C       update for transformation deformation gradient and the work 
C       conjugate	
C   B9. Update Cauchy stress components
C   B10.Update elastic and transformation deformation gradient
C   B11.Update some algorithmic variables before exit
C
C***********************************************************************
      
C!!!!!!! storing things in G-matrix ordering allows double contractions
C        between two 4th-order tensors to be expressed as matrix products
C        -> explore in martensitic model SU and CT !!
      
      
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
      GAMMAN=RSTAVA(19)
      FTOTN=RESHAPE(RSTAVA(29:37),[3,3])
C
      
C7
      
C First, perform plastic return-mapping (only if transformation has not already started)
C      IF(.NOT.IFTRN2)THEN
C      IF(.NOT.IFTRN2.AND.GAMMA==R0)THEN
      GAMMA=RSTAVA(19) ! Previous state of transformation multiplier
      IF(GAMMA==R0)THEN
C Check plastic admissibility 
C ===========================
C... and if PLASTIC then compute the viscoplastic contribution.         
C5        IF(CKPLAS.AND.IFPRTR.EQ.R0)THEN
C5        CALL SUVCPL
C5     1    ( DTIME      ,FINCR      ,IPROPS     ,LALGVA     ,NTYPE     ,
C5     2      RPROPS     ,RSTAVA     ,STRES      )
        IFPLUP=.TRUE.
        DTEMP=RSTAVA(10:18)
        RSTAVA(10:18)=RSTAVA(42:50) ! copy FP to its place
        DTEMP2=RSTAVA(19)
        RSTAVA(19)=RSTAVA(41) ! copy HRVAR
        CALL SUVSC2
     1(   DTIME      ,FINCR      ,IPROPS     ,LALGVA     ,NTYPE      ,
     2    RPROPS     ,RSTAVA     ,STRES      )
        RSTAVA(42:50)=RSTAVA(10:18) ! new FP
        RSTAVA(41)=RSTAVA(19) ! new HRVAR
        RSTAVA(10:18)=DTEMP
        RSTAVA(19)=DTEMP2
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
      GAMMA=RSTAVA(19) ! Previous state of transformation multiplier
      IACTRS=INT(RSTAVA(39)) ! Activated transformation system
C7      IFVARS=INT(RSTAVA(40)) ! Identify whether transformation system selection is required
C7      IFPRTR=INT(RSTAVA(51)) ! State of transformation in the previous converged solution
C
C Trial deformation gradients:
C ----------------------------
C Retrieve previous elastic, austenite plastic and transformation
C deformation gradients
      FEN=RESHAPE(RSTAVA(1:9),[3,3])
      FPAUSN=RESHAPE(RSTAVA(42:50),[3,3])
      FTRAN=RESHAPE(RSTAVA(10:18),[3,3])
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
      ELSE
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
      IF(GAMMA==R1)GOTO 666
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
      IF(IACTRS==0)THEN !variant selection still needed
        DO ISYST=1,NTRSYS
          PHIVAR(ISYST)=SUM(TNTRL*HVD0M(:,:,ISYST))-DGCR
        ENDDO
C5        IACTRS=MAXLOC(PHIVAR,1)
C5        PHI=PHIVAR(IACTRS)
        PHI=MAXVAL(PHIVAR,1)
      ELSE
        PHI=SUM(TNTRL*HVD0M(:,:,IACTRS))-DGCR
      ENDIF
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
      IF(PHI/DGCR.GT.TOL)THEN
        IFTRAN=.TRUE.
C If transforming variant not yet defined, select it
        IF(IACTRS==0)THEN
          IACTRS=MAXLOC(PHIVAR,1)
        ENDIF
      ENDIF
      
      
C!!! add check if transformation happens directly from elasticity, if
C    plasticity would happen?? (and activate increment cutting...)
C    -> compare results using both approaches
                   
C Transformation step: Apply return-mapping
C =========================================        
C Setting initial states at the beginning of the return-mapping loop 
C ------------------------------------------------------------------
      IF(IFTRAN)THEN
        FE=FETRL
        FTRA=FTRAN
        TAUC=TAUTRC
C
C Start Newton-Raphson iterations to find solution of the return-
C mapping equation as a function of transformation multiplier 
        nr: DO ITER=1,MXITER
C Compute Jacobian matrix of the return-mapping (DPHI)
C ---------------------------------------------
C ... compute inverse of transformation trial deformation gradient
          CALL INVMT3(FTRA, FTRINV, DETFPN)
C DFEDDG: derivative of the elastic deformation gradient with respect to
C the increment of gamma (transformation multiplier)
C (F.PhD eq. 5.71 - or, equivalently, transpose of tensor M on F.PhD 
C  eq. 5.42)
C ---------------------------------------------------------------
          DFEDDG=-MATMUL(MATMUL(FE,HVD0M(:,:,IACTRS)),FTRINV)
C DTDFE: derivative of the Kirchoff stress with respect to the elastic
C deformation gradient (only part of jacobian that elastic constitutive 
C model - in this case, the regularised Neo-Hookean model)
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
C8          DPHIB=MATMUL(
C8     1        MATMUL(TRANSPOSE(DFEDDG),TAUC)+MATMUL(TRANSPOSE(FE),DPHIA)
C8     2        ,TRANSPOSE(FTOINV))
C8
          DPHIB=MATMUL(MATMUL(
     1        MATMUL(TRANSPOSE(DFEDDG),TAUC)+MATMUL(TRANSPOSE(FE),DPHIA)
     2        ,TRANSPOSE(FTOINV)),TRANSPOSE(FPAUSN))
C8
          
C DPHI: double contraction of DPHIB and HVD0M (of active transformation
C system)
          DPHI=SUM(DPHIB*HVD0M(:,:,IACTRS))
C Apply Newton-Raphson correction to the transformation multiplier
C ----------------------------------------------------------------
          DGAMMA=DGAMMA-PHI/DPHI
          GAMMA=RSTAVA(19)+DGAMMA
C Update return-mapping function after the correction
C ---------------------------------------------------
C Update transformation deformation gradient
          FTRA=DELTA+GAMMA*HVD0M(:,:,IACTRS)
C... update inelastic deformation gradient
          FINEL=MATMUL(FTRA,FPAUSN)
C... compute inverse of inelastic deformation gradient
          CALL INVMT3(FINEL, FININV, DETFNE)
C... update elastic deformation gradient
          FE=MATMUL(FTOTA,FININV)
C
C Update work conjugate of the transformation deformation gradient (TN):
C ----------------------------------------------------------------------
          CALL TAUMEP(FE, NDIM, RPROPS, TAUC) ! Kirchhoff stress
          
          
C8
C8          TN=MATMUL(TRANSPOSE(FE),MATMUL(TAUC,TRANSPOSE(FTOINV)))
C8
          TN=MATMUL(MATMUL(TRANSPOSE(FE),
     1              MATMUL(TAUC,TRANSPOSE(FTOINV))),TRANSPOSE(FPAUSN))
C8
      
      
C Compute new yield function value at the current iteration
C ---------------------------------------------------------
          PHI=SUM(TN*HVD0M(:,:,IACTRS))-DGCR
C Check for convergence
C ---------------------
          RESNOR=ABS(PHI)/DGCR
          IF(RESNOR<=TOL)THEN
C... N-R loop converged: check validity of the current solution
            IF(PHI/DGCR-TOL>R0)THEN
              WRITE(*,*)'PHI: N-R loop converged but not valid'
            ENDIF
            IF(DGAMMA<R0)THEN
             WRITE(*,*)'DGAMMA: N-R loop converged but not valid',DGAMMA
            ENDIF
C... stress updated converged: Break loop    
            GOTO 200
          ENDIF 
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
        IF(GAMMA > R1)THEN
          WRITE(*,*)'Warning: GAMMA > 1 in SUMEPC - readjusting...'
C??? activate increment cutting instead?
C!!!! 
          GAMMA=R1       
C Re-update transformation deformation gradient
C --------------------------------------
          FTRA=DELTA+GAMMA*HVD0M(:,:,IACTRS)
C... update inelastic deformation gradient
          FINEL=MATMUL(FTRA,FPAUSN)
C... compute inverse of inelastic deformation gradient
          CALL INVMT3(FINEL, FININV, DETFNE)
C          
C Re-update elastic deformation gradient
C --------------------------------------
          FE=MATMUL(FTOTA,FININV)
C Update work conjugate of transformation deformation gradient
C ------------------------------------------------------------
          CALL TAUMEP(FE, NDIM, RPROPS, TAUC) ! Update Kirchhoff stress
C8
C8        TN=MATMUL(TRANSPOSE(FE),MATMUL(TAUC,TRANSPOSE(FTOINV)))
C8
          TN=MATMUL(MATMUL(TRANSPOSE(FE),
     1              MATMUL(TAUC,TRANSPOSE(FTOINV))),TRANSPOSE(FPAUSN))
C8
        ENDIF     
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
      RSTAVA(10:18)=RESHAPE(FTRA,[9])
C... update transformation multiplier
      RSTAVA(19)=GAMMA
C... update work conjugate of transformation deformation gradient
      RSTAVA(20:28)=RESHAPE(TN,[9])
C!!!!!!!! Not used anywhere, except for ORMEPC (for dissipation calculation)
C... update total deformation gradient components
      RSTAVA(29:37)=RESHAPE(FTOTA,[9])
C!!! need to store it? can be calculated from the others
C... update differential volume (for volume fraction calculation)
      RSTAVA(38)=DVOLU
C... update active variant of the transformation system
      RSTAVA(39)=IACTRS
C 
      
C Calculate transformation dissipation
      RSTAVA(52)=RSTAVA(52)+(GAMMA-GAMMAN)*SUM(TN*HVD0M(:,:,IACTRS))
C Calculate total dissipation
      CALL ATSYM(STRES,SIGMA,NDIM,.TRUE.)
      CALL INVMT3(FTOTA, FTOTINV, DETFTOT)
      PIOLA=DETFTOT*MATMUL(SIGMA, TRANSPOSE(FTOTINV))!Piola-Kirchhoff stress: P = J sigma F^-T
      RSTAVA(51)=RSTAVA(51)+SUM(PIOLA*(FTOTA-FTOTN))
      
      
      WRITE(16,'(A,2E20.8)')'Dissipation:', RSTAVA(51:52)
      
      
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
CDOC END_SUBROUTINE SUMEPC
