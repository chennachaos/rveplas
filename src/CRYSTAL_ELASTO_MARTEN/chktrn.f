CDOC BEGIN_SUBROUTINE CHKTRN
CDOC Check for transformation admissibility 
CDOC
CDOC BEGIN_PARAMETERS    
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
CDOC END_PARAMETERS
CHST
CHST M. Fauzan Adziman, April 2013: Initial coding
CHST
      SUBROUTINE CHKTRN
     1(   IPROPS     ,IFTRAN     ,NTYPE      ,RPROPS     ,RSTAVA     )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)	
      PARAMETER
     1(   NIPROP=18  ,NRSTAV=52    )
C Arguments
      LOGICAL IFTRAN
      DIMENSION
     1    IPROPS(NIPROP)     ,
     2    RPROPS(*)          ,RSTAVA(NRSTAV)
C Local arrays and variables for martensitic transformation
      DIMENSION
     1    FTOTA(3,3) ,FPAUSN(3,3),
     2    FETRL(3,3) ,
     3    FTRAN(3,3) ,FTOINV(3,3)
      DIMENSION
     2    HVD0M(3,3,IPROPS(5))
      DIMENSION
     1    PAUXTR(3,IPROPS(5)),
     2    PHIVAR(IPROPS(5))  ,PHIAUX(IPROPS(5))  ,
     3    PHIASC(IPROPS(5))  ,TNTRL(3,3) ,
     4    TAUTRC(3,3)
C Predefined data used in this routine
      DOUBLE PRECISION, PARAMETER :: TOL=1.0D-8
C
C***********************************************************************
C CHECK TRANSFORMATION ADMISSIBILITY
C***********************************************************************
C
C Check analysis type is valid for current model
      IF(NTYPE==2)THEN
        NDIM=2
      ELSEIF(NTYPE==4)THEN
        NDIM=3
      ELSE
        STOP 'Invalid NTYPE - CHKTRN' !CALL ERRPRT('EI0036')
      ENDIF
C Initialise some algorithmic and internal variables
      IFTRAN=.FALSE.
C... input parameters for the transformation
C... the required amount of work to start a transformation	  
      DGCR=RPROPS(5)
C... retrieve activated variant in the transformation system
      IACTRS=INT(RSTAVA(39)) 
C... flag to identify whether further variant selection is required or 
C not	
      IFVARS=INT(RSTAVA(40))	
C Recover systems information from RPROPS
      NTRSYS=IPROPS(5)
      HVD0M=RESHAPE(RPROPS(IPROPS(10):IPROPS(11)),[3,3,NTRSYS]) ! System constant matrix: d (x) m
C	
C Compute trial states of deformation gradients
C ---------------------------------------------
C... retrieve updated elastic transformation deformation gradient from
C    the plastic update state
      FETRL=RESHAPE(RSTAVA(1:9),[3,3])
C... retrieve previous state of austenite plastic deformation gradient
      FPAUSN=RESHAPE(RSTAVA(42:50),[3,3])
C... retrieve previous state of transformation deformation gradient
      FTRAN=RESHAPE(RSTAVA(10:18),[3,3])
C... compute total deformation gradient
      FTOTA=MATMUL(FETRL,MATMUL(FTRAN,FPAUSN))
C... compute inverse of total deformation gradient
      CALL INVMT3(FTOTA, FTOINV, DETFTO)
C Compute work conjugate trial of transformation deformation gradient
C -------------------------------------------------------------------
      CALL TAUMEP(FETRL, NDIM, RPROPS, TAUTRC)
C ... assemble work conjugate trial of transformation deformation 
C     gradient
      TNTRL=MATMUL(TRANSPOSE(FETRL),MATMUL(TAUTRC,TRANSPOSE(FTOINV)))
C      
C Check martensitic transformation admissibility	
C ==============================================       
C Compute transformation function value at trial state  
C ----------------------------------------------------
      DO ISYST=1,NTRSYS
        PHIVAR(ISYST)=SUM(TNTRL*HVD0M(:,:,ISYST))-DGCR
      ENDDO
      PHIAUX=PHIVAR/DGCR-TOL
      PHIASC=PHIAUX
C... select the most favourable variant for transformation
      IF(IFVARS.EQ.0)THEN
C... sort variants (in an ascending order) according to their potential   
C    to initiate a transformation
          DO I=1,NTRSYS-1
            DO J=I+1,NTRSYS
              IF(PHIAUX(I).GT.PHIAUX(J))THEN
                 AUXLIST=PHIASC(I)
                 PHIASC(I)=PHIASC(J)
                 PHIASC(J)=AUXLIST		  
              ENDIF
            END DO
          END DO
C... the bottom data in the sorted array is the most positive value, 
C    hence the most favoured variant to transform	    
          PHIDEL=PHIASC(NTRSYS)      
        DO ISYST=1,NTRSYS
          IF(PHIDEL.EQ.PHIAUX(ISYST))THEN
            IACTRS=ISYST
            PHI=PHIVAR(IACTRS)
          ENDIF
        END DO
C... once a transformation has been initiated at a gauss point, the
C    selected variant at the gauss point to be preserved throughout the
C    transformation process
      ELSE
        PHI=PHIVAR(IACTRS)
      ENDIF
C... check transformation consistency
      IF(PHI/DGCR.GT.TOL)IFTRAN=.TRUE.
C 
      RETURN
      END
CDOC END_SUBROUTINE CHKTRN