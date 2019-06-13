CDOC BEGIN_SUBROUTINE CTDP
CDOC Consistent tangent matrix for the Drucker-Prager model
CDOC
CDOC The tangent matrix computed in this routine is consistent with
CDOC fully implicit elastic predictor/return mapping algorithm for the
CDOC Drucker-Prager elasto-plastic material model with piece-wise linear
CDOC isotropic hardening carried out in subroutine SUDP.
CDOC It contains only the plane strain and axisymmetric implementations.
CDOC It returns either the elastic tangent or the elasto-plastic
CDOC consistent tangent matrix depending on the input value of
CDOC the logical argument EPFLAG.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DGAM   >  Incremental plastic multipliers obtained
CDOC C                          in routine SUDP.
CDOC DOUBLE_PRECISION DMATX  <  Consistent tangent matrix.
CDOC LOGICAL          EPFLAG >  Elasto-plastic flag.
CDOC C                          If .FALSE., DMATX returns as the elastic
CDOC C                          matrix. If .TRUE., DMATX returns as the
CDOC C                          elasto-plastic tangent consistent with
CDOC C                          the return mapping algorithm implemented
CDOC C                          in routine SUDP.
CDOC INTEGER          IPROPS >  Array of integer material properties.
CDOC C                          This array is set in routines
CDOC C                          INDATA and RDDP.
CDOC C                          The number of points on the piece-wise
CDOC C                          linear hardening curve is the only
CDOC C                          element of this array used here.
CDOC LOGICAL          LALGVA >  Array of logical algorithmic flags.
CDOC C                          See list of arguments of SUDP.
CDOC INTEGER          NTYPE  >  Stress state type.
CDOC DOUBLE_PRECISION RPROPS >  Array of real material properties.
CDOC C                          Same as in the argument list of
CDOC C                          subroutine SUDP.
CDOC DOUBLE_PRECISION RSTAVA >  Array of real state variables other
CDOC C                          than the stress tensor components.
CDOC C                          Output of SUDP.
CDOC DOUBLE_PRECISION STRAT  >  Array of elastic trial (engineering)
CDOC C                          strain components last used as input of
CDOC C                          subroutine SUDP.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto and P.H.Saksono, June 1996: Initial coding
CHST
      SUBROUTINE CTDP3D
     1(   DGAM       ,DMATX      ,EPFLAG     ,IPROPS     ,LALGVA     ,
     2    NTYPE      ,RPROPS     ,RSTAVA     ,STRAT      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(IPHARD=7  ,MSTRE=6 ,NSTRE=6)
      LOGICAL APEX, EPFLAG, LALGVA(3)
      DIMENSION
     1    DMATX(MSTRE,MSTRE),IPROPS(*)           ,RPROPS(*)          ,
     2    RSTAVA(MSTRE+1)   ,STRAT(MSTRE)
      DIMENSION
     1    EETD(MSTRE)        ,FOID(MSTRE,MSTRE)  ,SOID(MSTRE)        ,
     2    UNIDEV(MSTRE)
      DATA
     1    FOID(1,1),FOID(1,2),FOID(1,3),FOID(1,4),FOID(1,5),FOID(1,6)/
     2    1.0D0    ,0.0D0    ,0.0D0    ,0.0D0    ,0.0D0    ,0.0D0 /
     3    FOID(2,1),FOID(2,2),FOID(2,3),FOID(2,4),FOID(2,5),FOID(2,6)/
     4    0.0D0    ,1.0D0    ,0.0D0    ,0.0D0    ,0.0D0    ,0.0D0 /
     5    FOID(3,1),FOID(3,2),FOID(3,3),FOID(3,4),FOID(3,5),FOID(3,6)/
     6    0.0D0    ,0.0D0    ,1.0D0    ,0.0D0    ,0.0D0    ,0.0D0 /
     7    FOID(4,1),FOID(4,2),FOID(4,3),FOID(4,4),FOID(4,5),FOID(4,6)/
     8    0.0D0    ,0.0D0    ,0.0D0    ,0.5D0    ,0.0D0    ,0.0D0/
     9    FOID(5,1),FOID(5,2),FOID(5,3),FOID(5,4),FOID(5,5),FOID(5,6)/
     O    0.0D0    ,0.0D0    ,0.0D0    ,0.0D0    ,0.5D0    ,0.0D0/
     1    FOID(6,1),FOID(6,2),FOID(6,3),FOID(6,4),FOID(6,5),FOID(6,6)/
     2    0.0D0    ,0.0D0    ,0.0D0    ,0.0D0    ,0.0D0    ,0.5D0/
      DATA
     1    SOID(1)  ,SOID(2)  ,SOID(3)  ,SOID(4)   ,SOID(5)  ,SOID(6) /
     2    1.0D0    ,1.0D0    ,1.0D0    ,0.0D0     ,0.0D0    ,0.0D0 /
      DATA
     1    R0   ,R1   ,RP5  ,R2   ,R3   /
     2    0.0D0,1.0D0,0.5D0,2.0D0,3.0D0/
C***********************************************************************
C COMPUTATION OF CONSISTENT TANGENT MODULUS FOR DRUCKER-PRAGER TYPE
C ELASTO-PLASTIC MATERIAL WITH ASSOCIATIVE/NON-ASSOCIATIVE FLOW RULE AND
C PIECE-WISE LINEAR ISOTROPIC HARDENING
C
C REFERENCE: Section 8.3.5
C***********************************************************************
      IF(NTYPE.NE.4)CALL ERRPRT('EI0039')
C Retrieve accumulated plastic strain, DGAMA and APEX algorithm flag
      EPBAR=RSTAVA(MSTRE+1)
      DGAMA=DGAM
      APEX=LALGVA(3)
C Set some material properties
      YOUNG=RPROPS(2)
      POISS=RPROPS(3)
      ETA=RPROPS(4)
      XI=RPROPS(5)
      ETABAR=RPROPS(6)
      NHARD=IPROPS(3)
C and some constants
      GMODU=YOUNG/(R2*(R1+POISS))
      BULK=YOUNG/(R3*(R1-R2*POISS))
      R2G=R2*GMODU
      R1D3=R1/R3
      ROOT2=SQRT(R2)
C
      IF(EPFLAG)THEN
C Compute elastoplastic consistent tangent
C ========================================
C Hardening slope
        HSLOPE=DPLFUN(EPBAR,NHARD,RPROPS(IPHARD))
        IF(APEX)THEN
C Elastoplastic tangent consistent with apex return
C -------------------------------------------------
          ALPHA=XI/ETABAR
          BETA=XI/ETA
          AFACT=BULK*(R1-BULK/(BULK+ALPHA*BETA*HSLOPE))
          DO 20 I=1,NSTRE
            DO 10 J=1,NSTRE
              DMATX(I,J)=AFACT*SOID(I)*SOID(J)
   10       CONTINUE       
   20     CONTINUE
        ELSE
C Elastoplastic tangent consistent with smooth cone wall return
C -------------------------------------------------------------
C Elastic trial deviatoric (physical) strain
          EEVD3=(STRAT(1)+STRAT(2)+STRAT(3))*R1D3
          EETD(1)=STRAT(1)-EEVD3
          EETD(2)=STRAT(2)-EEVD3
          EETD(4)=STRAT(4)*RP5
          EETD(5)=STRAT(5)*RP5
          EETD(6)=STRAT(6)*RP5
          EETD(3)=STRAT(3)-EEVD3
          ETDNOR=SQRT(EETD(1)*EETD(1)+EETD(2)*EETD(2)+
     1             R2*EETD(4)*EETD(4)+R2*EETD(5)*EETD(5)+
     2             R2*EETD(6)*EETD(6)+EETD(3)*EETD(3))
C Unit deviatoric flow vector
          IF(ETDNOR.NE.R0)THEN
            EDNINV=R1/ETDNOR
          ELSE
            EDNINV=R0
          ENDIF
          DO 30 I=1,NSTRE
            UNIDEV(I)=EETD(I)*EDNINV
   30     CONTINUE
C Assemble tangent
          AUX=R1/(GMODU+BULK*ETA*ETABAR+XI*XI*HSLOPE)
          AFACT=R2G*(R1-DGAMA/(ROOT2*ETDNOR))
          AFACD3=AFACT*R1D3
          BFACT=R2G*(DGAMA/(ROOT2*ETDNOR)-GMODU*AUX)
          CFACT=-ROOT2*GMODU*BULK*AUX
          DFACT=BULK*(R1-BULK*ETA*ETABAR*AUX)
          DO 50 I=1,NSTRE
            DO 40 J=1,NSTRE
              DMATX(I,J)=AFACT*FOID(I,J)+BFACT*UNIDEV(I)*UNIDEV(J)+
     1                   CFACT*(ETA*UNIDEV(I)*SOID(J)+
     2                           ETABAR*SOID(I)*UNIDEV(J))+
     3                   (DFACT-AFACD3)*SOID(I)*SOID(J)
   40       CONTINUE       
   50     CONTINUE
        ENDIF
      ELSE
C Compute elasticity matrix
C =========================
        FACTOR=BULK-R2G*R1D3
        DO 70 I=1,NSTRE
          DO 60 J=I,NSTRE
            DMATX(I,J)=R2G*FOID(I,J)+FACTOR*SOID(I)*SOID(J)
   60     CONTINUE       
   70   CONTINUE
        DO 90 J=1,NSTRE-1
          DO 80 I=J+1,NSTRE
            DMATX(I,J)=DMATX(J,I)
   80     CONTINUE
   90   CONTINUE
      ENDIF
      RETURN
      END
CDOC END_SUBROUTINE CTDP
