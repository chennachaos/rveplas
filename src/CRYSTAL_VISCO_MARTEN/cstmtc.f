CDOC BEGIN_SUBROUTINE CSTMTC
      SUBROUTINE CSTMTC 
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
      
      
      
      
      DOUBLE PRECISION GAMMA(IPROPS(5)), FINCRP(3,3), FINCINV(3,3), 
     1       FTOTAP(3,3), STRESP(NSTRE), SIGMAP(3,3),
     2       B4TH(3,3,3,3)
      DOUBLE PRECISION RSTAVAP(NRSTAV)
      LOGICAL LALGVAP(NLALGV)
      
      
      
      
      
C Local parameters
      DOUBLE PRECISION, PARAMETER :: R0=0.0D0        , RP5=0.5D0       ,
     1                               R1=1.0D0        , R2=2.0D0        ,
     2                               R1D3=1.0D0/3.0D0, R2D3=2.0D0/3.0D0,
     3                               EPS=1.0D-7 ! perturbation factor
      DOUBLE PRECISION, PARAMETER, DIMENSION(3,3) :: ! delta_ij
     1     DELTA=reshape((/R1,R0,R0,R0,R1,R0,R0,R0,R1/), (/3,3/) )
C***********************************************************************
      
C      RETURN
      
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
        STOP 'invalid NTYPE - CSTMTC' !CALL ERRPRT('EI0082')
      ENDIF
      
C Numerical tangent evaluation
C=============================
C Matrix form of current stress
      CALL ATSYM( STRES, SIGMA, NDIM, .TRUE.)
C F total
      FE=RESHAPE(RSTAVA(1:9),[3,3])      ! Elastic
      FPAUS=RESHAPE(RSTAVA(10:18),[3,3]) ! Austenite plastic
      
      NTRSYS=IPROPS(5)
      GAMMA=RSTAVA(20:43)
      HVD0M=RESHAPE(RPROPS(IPROPS(10):IPROPS(11)),[3,3,NTRSYS])
      FTRA=DELTA
      DO I=1,NTRSYS
        FTRA=FTRA + GAMMA(I)*HVD0M(:,:,I)
      ENDDO
      
      FTOTA=MATMUL(MATMUL(FE,FTRA),FPAUS)
      DETF=DETM23(3, FTOTA, NDIM, .FALSE.)
      
      DETFINCR=DETM23(3, FINCR, NDIM, .FALSE.)
      
c      WRITE(*,*)'GAMMA total: ', SUM(GAMMA)
c      WRITE(*,'(*(G20.6))')GAMMA
c      DO I=1,3
c        WRITE(*,'(*(G20.6))')FTRA(I,:)
c      ENDDO
c      WRITE(*,*)DETM23(3, FTRA, NDIM, .FALSE.)
c      
c      WRITE(*,*)'GAMMA slip total: ', RSTAVA(19)
c      DO I=1,3
c        WRITE(*,'(*(G20.6))')FPAUS(I,:)
c      ENDDO
c      WRITE(*,*)DETM23(3, FPAUS, NDIM, .FALSE.)
      
      
      
      
      
      DO K=1,3
        DO L=1,3
C1          FINCRP=DELTA
          FINCRP=FINCR
          FINCRP(K,L)=FINCRP(K,L)+EPS
          
          DETFINCRP=DETM23(3, FINCRP, NDIM, .FALSE.)
          DETFP=DETF*DETFINCRP/DETFINCR
          
          
C1          RSTAVAP=RSTAVA(1:NRSTAV)
          RSTAVAP=RSTAVN(1:NRSTAV)
          LALGVAP=LALGVA(1:NLALGV)
          DGAMMAP=DGAMMA
          
          STRESP=R0
          
          CALL SUMTSC
     1(   DGAMMAP    ,FINCRP     ,IPROPS     ,LALGVAP    ,NTYPE      ,
     2    RPROPS     ,RSTAVAP    ,STRESP     ,DVOLU      ,IELEM      ,
     3    IGAUSP     ,DTIME      )
      
      
          CALL ATSYM( STRESP, SIGMAP, NDIM, .TRUE.)
          
          B4TH(:,:,K,L)=(DETFP*SIGMAP-DETF*SIGMA)/EPS
      
        ENDDO
      ENDDO
      
      
      DO I=1,3
        DO J=1,3
          DO K=1,3
            DO L=1,3
C1              A4TH(I,J,K,L)=SUM(B4TH(I,J,K,:)*FTOTA(L,:))/DETF
              A4TH(I,J,K,L)=SUM(B4TH(I,J,K,:)*FINCR(L,:))/DETF
     1                      -SIGMA(I,L)*DELTA(J,K)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      
      CALL ARRGO3( A4TH, AMATX  )
      ! generalise to include 2D too!
      ! compare to numerical tangent for plasticity if(ifplas)
      
c      WRITE(*,*)'NUMERICAL=========================='
c      DO I=1,9
c        WRITE(*,'(*(G20.6))')AMATX(I,:)
c      ENDDO
c      WRITE(*,*)'ANALYTICAL========================='
      
      
c      AMATX=R0
c      CALL CSTVS2
c     1(   AMATX ,  DTIME ,  EPFLAG,  FINCR,  IPROPS,  NTYPE,  RPROPS,
c     2    RSTAVA,  RSTAVN,  STRES   )
c      DO I=1,9
c        WRITE(*,'(*(G20.6))')AMATX(I,:)
c      ENDDO
c      READ(*,*)
      
      RETURN
      END
CDOC END_SUBROUTINE CSTMTC
