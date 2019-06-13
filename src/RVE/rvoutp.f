CDOC BEGIN_SUBROUTINE RVOUTP
CDOC Compute and print homogenised stresses
CDOC
CDOC This routine computes homogenised stress using boundary nodal
CDOC reaction-based formula      
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION FBD    >  Boundary nodal reaction vector.      
CDOC INTEGER          NOUTP  >  Load output flag that determines what
CDOC C                          results are to be printed out.
CDOC DOUBLE_PRECISION TFACT  >  Current total load factor.   
CDOC END_PARAMETERS      
CDOC    
CHST M.F. Adziman, July 2013: Initial coding. 
CHST                          Coded in RVEPLAS as in MICROPLAST-Matlab  
CHST                          written by EA de Souza Neto.
CHST D. de Bortoli, February 2015:
CHST                          Added conversion from first 
CHST                          Piola-Kirchhoff stress to Cauchy stress in
CHST                          large strains
CHST      
      SUBROUTINE RVOUTP ( IINCS, FBD, TFACT )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C Hyplas database
      INCLUDE '../MAXDIM.INC'
      INCLUDE '../MATERIAL.INC'
      INCLUDE '../ELEMENTS.INC'
      INCLUDE '../GLBDBASE.INC'
      INCLUDE '../RVE.INC' 
C Arguments      
      DIMENSION
     1   FBD(MTOTV)
C Local arrays and variables
      DIMENSION
     1   FBOUND(2,MTOTV) ,YBOUND(MPOIN,2),NODESB(MPOIN) ,  
     2   SIGAUX(3,3)     ,SIGAU2(3,3)    ,SIGMAB(3,3)   ,
     3   STRSH(MSTRE)    ,PSTRS(3)       ,F(3,3)        ,
     4   HOMO(3,3)       ,FFF(6)         ,STRAIND(6)
      PARAMETER( MGDIM=9 ,MBDIM=6 )
C
      DOUBLE PRECISION 
     1        ELCOD(NDIME,MNODE)  ,SHAPE(MNODE)   ,DERIV(NDIME,MNODE)  ,
     2        EISCRD(NDIME)       ,GPCOD(NDIME)   ,CARTD(NDIME,MNODE)  ,
     3        ELCOD1(NDIME,MNODE) ,FMICRO(3,3),    FAVER(3,3),
     4        TELDIS(MDOFN,MNODE) ,PMICRO(3,3),    SMICRO(3,3),
     5        PAVER(3,3)          ,SAVER(3,3) ,    SIGAU3(3,3),
     6        FMICIN(3,3)         ,GMATX(MGDIM,MEVAB)         ,
     7        STRSM(1,4,788,4)   ,STRSF(6,1,1,1)              ,
     8        STRSMASH(3,3)      
      DOUBLE PRECISION, PARAMETER :: R0=0.0D0, R1=1.0D0, R2=2.0D0,
     1                               R3=3.0D0
C***********************************************************************
C COMPUTATION AND OUTPUTS OF HOMOGENISED STRAINS AND STRESSES
C***********************************************************************
 2000 FORMAT(' Homogenised strains and stresses for increment ',I0,':'/,
     1        ' ---------------------------------------------------')
 2010 FORMAT(' e-xx  = ',G12.5,6X,' e-yy  = ',G12.5,6X,
     1       ' e-xy  = ',G12.5) 
 2012 FORMAT(' e-xx  = ',G12.5,6X,' e-yy  = ',G12.5,6X,
     1       ' e-zz  = ',G12.5,6X/' e-xy  = ',G12.5,6X,
     2       ' e-yz  = ',G12.5,6X,' e-xz  = ',G12.5)
c 2015 FORMAT(' F-xx  = ',G12.5,6X,' F-xy  = ',G12.5,6X/
c     1       ' F-yx  = ',G12.5,6X,' F-yy  = ',G12.5)
 2015 FORMAT( G12.5,6X)
 2017 FORMAT(' F-xx  = ',G12.5,6X,' F-xy  = ',G12.5,6X,
     1       ' F-xz  = ',G12.5,6X/' F-yx  = ',G12.5,6X,
     2       ' F-yy  = ',G12.5,6X,' F-yz  = ',G12.5,6X/
     3       ' F-zx  = ',G12.5,6X,' F-zy  = ',G12.5,6X,
     4       ' F-zz  = ',G12.5)
 2020 FORMAT(' S-xx  = ',G12.5,6X,' S-yy  = ',G12.5,6X,
     1       ' S-xy  = ',G12.5)
 4000 FORMAT(' S-xx  = ',G12.5,6X,' S-yy =' ,G12.5,6X,
     1        ' S-zz=',G12.5,6X,'S-xy=',G12.5,6X)
 4003 FORMAT(' S-xx  = ',G12.5,6X,' S-yy =' ,G12.5,6X,
     1       ' S-zz=',G12.5,6X,'S-xy=',G12.5,6X,
     2       ' S-yz=',G12.5,6X,'S-xz=',G12.5,6X)
 2022 FORMAT(' S-xx  = ',G12.5,6X,' S-yy  = ',G12.5,6X,
     1       ' S-zz  = ',G12.5,6X/' S-xy  = ',G12.5,6X,
     2       ' S-yz  = ',G12.5,6X,' S-xz  = ',G12.5)
 2025 FORMAT(' S-max = ',G12.5,6X,' S-min = ',G12.5,6X,
     1       ' Angle = ',G12.5)
 2027 FORMAT(' S-max = ',G12.5,6X,' S-mid = ',G12.5,6X,
     1       ' S-min = ',G12.5)
 2030 FORMAT(' S-eff = ',G12.5/    
     1       ' ======================================================='
     2       '===================')
 4001 FORMAT('',G12.5,6X'',G12.5)
 4009 FORMAT('',G12.5,6X'',G12.5,6X'',G12.5)
 4004 FORMAT('',G12.5,6X'',G12.5,6x'',G12.5,6x'',G12.5,6x'',G12.5)
C
C Computation of homogenised stress
C***********************************************************************
 1020 FORMAT(//' Results for increment ',I4,' Time = ',G12.6,1X,
     1         ' Load factor = ',G12.6/
     2 1X,74('=')//' Converged solution at iteration number =',2I5/)
 1030 FORMAT(//' Displacement of structure from initial configuration'/,
     1         ' ====================================================')
 1040 FORMAT(/' Node       R-Disp         Z-Disp')
 1045 FORMAT(/' Node       X-Disp         R-Disp')
 1050 FORMAT(/' Node       X-Disp         Y-Disp')
 1055 FORMAT(/' Node       X-Disp         Y-Disp         Z-Disp')
 1060 FORMAT(I5,2X,6G15.6)
 1070 FORMAT(//' Reactions'/1X,9('='))
 1080 FORMAT(/' Node      R-Force          Z-Force')
 1081 FORMAT(/' Node      X-Force          R-Force')
 1085 FORMAT(/' Node      R-Force          Z-Force          X-Local',
     1'          Y-Local')
 1086 FORMAT(/' Node      X-Force          R-Force          X-Local',
     1'          Y-Local')
 1090 FORMAT(/' Node      X-Force          Y-Force')
 1091 FORMAT(/' Node      X-Force          Y-Force          Z-Force')
 1095 FORMAT(/' Node      X-Force          Y-Force          X-Local',
     1'          Y-Local')
 1100 FORMAT(I5,6(2X,G15.6))
 1105 FORMAT('       ---------------  ---------------'/' Totals',
     1G15.6,2X,G15.6,2X,G15.6)
 1110 FORMAT(//' Gauss point stresses and and other state variables'/
     1         ' ==================================================')
 4440 FORMAT(//'UPDATED RESULTS'/
     1         ' ==================================================')
 1115 FORMAT(//' Element number',I5)
 1150 FORMAT(/' Gauss point ',I2,6X,' X-Coord= ',G11.4,
     1        ' Y-Coord= ',G11.4)
 1152 FORMAT(/' Gauss point ',I2,6X,' R-Coord= ',G11.4,
     1        ' Z-Coord= ',G11.4)
 1153 FORMAT(/' Gauss point ',I2,6X,' X-Coord= ',G11.4,
     1        ' R-Coord= ',G11.4)
 1154 FORMAT(/' Gauss point ',I2,6X,' X-Coord= ',G11.4,
     1        ' Y-Coord= ',G11.4,' Z-Coord= ',G11.4)
 1160 FORMAT(' S-xx  = ',G12.4,' S-yy  = ',G12.4,' S-xy  = ',G12.4)
 1161 FORMAT(' S-xx  = ',G12.4,' S-yy  = ',G12.4,' S-xy  = ',G12.4,
     1       ' S-zz  = ',G12.4)
 1162 FORMAT(' S-rr  = ',G12.4,' S-zz  = ',G12.4,' S-rz  = ',G12.4,
     1       ' S-h   = ',G12.4)
 1163 FORMAT(' S-xx  = ',G12.4,' S-rr  = ',G12.4,' S-xr  = ',G12.4,
     1       ' S-h   = ',G12.4)
 1164 FORMAT(' S-max = ',G12.4,' S-min = ',G12.4,' Angle = ',G12.4)
 1165 FORMAT(' Thick = ',G12.4)
 1166 FORMAT(' S-xx  = ',G12.4,' S-yy  = ',G12.4,' S-zz  = ',G12.4/
     1       ' S-xy  = ',G12.4,' S-yz  = ',G12.4,' S-xz  = ',G12.4)
 1167 FORMAT(' S-max = ',G12.4,' S-mid = ',G12.4,' S-min = ',G12.4)
C
C =================================
C Total macroscopic deformation gradient at current increment
      IF(NLARGE==1)THEN
        F=R0
C 2-D analysis
        IF((NTYPE.EQ.1).OR.(NTYPE.EQ.2).OR.(NTYPE.EQ.3))THEN
          F(1,1)=R1+TFACT*(DRV_STRAIN(1)-R1)
          F(1,2)=TFACT*DRV_STRAIN(2)
          F(2,1)=TFACT*DRV_STRAIN(3)
          F(2,2)=R1+TFACT*(DRV_STRAIN(4)-R1)
          F(3,3)=R1
C 3-D analysis
        ELSEIF(NTYPE.EQ.4)THEN
          F(1,1)=R1+TFACT*(DRV_STRAIN(1)-R1)
          F(1,2)=TFACT*DRV_STRAIN(2)
          F(1,3)=TFACT*DRV_STRAIN(3)
          F(2,1)=TFACT*DRV_STRAIN(4)
          F(2,2)=R1+TFACT*(DRV_STRAIN(5)-R1)
          F(2,3)=TFACT*DRV_STRAIN(6)
          F(3,1)=TFACT*DRV_STRAIN(7)
          F(3,2)=TFACT*DRV_STRAIN(8)
          F(3,3)=R1+TFACT*(DRV_STRAIN(9)-R1)
        ENDIF
C Jacobian
        DETJAC=DETM23(   3      ,F       ,NDIME       ,.FALSE.       )
      ENDIF
C P = Fbound (x) Ybound (original coordinates)
C Loop through boundary nodes (plus, minus, corner), listing dofs and
C multiplying by coordinates
C      
C Plus nodes
      SIGAUX=R0
      DO IND=IRV_NDSPPT(3),IRV_NDSPPT(4)
        INODE=IRV_NDSP(IND)
        DO I=1,NDOFN
          IDOF=(INODE-1)*NDOFN+I
          DO J=1,NDOFN
            SIGAUX(I,J)=SIGAUX(I,J)+STFOR(IDOF)*COORD(J,INODE,0)
          ENDDO
        ENDDO
      ENDDO
C Minus nodes
      DO IND=IRV_NDSPPT(5),IRV_NDSPPT(6)
        INODE=IRV_NDSP(IND)
        DO I=1,NDOFN
          IDOF=(INODE-1)*NDOFN+I
          DO J=1,NDOFN
            SIGAUX(I,J)=SIGAUX(I,J)+STFOR(IDOF)*COORD(J,INODE,0)
          ENDDO
        ENDDO
      ENDDO
C Corner nodes
      DO IND=IRV_NDSPPT(7),IRV_NDSPPT(8)
        INODE=IRV_NDSP(IND)
        DO I=1,NDOFN
          IDOF=(INODE-1)*NDOFN+I
          DO J=1,NDOFN
            SIGAUX(I,J)=SIGAUX(I,J)+STFOR(IDOF)*COORD(J,INODE,0)
          ENDDO
        ENDDO
      ENDDO
      
C
C Small strains: Cauchy stress tensor is the symmetric part of SIGAUX
C divided by the RVE volume
      IF(NLARGE==0)THEN
        SIGMAB=(SIGAUX+TRANSPOSE(SIGAUX))/(R2*DRV_CELVOL)        
C
C Large strains: SIGAUX divided by the RVE volume is the first Piola-
C Kirchhoff stress. Convert it to Cauchy stress:
C                                     T
C                Sigma = (1/J) * P * F
C
      ELSEIF(NLARGE.EQ.1)THEN      
        SIGMAB=MATMUL(SIGAUX,TRANSPOSE(F))/DETJAC/DRV_CELVOL
      ENDIF

      SELECT CASE (NTYPE)
      CASE (1,2,3)
        NDIM = 2
      CASE(4)
        NDIM = 3
      END SELECT
      CALL SYMTA(   SIGMAB       ,STRSH       ,NDIM       ,.TRUE.)
C
C Computation of principal stresses 
C ---------------------------------
C 2-D analysis
      IF((NTYPE==1).OR.(NTYPE==2).OR.(NTYPE==3))THEN
        CALL PRINC2(PSTRS,STRSH(1:3))
C 3-D analysis
      ELSEIF(NTYPE==4)THEN
        CALL PRINC3(PSTRS,STRSH(1:6))
      ENDIF
C       
C Computation of effective von Mises stresses
C -------------------------------------------
C
      S1=PSTRS(1)
      S2=PSTRS(2)
      IF(NTYPE.EQ.1)THEN
C Plane stress: third principal stress is S-zz = 0
        S3=R0
      ELSEIF(NTYPE.EQ.2)THEN
C Plane strain: third principal stress is S-zz (out of plane)
        S3=STRSH(4)
C        S3=0.3*(S1+S2)
      ELSEIF(NTYPE.EQ.4)THEN
        S3=PSTRS(3)
      ENDIF
      EFFST=SQRT(((S1-S2)**2+(S2-S3)**2+(S3-S1)**2)/R2)
C
C Print and display the homogenised and effective stresses
C --------------------------------------------------------		
c       WRITE(*, 2000) IINCS 
c       WRITE(16,2000) IINCS
C 2-D analysis
       IF((NTYPE==1).OR.(NTYPE==2).OR.(NTYPE==3))THEN
         IF(NLARGE.EQ.0)THEN
           WRITE(*, 2010)TFACT*DRV_STRAIN(1:3)
           WRITE(16,2010)TFACT*DRV_STRAIN(1:3)
           VOLSTRA=TFACT*DRV_STRAIN(1)+TFACT*DRV_STRAIN(2)
         ELSEIF(NLARGE.EQ.1)THEN
           WRITE(*, 2015)F(1,1:2),F(2,1:2)
           WRITE(16,2015)F(1,2)
         ENDIF
          WRITE(*, 2020) STRSH(1:3)
          WRITE(16,2020) STRSH(3)
          WRITE(*, 2025) PSTRS(1:3)
c          WRITE(16,2025) PSTRS(1:3)
C 3-D analysis
       ELSEIF(NTYPE==4)THEN
         IF(NLARGE.EQ.0)THEN
           WRITE(*, 2012)TFACT*DRV_STRAIN(1:6)
           WRITE(16,2012)TFACT*DRV_STRAIN(1:6)
           VOLSTRA=TFACT*DRV_STRAIN(1)+TFACT*DRV_STRAIN(2)+
     1             TFACT*DRV_STRAIN(3)
         ELSEIF(NLARGE.EQ.1)THEN
           WRITE(*, 2017)F(1,:),F(2,:),F(3,:)
           WRITE(16,2017)F(1,:),F(2,:),F(3,:)
         ENDIF
         WRITE(*, 2022) STRSH(1:6)
         WRITE(16,2022) STRSH(1:6)
         WRITE(*, 2027) PSTRS(1:3)
         WRITE(16,2027) PSTRS(1:3)
       ENDIF
C
      WRITE(*, 2030) EFFST
c      WRITE(16,2030) EFFST
C
CI HAVE ADDED OUTPUT HERE
C
C Stresses and other state and algorithmic variables at gauss points
C ==================================================================
C
C      IF(N3.NE.1)THEN
        WRITE(16,1110)
        STRSF=0
        RSTAVF=R0
        DEP=R0
        DP=R0
        WORKDP=R0
        WORKGOD=R0
        VOLUMEE=R0
        STRAINWO=R0
        DO 120 IELEM=1,NELEM
          IGRUP=IGRPID(IELEM)
          IELIDN=IELTID(IGRUP)
          IELTYP=IELPRP(1,IELIDN)
          NNODE =IELPRP(3,IELIDN)
          NGAUSP=IELPRP(4,IELIDN)
C
          WRITE(16,1115)IELEM
C Evaluate Gauss point coordinates
          DO 91 INODE=1,NNODE
            LNODE=IABS(LNODS(IELEM,INODE))
            DO 90 IDIME=1,NDIME
              ELCOD(IDIME,INODE)=COORD(IDIME,LNODE,1)
   90       CONTINUE
   91     CONTINUE
          IF(NTYPE.EQ.1)THEN
            LSTRE=3
          ELSEIF(NTYPE.EQ.2)THEN
            LSTRE=4
          ELSEIF(NTYPE.EQ.3)THEN
            LSTRE=4
          ELSEIF(NTYPE.EQ.4)THEN
            LSTRE=6
          ENDIF
          DO 110 IGAUSP=1,NGAUSP
            IPPOS=NDIME*(IGAUSP-1)
C Set Gauss points coordinates
            DO IDIME=1,NDIME
              EISCRD(IDIME)=RELPRP(IPPOS+IDIME,IELIDN)
            ENDDO
C
            CALL SHPFUN
     1(   DERIV      ,EISCRD      ,0          ,IELTYP     ,
     2    NDIME      ,SHAPE      )
C
            CALL JACDET
     1(   CARTD      ,DERIV      ,DETJAC     ,ELCOD      ,IELEM      ,
     2    NDIME      ,NDIME      ,NNODE      )
            CALL GETGCO
     1(   GPCOD      ,ELCOD      ,MDIME      ,NDIME      ,NNODE      ,
     2    SHAPE      )
C Output gauss points stresses (common to all materials)
C ------------------------------------------------------
            IF(NTYPE.EQ.1)THEN
              WRITE(16,1150)IGAUSP,(GPCOD(I),I=1,NDIME)
              WRITE(16,1160)(STRSG(I,IGAUSP,IELEM,1),I=1,LSTRE)
            ELSEIF(NTYPE.EQ.2)THEN
              WRITE(16,1150)IGAUSP,(GPCOD(I),I=1,NDIME)
              WRITE(16,1161)(STRSG(I,IGAUSP,IELEM,1),I=1,LSTRE)
            ELSEIF(NTYPE.EQ.3)THEN
              IF(NAXIS.EQ.1)THEN
                WRITE(16,1152)IGAUSP,(GPCOD(I),I=1,NDIME)
                WRITE(16,1162)(STRSG(I,IGAUSP,IELEM,1),I=1,LSTRE)
              ELSE
                WRITE(16,1153)IGAUSP,(GPCOD(I),I=1,NDIME)
                WRITE(16,1163)(STRSG(I,IGAUSP,IELEM,1),I=1,LSTRE)
              ENDIF
            ELSEIF(NTYPE.EQ.4)THEN
              WRITE(16,1154)IGAUSP,(GPCOD(I),I=1,NDIME)
              WRITE(16,1166)(STRSG(I,IGAUSP,IELEM,1),I=1,LSTRE)
            ENDIF
C and principal stresses
C 2-D analysis
            IF((NTYPE.EQ.1).OR.(NTYPE.EQ.2).OR.(NTYPE.EQ.3))THEN
              CALL PRINC2(PSTRS,STRSG(1,IGAUSP,IELEM,1))
              WRITE(16,1164)(PSTRS(I),I=1,3)
C 3-D analysis
            ELSEIF(NTYPE.EQ.4)THEN
              CALL PRINC3(PSTRS,STRSG(1,IGAUSP,IELEM,1))
              WRITE(16,1167)(PSTRS(I),I=1,3)
            ENDIF
C
C output current thickness (for large strains in plane stress only)
            IF(NLARGE.EQ.1.AND.NTYPE.EQ.1)THEN
              WRITE(16,1165)THKGP(IGAUSP,IELEM,1)
            ENDIF
C Output other (material-specific) state and algorithmic variables
C ----------------------------------------------------------------
            CALL MATIOR
     1(   NTYPE                    ,IPROPS(1,MATTID(IGRPID(IELEM)))  ,
     2    RALGVA(1,IGAUSP,IELEM,1) ,RPROPS(1,MATTID(IGRPID(IELEM)))  ,
     3    RSTAVA(1,IGAUSP,IELEM,1) ,STRSG(1,IGAUSP,IELEM,1)          ,
     4    STREPG(1,IGAUSP,IELEM,1)                                   )
           DVOLU=DETJAC
           T_BAR=0.1-(0.1-0)*EXP(-100*RSTAVA(7,IGAUSP,IELEM,1))*1
           RSTAVF=RSTAVF+RSTAVA(7,IGAUSP,IELEM,1)*DVOLU
            DO I=1,LSTRE
             STRSF(I,1,1,1)= STRSF(I,1,1,1)+
     1                       STRSG(I,IGAUSP,IELEM,1)*DVOLU
            ENDDO
C            CCCCCCCCCCCCCCCCCCCC
c     1              STREPG(JJ,IGAUSP,IELEM,1)*STRSG(JJ,IGAUSP,IELEM,1)
C                    TFACT*DRV_STRAIN(JJ)*STRSG(JJ,IGAUSP,IELEM,1)
C            WORKDP=R0
            VOLUMEE=VOLUMEE+DVOLU
!             FFF(1)=F(1,1)
!             FFF(2)=F(2,2)
!             FFF(3)=F(3,3)
!             FFF(4)=F(1,2)
!             FFF(5)=F(1,3)
!             FFF(6)=F(2,3)
!             DO JJ=1,6
!              WORKDP=WORKDP+
!      1              TFACT*DRV_STRAIN(JJ)*STRSG(JJ,IGAUSP,IELEM,1)*DVOLU
!             ENDDO
C            WORKDP=WORKDP+WORKDP*DVOLU
C            WORKGOD=WORKGOD+WORKDP
!  CDISIPATION NEED TO BE HERE INSIDE GPANDELE LOOP!!!!!!!! 
             ETA=RPROPS(1,4)
             MU=ETA/1.73
             DENOM=SQRT(1+(ETA/1.73)**2)
             KK=490*1.73/(1+(MU/1.41))
             DP=R0
             DO II=1,6
! !              DO JJ=1,6
              DP=DP+STREPG(II,IGAUSP,IELEM,1)*STREPG(II,IGAUSP,IELEM,1)
! !              ENDDO 
            ENDDO
            DEP=DEP+(KK*SQRT(2*DP)/DENOM)*DVOLU
! CCCCCCCCCCCCCCCCCCCC
!           WORKDIS=R0
!           DO JJ=1,6
!             WORKDIS=WORKDIS+
!                     STREPG(II,IGAUSP,IELEM,1)*STRSG(I,IGAUSP,IELEM,1)
!           ENDDO
  110     CONTINUE
  120   CONTINUE
C      ENDIF
cccADDING DISSIOPATION
C         DO II=1,6
C          DP=SQRT(2*STREPG(II,IGAUSP,IELEM,1)*STREPG(II,IGAUSP,IELEM,1))
C         ENDDO
         
C         VOLSTRA=TFACT*DRV_STRAIN(1)+TFACT*DRV_STRAIN(2)+TFACT*DRV_STRAIN(3)
         STRAIND(1)=TFACT*DRV_STRAIN(1)-VOLSTRA/R3
         STRAIND(2)=TFACT*DRV_STRAIN(2)-VOLSTRA/R3
         STRAIND(3)=TFACT*DRV_STRAIN(3)-VOLSTRA/R3
         STRAIND(4)=TFACT*DRV_STRAIN(4)
         STRAIND(5)=TFACT*DRV_STRAIN(5)
         STRAIND(6)=TFACT*DRV_STRAIN(6)
         STRAINWO = TFACT*DRV_STRAIN(4)
         SD1=STRAIND(1)
         SD2=STRAIND(2)
         SD3=STRAIND(3)
         SD4=STRAIND(4)
         SD5=STRAIND(5)
         SD6=STRAIND(6)
         DEVSTR=SQRT(((SD1-SD2)**2+(SD2-SD3)**2+
     1       (SD3-SD1)**2+6*(SD4**2+SD5**2+SD6**2))/R2)
!          WORKGOD=WORKGOD+WORKDP 
!          HOMOWORK=R0
!          DO MM=1,6
!          HOMOWORK=HOMOWORK+STRSF(MM,1,1,1)*TFACT*DRV_STRAIN(MM)
!          ENDDO 
         IF(NTYPE.EQ.2)THEN
           S11=STRSF(1,1,1,1)
           S22=STRSF(2,1,1,1)
           S12=STRSF(3,1,1,1)
           S21=STRSMASH(1,2)
           S33=STRSF(4,1,1,1)
           WRITE(*,4440)  
           PRESSURE=(S11+S22+S33)/3
       EFFSTU=SQRT(((S11-S22)**2+(S22-S33)**2+(S33-S11)**2+6*S12**2)/R2)
           WRITE(*, 4000) S11,S22,S33,S12
           WRITE(*, 4001) PRESSURE,EFFSTU
           WRITE(*,*)'======================================'
           WRITE(18,4004) VOLSTRA,DEVSTR,PRESSURE,EFFSTU,DEP
c3d the strsf(3) is sigmazz!!
         ELSE
      S11=STRSF(1,1,1,1)
      S22=STRSF(2,1,1,1)
      S12=STRSF(4,1,1,1)
      S13=STRSF(5,1,1,1)
      S23=STRSF(6,1,1,1)
      S33=STRSF(3,1,1,1)
c      HOMO=MATMUL(STRSMASH,TRANSPOSE(F))/DRV_CELVOL
      WRITE(*,4440)  
      PRESSURE=(S11+S22+S33)/3
      EFFSTU=SQRT(((S11-S22)**2+(S22-S33)**2+
     1       (S33-S11)**2+6*(S12**2+S13**2+S23**2))/R2)
      WRITE(*, 4003) S11,S22,S33,S12,S13,S23
!       IF(RSTAVA(1,IGAUSP,IELEM,1).NE.0)THEN
!          WRITE(*, *) 'STOP'
!       END IF
      WRITE(*, 4001) PRESSURE,EFFSTU
      WRITE(*,*)'======================================'
      WRITE(18,4004) STRAINWO,DEVSTR,PRESSURE,EFFSTU,VOLUMEE
C      PRESSURE,EFFSTU,WORKGOD,DEP
C      WRITE(18,4009) VOLSTRA,PRESSURE,EFFSTU
C      WRITE(18,*) RSTAVF,T_BAR
      ENDIF
      RETURN
      END
CDOC END_SUBROUTINE RVOUTP
