CDOC BEGIN_SUBROUTINE NODAVE
CDOC Computes and prints averaged nodal values of Gauss-point variables.
CDOC
CDOC This routine computes and prints the averaged (smoothed) nodal
CDOC values of stresses and other state and algorithmic variables. The
CDOC smoothed value of each variable at a node is computed by firstly
CDOC extrapolating the value of the variable from the Gauss-points to
CDOC the node in question and then averaging the extrapolated values
CDOC obtained from all elements sharing that node. The values computed
CDOC here are used in HYPLAS only for the purpose of output.
CHST
CHST  M.E.Honnor, February 1987
CHST
CHST  E.de Souza Neto  July  1996
CHST
CHST  E.de Souza Neto  October 2008: Nodal principal stresses and angle
CHST                                 now computed from nodal averaged
CHST                                 stresses (previous version computed
CHST                                 nodal principal stresses and angle
CHST                                 by extrapolation from Gauss point
CHST                                 values followed by averaging).
CHST
CHST  D. de Bortoli     March  2015: 3-D case added
CHST
      SUBROUTINE NODAVE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C Hyplas database
      INCLUDE '../MAXDIM.INC'
      INCLUDE '../MATERIAL.INC'
      INCLUDE '../ELEMENTS.INC'
      INCLUDE '../GLBDBASE.INC'
      DIMENSION
     1    RALGN(MRALGV,MNODE),RSTAN(MRSTAV,MNODE),STRSN(MSTRE,MNODE) ,
     2    THKN(MNODE)
C      DIMENSION
C     1    STRSA(MSTRE,MGRUP,MPOIN)     ,PSTRA(3)                     ,
C     2    RALGA(MRALGV,MGRUP,MPOIN)    ,RSTAA(MRSTAV,MGRUP,MPOIN)    ,
C     3    THKA(MGRUP,MPOIN)
      DIMENSION
     1   PSTRA(3), THKA(MGRUP,MPOIN)
      REAL*8, ALLOCATABLE :: STRSA(:,:,:), RALGA(:,:,:), RSTAA(:,:,:)
     
      DATA R0   ,R1   /0.0D0,1.0D0/
C***********************************************************************
C PRINTS AVERAGED (SMOOTHED) NODAL STRESSES AND OTHER STATE AND
C ALGORITHMIC VARIABLES. THE SMOOTHED VARIABLE AT EACH NODE IS OBTAINED
C BY EXTRAPOLATING THE VALUE OF THE VARIABLE FROM THE GAUSS POINTS TO
C THE NODE AND THEN AVERAGING THE VALUES OBTAINED FROM ALL ELEMENTS
C SHARING THAT NODE.
C
C REFERENCE: E Hinton & JS Campbel. Local and global smoothing of
C            discontinuous finite element functions using a least
C            squares method. Int. J. Num. Meth. Engng., 8:461-480, 1974.
C            E Hinton & DRJ Owen. An introduction to finite element
C            computations. Pineridge Press, Swansea, 1979.
C***********************************************************************
 1000 FORMAT(//' Averaged nodal stresses and other state',
     1         ' variables for group number ',I2/
     1         ' =======================================',
     3         '==============================')
 1010 FORMAT(/' Node number ',I5,3X,' X-Coord= ',G11.4,
     1        ' Y-Coord= ',G11.4)
 1020 FORMAT(/' Node number ',I5,3X,' R-Coord= ',G11.4,
     1        ' Z-Coord= ',G11.4)
 1030 FORMAT(/' Node number ',I5,3X,' X-Coord= ',G11.4,
     1        ' R-Coord= ',G11.4)
 1035 FORMAT(/' Node number ',I5,3X,' X-Coord= ',G11.4,
     1        ' Y-Coord= ',G11.4,' Z-Coord= ',G11.4)
 1040 FORMAT(' S-xx  = ',G12.4,' S-yy  = ',G12.4,' S-xy  = ',G12.4)
 1050 FORMAT(' S-xx  = ',G12.4,' S-yy  = ',G12.4,' S-xy  = ',G12.4,
     1       ' S-zz  = ',G12.4)
 1060 FORMAT(' S-rr  = ',G12.4,' S-zz  = ',G12.4,' S-rz  = ',G12.4,
     1       ' S-h   = ',G12.4)
 1070 FORMAT(' S-xx  = ',G12.4,' S-rr  = ',G12.4,' S-xr  = ',G12.4,
     1       ' S-h   = ',G12.4)
 1080 FORMAT(' S-max = ',G12.4,' S-min = ',G12.4,' Angle = ',G12.4)
 1090 FORMAT(' Thick = ',G12.4)
 1100 FORMAT(' S-xx  = ',G12.4,' S-yy  = ',G12.4,' S-zz  = ',G12.4/
     1       ' S-xy  = ',G12.4,' S-yz  = ',G12.4,' S-xz  = ',G12.4)
 1110 FORMAT(' S-max = ',G12.4,' S-mid = ',G12.4,' S-min = ',G12.4)
C
      ALLOCATE( STRSA(MSTRE,MGRUP,MPOIN) )
      ALLOCATE( RALGA(MRALGV,MGRUP,MPOIN) )
      ALLOCATE( RSTAA(MRSTAV,MGRUP,MPOIN) )
C Loop over element groups
C ========================
      DO 90 IGRUP=1,NGRUP
        WRITE(16,1000)IGRUP
        DO 10 IPOIN=1,NPOIN
          CALL RVZERO(STRSA(1,IGRUP,IPOIN),MSTRE)
          CALL RVZERO(PSTRA,3)
          CALL RVZERO(RSTAA(1,IGRUP,IPOIN),MRSTAV)
          CALL RVZERO(RALGA(1,IGRUP,IPOIN),MRALGV)
          THKA(IGRUP,IPOIN)=R0
   10   CONTINUE
C Loop over elements
C ------------------
        DO 70 IELEM=1,NELEM
          LGRUP=IGRPID(IELEM)
          IF(LGRUP.NE.IGRUP)GOTO 70
          IELIDN=IELTID(IGRUP)
          NNODE =IELPRP(3,IELIDN)
          NGAUSP=IELPRP(4,IELIDN)
C Extrapolate stresses and other state and algorithmic variables from
C gauss points to nodes
          IPOS=NGAUSP*NDIME+NGAUSP+1
          CALL EXTNOD
     1(   RELPRP(IPOS,IELIDN),
     2    STRSG(1,1,IELEM,1) ,STRSN    ,MSTRE     ,NGAUSP     ,NNODE   )
          CALL EXTNOD
     1(   RELPRP(IPOS,IELIDN),
     2    RSTAVA(1,1,IELEM,1),RSTAN    ,MRSTAV    ,NGAUSP     ,NNODE   )
          CALL EXTNOD
     1(   RELPRP(IPOS,IELIDN),
     2    RALGVA(1,1,IELEM,1),RALGN    ,MRALGV    ,NGAUSP     ,NNODE   )
C thickness (for large strains in plane stress only)
          IF(NLARGE.EQ.1.AND.NTYPE.EQ.1)CALL EXTNOD
     1(   RELPRP(IPOS,IELIDN),
     2    THKGP(1,IELEM,1)   ,THKN     ,1         ,NGAUSP     ,NNODE   )
C Nodal averaging
          DO 60 INODE=1,NNODE
            IPOIN=IABS(LNODS(IELEM,INODE))
            R1DVAL=R1/DBLE(NVALEN(IPOIN,IGRUP))
            DO 30 ISTRE=1,MSTRE
              STRSA(ISTRE,IGRUP,IPOIN)=STRSA(ISTRE,IGRUP,IPOIN)+
     1                                 STRSN(ISTRE,INODE)*R1DVAL
   30       CONTINUE
            DO 40 INTV=1,MRSTAV
              RSTAA(INTV,IGRUP,IPOIN)=RSTAA(INTV,IGRUP,IPOIN)+
     1                                RSTAN(INTV,INODE)*R1DVAL
   40       CONTINUE
            DO 50 IALGV=1,MRALGV
              RALGA(IALGV,IGRUP,IPOIN)=RALGA(IALGV,IGRUP,IPOIN)+
     1                                 RALGN(IALGV,INODE)*R1DVAL
   50       CONTINUE
            IF(NLARGE.EQ.1.AND.NTYPE.EQ.1)
     1               THKA(IGRUP,IPOIN)=THKA(IGRUP,IPOIN)+
     2                                 THKN(INODE)*R1DVAL
   60     CONTINUE
   70   CONTINUE
C
C Output average stresses to results file (common to all materials)
C -----------------------------------------------------------------
C
        DO 80 IPOIN=1,NPOIN
          IF(NVALEN(IPOIN,IGRUP).EQ.0)GOTO 80
          IF(NTYPE.EQ.1)THEN
            LSTRE=3
            WRITE(16,1010)IPOIN,(COORD(I,IPOIN,1),I=1,NDIME)
            WRITE(16,1040)(STRSA(I,IGRUP,IPOIN),I=1,LSTRE)
          ELSEIF(NTYPE.EQ.2)THEN
            LSTRE=4
            WRITE(16,1010)IPOIN,(COORD(I,IPOIN,1),I=1,NDIME)
            WRITE(16,1050)(STRSA(I,IGRUP,IPOIN),I=1,LSTRE)
          ELSEIF(NTYPE.EQ.3)THEN
            LSTRE=4
            IF(NAXIS.EQ.1)THEN
              WRITE(16,1020)IPOIN,(COORD(I,IPOIN,1),I=1,NDIME)
              WRITE(16,1060)(STRSA(I,IGRUP,IPOIN),I=1,LSTRE)
            ELSE
              WRITE(16,1030)IPOIN,(COORD(I,IPOIN,1),I=1,NDIME)
              WRITE(16,1070)(STRSA(I,IGRUP,IPOIN),I=1,LSTRE)
            ENDIF
          ELSEIF(NTYPE.EQ.4)THEN
            LSTRE=6
            WRITE(16,1035)IPOIN,(COORD(I,IPOIN,1),I=1,NDIME)
            WRITE(16,1100)(STRSA(I,IGRUP,IPOIN),I=1,LSTRE)
          ENDIF
C Compute and output principal stresses (and angle in 2-D)
          IF((NTYPE.EQ.1).OR.(NTYPE.EQ.2).OR.(NTYPE.EQ.3))THEN
            CALL PRINC2(PSTRA,STRSA(1,IGRUP,IPOIN))
            WRITE(16,1080)(PSTRA(I),I=1,3)
          ELSEIF(NTYPE.EQ.4)THEN
            CALL PRINC3(PSTRA,STRSA(1,IGRUP,IPOIN))
            WRITE(16,1110)(PSTRA(I),I=1,3)
          ENDIF
C current thickness (for large strains in plane stress only)
          IF(NLARGE.EQ.1.AND.NTYPE.EQ.1)WRITE(16,1090)THKA(IGRUP,IPOIN)
C
C and other (material-specific) state and algorithmic variables
C -------------------------------------------------------------
C
          CALL MATIOR
     1(   NTYPE      ,IPROPS(1,MATTID(IGRUP))  ,RALGA(1,IGRUP,IPOIN)  ,
     2    RPROPS(1,MATTID(IGRUP))              ,RSTAA(1,IGRUP,IPOIN)  ,
     3    STRSA(1,IGRUP,IPOIN)                 )
C
   80   CONTINUE
   90 CONTINUE
      DEALLOCATE( STRSA )
      DEALLOCATE( RALGA )
      DEALLOCATE( RSTAA )
      RETURN
      END
CDOC END_SUBROUTINE NODAVE
