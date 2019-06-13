CDOC BEGIN_SUBROUTINE OUTPUT
CDOC Prints displacements, reactions and other variables to results file
CDOC
CDOC This routine prints the nodal displacements and reactions and
CDOC Gauss-point and nodal averaged values of stresses and other state
CDOC and algorithmic variables to the results file. The argument NOUTP
CDOC controls what results are to be printed out.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION TFACT  >  Total proportional load factor.
CDOC DOUBLE_PRECISION TTIME  >  Current time.
CDOC INTEGER          IINCS  >  Current load increment number.
CDOC INTEGER          IITER  >  Equilibrium iteration at which
CDOC C                          convergence was achieved in the current
CDOC C                          load step.
CDOC INTEGER          NOUTP  >  Load output flag that determines what
CDOC C                          results are to be printed out.
CDOC END_PARAMETERS
CDOC
      SUBROUTINE OUTPUT
     1(   TFACT      ,TTIME      ,IINCS      ,IITER      ,NOUTP      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C Hyplas database
      INCLUDE '../MAXDIM.INC'
      INCLUDE '../MATERIAL.INC'
      INCLUDE '../ELEMENTS.INC'
      INCLUDE '../GLBDBASE.INC'
C
      DIMENSION NOUTP(5)
C
      DIMENSION
     1    DERIV(MDIME,MNODE) ,ELCOD(MDIME,MNODE) ,GPCOD(MDIME)       ,
     2    PSTRS(3)           ,SHAPE(MNODE)       ,EISCRD(MDIME)
      DATA R0/0.0D0/
C
C***********************************************************************
C OUTPUTS DISPLACEMENTS, REACTIONS, STRESSES AND OTHER STATE AND
C ALGORITHMIC VARIABLES TO RESULTS FILE
C
C REFERENCE: Section 5.4.7
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
C
C Set output flags
C ================
C
      N1=NOUTP(1)
      N2=NOUTP(2)
      N3=NOUTP(3)
      N4=NOUTP(4)
      IF(NALGO.LT.0)THEN
        IF(N1.NE.0)THEN
          IF(MOD(IINCS,N1).EQ.0)THEN
            N1=1
          ELSE
            N1=0
          ENDIF
        ENDIF
        IF(N2.NE.0)THEN
          IF(MOD(IINCS,N2).EQ.0)THEN
            N2=1
          ELSE
            N2=0
          ENDIF
        ENDIF
        IF(N3.NE.0)THEN
          IF(MOD(IINCS,N3).EQ.0)THEN
            N3=1
          ELSE
            N3=0
          ENDIF
        ENDIF
        IF(N4.NE.0)THEN
          IF(MOD(IINCS,N4).EQ.0)THEN
            N4=1
          ELSE
            N4=0
          ENDIF
        ENDIF
      ENDIF
      IF(N1.EQ.0.AND.N2.EQ.0.AND.N3.EQ.0.AND.N4.EQ.0)RETURN
C
C
      WRITE(16,1020)IINCS,TTIME,TFACT,IITER
C
C Output displacements
C ====================
C
      IF(N1.NE.0)THEN
        WRITE(16,1030)
C 2-D analysis (axisymmetric)
        IF(NTYPE.EQ.3)THEN
          IF(NAXIS.EQ.1)THEN
            WRITE(16,1040)
          ELSE
            WRITE(16,1045)
          ENDIF
C 2-D analysis (plane stress/strain)
        ELSEIF((NTYPE.EQ.1).OR.(NTYPE.EQ.2))THEN
          WRITE(16,1050)
C 3-D analysis
        ELSEIF(NTYPE.EQ.4)THEN
          WRITE(16,1055)
        ENDIF
        DO 10 IPOIN=1,NPOIN
          NGASH=IPOIN*NDOFN
          NGISH=NGASH-NDOFN+1
          WRITE(16,1060)IPOIN,(TDISP(IGASH),IGASH=NGISH,NGASH)
   10   CONTINUE
      ENDIF
C
C Output reactions
C ================
C
      IF(N2.NE.0)THEN
        WRITE(16,1070)
        DO 45 IVFIX=1,NVFIX
          IF(ANGLE(IVFIX).NE.R0)THEN
            IF(NTYPE.EQ.3)THEN
              IF(NAXIS.EQ.1)THEN
                WRITE(16,1085)
              ELSE
                WRITE(16,1086)
              ENDIF
            ELSE
              WRITE(16,1095)
            ENDIF
            GOTO 47
          ENDIF
   45   CONTINUE
        IF(NTYPE.EQ.3)THEN
C 2-D analysis (axisymmetric)
          IF(NAXIS.EQ.1)THEN
            WRITE(16,1080)
          ELSE
            WRITE(16,1081)
          ENDIF
C 2-D analysis (plane stress/strain)
        ELSEIF((NTYPE.EQ.1).OR.(NTYPE.EQ.2))THEN
          WRITE(16,1090)
C 3-D analysis
        ELSEIF(NTYPE.EQ.4)THEN
          WRITE(16,1091)
        ENDIF
   47   CONTINUE
        TRX=R0
        TRY=R0
        TRZ=R0
        DO 70 IPOIN=1,NPOIN
          ISVAB=(IPOIN-1)*NDOFN
          DO 50 IVFIX=1,NVFIX
            IF(NOFIX(IVFIX).EQ.IPOIN)GOTO 60
   50     CONTINUE
          GOTO 70
   60     CONTINUE
          IF(ANGLE(IVFIX).NE.R0)THEN
            C=COS(ANGLE(IVFIX))
            S=SIN(ANGLE(IVFIX))
            GASHI= C*TREAC(IVFIX,1)+S*TREAC(IVFIX,2)
            GASHJ=-S*TREAC(IVFIX,1)+C*TREAC(IVFIX,2)
            IF(IFFIX(ISVAB+1).EQ.0)GASHI=R0
            IF(IFFIX(ISVAB+2).EQ.0)GASHJ=R0
            WRITE(16,1100)IPOIN,(TREAC(IVFIX,IDOFN),IDOFN=1,NDOFN),
     1                    GASHI,GASHJ
          ELSE
            WRITE(16,1100)IPOIN,(TREAC(IVFIX,IDOFN),IDOFN=1,NDOFN)
          ENDIF
          TRX=TRX+TREAC(IVFIX,1)
          TRY=TRY+TREAC(IVFIX,2)
          IF(NTYPE.EQ.4)TRZ=TRZ+TREAC(IVFIX,3)
   70   CONTINUE
        IF(NTYPE.EQ.4)THEN
          WRITE(16,1105)TRX,TRY,TRZ
        ELSE
          WRITE(16,1105)TRX,TRY
        ENDIF
      ENDIF
C
C Stresses and other state and algorithmic variables at gauss points
C ==================================================================
C
      IF(N3.NE.0)THEN
        WRITE(16,1110)
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
     2    MDIME      ,SHAPE      )
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
     3    RSTAVA(1,IGAUSP,IELEM,1) ,STRSG(1,IGAUSP,IELEM,1)          )
C
  110     CONTINUE
  120   CONTINUE
      ENDIF
C
C Stresses and other state and algorithmic variables at nodes
C ===========================================================
C
      IF(N4.NE.0)THEN
        CALL NODAVE
      ENDIF      
C
      RETURN
      END
CDOC END_SUBROUTINE OUTPUT
