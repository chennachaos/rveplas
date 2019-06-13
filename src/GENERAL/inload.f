CDOC BEGIN_SUBROUTINE INLOAD
CDOC Reads external loads and assembles the external force vector.
CDOC
CDOC This routine reads the prescribed body forces, surface tractions
CDOC and point (nodal) loads from the input data file and assembles the
CDOC corresponding global external force vector.
CDOC
CHST
CHST E. de Souza Neto, April 2011: I/O error messages added
CHST
CHST    D. de Bortoli, March 2015: Added 3-D load cases (point and
CHST                               pressure loads);
CHST                               Changed 'JACOB2' to 'JACDET' in
CHST                               function calls;
CHST
      SUBROUTINE INLOAD
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C Hyplas global database
      INCLUDE '../MAXDIM.INC'
      INCLUDE '../MATERIAL.INC'
      INCLUDE '../ELEMENTS.INC'
      INCLUDE '../GLBDBASE.INC'
C
      LOGICAL FOUND
      CHARACTER*80 INLINE
      DIMENSION
     1    SHAPE(MNODE)       ,DERIV(MDIME,MNODE) ,CARTD(MDIME,MNODE) ,
     2    GPCOD(MDIME)       ,EISCRD(MDIME)
      DIMENSION
     1    GRVDIR(3)
      DIMENSION
     1    DGASH(MDOFN)       ,EGASH(MDOFN)       ,EXTNOR(MDOFN)      ,
     2    ELCOD(MDIME,MNODE) ,NOPRS(MNODE)       ,PGASH(MDOFN)       ,
     3    POINT(MDOFN)       ,PRESS(MNODE,MDOFN)
      DIMENSION
     1    IWBEG(40)          ,IWEND(40)          ,NODCHK(MNODE)      ,
     2    NODAUX(MPOIN)      ,THKN(MNODE)
C      
      DATA R0,R1,R8,R45/0.0D0,1.0D0,8.0D0,45.0D0/
C***********************************************************************
C READS EXTERNAL LOADINGS (BODY FORCE AND SURFACE TRACTIONS) FROM INPUT
C DATA FILE AND ASSEMBLES THE GLOBAL EXTERNAL FORCE VECTOR
C
C REFERENCE: Figure 5.1
C            Section 5.3.3
C***********************************************************************
 1040 FORMAT(///
     1' Loading specification (other than prescribed displacements)'/
     1' ==========================================================='//
     2'    If any of the flags below is set to 1,  then'/
     3'    the corresponding type of loading is applied'/
     4'    to the structure.'//
     5' Point loading flag ...........................=',I3/
     6' Gravity loading flag .........................=',I3/
     7' Distributed edge loading flag ................=',I3)
 1050 FORMAT(//' Point load applied in',I6,' nodes'/
     1         ' ---------------------------------'//
     2         ' Node        X-Component    Y-Component')
 1055 FORMAT(//' Point load applied in',I6,' nodes'/
     1         ' ---------------------------------'//
     2         ' Node        R-Component    Z-Component')
 1060 FORMAT(//' Point load applied in',I6,' nodes'/
     1         ' ---------------------------------'//
     2         ' Node        X-Component    Y-Component    Z-Component')
 1070 FORMAT(I5,5X,3G15.6)
 1080 FORMAT(//' Gravity load'/
     1         ' ------------'//
     2         ' Gravity angle (degrees) =',G15.6/
     3         ' Gravity constant .......=',G15.6)
 1085 FORMAT(//' Gravity load'/
     1         ' ------------'//
     2         ' Gravity direction (not normalised):'/
     3         '     Component X  .......=',G15.6/
     4         '     Component Y  .......=',G15.6/
     5         '     Component Z  .......=',G15.6/
     6         ' Gravity constant .......=',G15.6)
 1110 FORMAT(//' Edge load applied in ',I6,' edges'/
     1         ' ---------------------------------')
 1130 FORMAT(/' Element Number =',I5/' Node      X-Coord.  ',
     1'     Y-Coord.       Norm. Load         Tang. Load')
 1135 FORMAT(/' Element Number =',I5/' Node      R-Coord.  ',
     1'     Z-Coord.       Norm. Load         Tang. Load')
 1137 FORMAT(/' Element Number =',I5/' Node      X-Coord.  ',
     1'     Y-Coord.       Z-Coord.        Norm. Load',
     2'     Tang. Load 1   Tang. Load 2')
 1140 FORMAT(I5,2X,2G15.6,3X,2G15.6)
 1147 FORMAT(I5,2X,3G15.6,3X,3G15.6)
C
C
      IF(NTYPE.EQ.3)TWOPI=R8*ATAN(R1)
C
C Initialize load vector of all elements
C ======================================
C
      DO 10 IELEM=1,NELEM
        IGRUP=IGRPID(IELEM)
        IELIDN=IELTID(IGRUP)
        NEVAB=IELPRP(5,IELIDN)
        CALL RVZERO(RLOAD(1,IELEM),NEVAB)
   10 CONTINUE
C
C Read data controlling loading types to be inputted
C ==================================================
C
      IPLOD=0
      IGRAV=0
      IEDGE=0
C
      CALL FNDKEY
     1(   FOUND      ,IWBEG      ,IWEND      ,'LOADINGS',
     2    INLINE     ,15         ,NWRD       )
      IF(.NOT.FOUND)CALL ERRPRT('ED0092')
      NLOADS=NWRD-1
      IF(NLOADS.NE.0)THEN
        DO 12 I=1,NLOADS
          IF(INLINE(IWBEG(1+I):IWEND(1+I)).EQ.'POINT')THEN
            IPLOD=1
          ELSEIF(INLINE(IWBEG(1+I):IWEND(1+I)).EQ.'EDGE')THEN 
            IEDGE=1
          ELSEIF(INLINE(IWBEG(1+I):IWEND(1+I)).EQ.'GRAVITY')THEN
            IGRAV=1
          ELSEIF((INLINE(IWBEG(1+I):IWEND(1+I)).EQ.'0').OR.
     1           (INLINE(IWBEG(1+I):IWEND(1+I)).EQ.'NONE'))THEN
            CONTINUE
          ELSE
            CALL ERRPRT('ED0030')
          ENDIF
   12   CONTINUE
      ENDIF
C
      WRITE(16,1040)IPLOD,IGRAV,IEDGE
C
C Read nodal point loads
C ======================
C
      IF(IPLOD.NE.0)THEN
        CALL FNDKEY
     1(   FOUND      ,IWBEG      ,IWEND      ,'POINT_LOADS',
     2    INLINE     ,15         ,NWRD       )
        IF(.NOT.FOUND)CALL ERRPRT('ED0093')
        IF(NWRD.EQ.1)CALL ERRPRT('ED0031')
        NPLOAD=INTNUM(INLINE(IWBEG(2):IWEND(2)))
        IF(NTYPE.EQ.1.OR.NTYPE.EQ.2)THEN
          WRITE(16,1050)NPLOAD
        ELSEIF(NTYPE.EQ.3)THEN
          WRITE(16,1055)NPLOAD
        ELSEIF(NTYPE.EQ.4)THEN
          WRITE(16,1060)NPLOAD
        ENDIF
        DO 55 IPLOAD=1,NPLOAD
          READ(15,*,ERR=900,END=900)LODPT,(POINT(IDOFN),IDOFN=1,NDOFN)
          WRITE(16,1070)LODPT,(POINT(IDOFN),IDOFN=1,NDOFN)
          IF(LODPT.LE.0.OR.LODPT.GT.NPOIN)CALL ERRPRT('ED0134')
C
C Associate the nodal point loads with an element
C
          DO 35 IELEM=1,NELEM
            IGRUP=IGRPID(IELEM)
            IELIDN=IELTID(IGRUP)
            NNODE=IELPRP(3,IELIDN)
            DO 30 INODE=1,NNODE
              NLOCA=IABS(LNODS(IELEM,INODE))
              IF(LODPT.EQ.NLOCA)GOTO 40
   30       CONTINUE
   35     CONTINUE
   40     CONTINUE
          DO 50 IDOFN=1,NDOFN
            NGASH=(INODE-1)*NDOFN+IDOFN
            RLOAD(NGASH,IELEM)=RLOAD(NGASH,IELEM)+POINT(IDOFN)
   50     CONTINUE
   55   CONTINUE
      ENDIF
C
C Gravity loading
C ===============
C
      IF(IGRAV.NE.0)THEN
C
        CALL FNDKEY
     1(   FOUND      ,IWBEG      ,IWEND      ,'GRAVITY_LOAD',
     2    INLINE     ,15         ,NWRD       )
        IF(.NOT.FOUND)CALL ERRPRT('ED0094')
        IF(NTYPE.EQ.4)THEN
C 3-D: Read gravity direction vector and gravitational constant
          READ(15,*,ERR=910,END=910)(GRVDIR(IDIME),IDIME=1,3),GRAVY
          WRITE(16,1085)(GRVDIR(IDIME),IDIME=1,3),GRAVY
C Normalise gravity direction vector
          IF(NORM2(GRVDIR).EQ.R0)THEN
            CALL ERRPRT('ED0133')
          ENDIF
          GRVDIR=GRVDIR/NORM2(GRVDIR)
        ELSE
C 2-D: Read gravity angle and gravitational constant
          READ(15,*,ERR=910,END=910)THETA,GRAVY
          WRITE(16,1080)THETA,GRAVY
          THETA=THETA*ATAN(R1)/R45
        ENDIF
C
C Loop over elements
C
        DO 90 IELEM=1,NELEM
          IGRUP =IGRPID(IELEM)
          IELIDN=IELTID(IGRUP)
          IELTYP=IELPRP(1,IELIDN)
          NNODE =IELPRP(3,IELIDN)
          NGAUSP=IELPRP(4,IELIDN)
C Set up preliminary constants
          MATIDN=MATTID(IGRPID(IELEM))
          DENSE=RPROPS(1,MATIDN)
          IF(DENSE.EQ.R0) GOTO 90
          IF(NTYPE.EQ.4)THEN
            GXCOM=DENSE*GRAVY*GRVDIR(1)
            GYCOM=DENSE*GRAVY*GRVDIR(2)
            GZCOM=DENSE*GRAVY*GRVDIR(3)
          ELSE
            GXCOM=DENSE*GRAVY*SIN(THETA)
            GYCOM=-DENSE*GRAVY*COS(THETA)
          ENDIF
C Coordinates of the element nodal points
          DO 65 INODE=1,NNODE
            LNODE=IABS(LNODS(IELEM,INODE))
            DO 60 IDIME=1,NDIME
              ELCOD(IDIME,INODE)=COORD(IDIME,LNODE,1)
   60       CONTINUE
   65     CONTINUE
C
C Loop for numerical integration over element domain
C
          IPWEI=NGAUSP*NDIME
          DO 85 IGAUSP=1,NGAUSP
            IPPOS=NDIME*(IGAUSP-1)
            DO IDIME=1,NDIME
              EISCRD(IDIME)=RELPRP(IPPOS+IDIME,IELIDN)
            ENDDO
            WEIGP=RELPRP(IPWEI+IGAUSP,IELIDN)
C
C Compute the shape functions at the sampling points and elemental
C volume
            CALL SHPFUN
     1(   DERIV      ,EISCRD      ,0          ,IELTYP     ,
     2    MDIME      ,SHAPE      )
            CALL JACDET
     1(   CARTD      ,DERIV      ,DETJAC     ,ELCOD      ,IELEM      ,
     2    MDIME      ,NDIME      ,NNODE      )
            CALL GETGCO
     1(   GPCOD      ,ELCOD      ,MDIME      ,NDIME      ,NNODE      ,
     2    SHAPE      )
C
            DVOLU=DETJAC*WEIGP
            IF(NTYPE.EQ.1)THEN
              DVOLU=DVOLU*THKGP(IGAUSP,IELEM,1)
            ELSEIF(NTYPE.EQ.3)THEN
              DVOLU=DVOLU*TWOPI*GPCOD(NAXIS)
            ENDIF
C Calculate equivalent nodal loads and add them to element force vector
            DO 70 INODE=1,NNODE
              NGASH=(INODE-1)*NDOFN+1
              MGASH=(INODE-1)*NDOFN+2
              RLOAD(NGASH,IELEM)=RLOAD(NGASH,IELEM)+
     1                           GXCOM*SHAPE(INODE)*DVOLU
              RLOAD(MGASH,IELEM)=RLOAD(MGASH,IELEM)+
     1                           GYCOM*SHAPE(INODE)*DVOLU
              IF(NTYPE.EQ.4)THEN
                LGASH=(INODE-1)*NDOFN+3
                RLOAD(LGASH,IELEM)=RLOAD(LGASH,IELEM)+
     1                             GZCOM*SHAPE(INODE)*DVOLU
              ENDIF
   70       CONTINUE
   85     CONTINUE
   90   CONTINUE
      ENDIF
C
C Distributed edge loads (pressure)
C =================================
C
      IF(IEDGE.NE.0)THEN
        CALL FNDKEY
     1(   FOUND      ,IWBEG      ,IWEND      ,'EDGE_LOADS',
     2    INLINE     ,15         ,NWRD       )
        IF(.NOT.FOUND)CALL ERRPRT('ED0095')
        IF(NWRD.EQ.1)CALL ERRPRT('ED0032')
        NLOADE=INTNUM(INLINE(IWBEG(2):IWEND(2)))
        WRITE(16,1110)NLOADE
C
C Loop over loaded edges
C
        DO 160 ILOADE=1,NLOADE
C Read and echo the element number and corresponding global node numbers
C with prescribed pressure
          READ(15,*,ERR=920,END=920)
     1                      IELEM,NNODEG,(NOPRS(INODEG),INODEG=1,NNODEG)
          IF((NTYPE.EQ.1).OR.(NTYPE.EQ.2))THEN
            WRITE(16,1130)IELEM
          ELSEIF(NTYPE.EQ.3)THEN
            WRITE(16,1135)IELEM
          ELSEIF(NTYPE.EQ.4)THEN
            WRITE(16,1137)IELEM
          ENDIF
          IF(IELEM.LE.0.OR.IELEM.GT.NELEM)CALL ERRPRT('ED0019')
          DO 80 INODEG=1,NNODEG
            IPOIN=NOPRS(INODEG)
            IF(IPOIN.LE.0.OR.IPOIN.GT.NPOIN)CALL ERRPRT('ED0020')
   80     CONTINUE 
C Set properties of the current element
          IGRUP=IGRPID(IELEM)
          IELIDN=IELTID(IGRUP)
          IELTYP=IELPRP(1,IELIDN)
          NNODE =IELPRP(3,IELIDN)
          NGAUSP=IELPRP(4,IELIDN)
          NEDGEL=IELPRP(6,IELIDN)
          MNODEG=IELPRP(7,IELIDN)
          NGAUSB=IELPRP(8,IELIDN)
          IPOS=9
C Read and echo pressures
          READ(15,*,ERR=930,END=930)
     1             ((PRESS(INODEG,IDOFN),INODEG=1,NNODEG),IDOFN=1,NDOFN)
          DO 95 INODEG=1,NNODEG
            IPOIN=NOPRS(INODEG)
            IF(NTYPE.NE.4)THEN
              WRITE(16,1140)IPOIN,(COORD(I,IPOIN,1),I=1,NDIME),
     1                                       (PRESS(INODEG,I),I=1,NDIME)
            ELSE
              WRITE(16,1147)IPOIN,(COORD(I,IPOIN,1),I=1,NDIME),
     1                                       (PRESS(INODEG,I),I=1,NDIME)
            ENDIF
   95     CONTINUE
          IF(NNODEG.GT.MNODEG)CALL ERRPRT('ED0011')
C Check that global node numbers supplied correspond exactly to an edge
C of the current element
          DO 96 INODE=1,NNODE
            NODCHK(INODE)=0
   96     CONTINUE
          DO 98 INODEG=1,NNODEG
            DO 97 INODE=1,NNODE
              IPOIN=IABS(LNODS(IELEM,INODE))
              IF(IPOIN.EQ.NOPRS(INODEG))NODCHK(INODE)=1
   97       CONTINUE
   98     CONTINUE
C
          CALL CHKNDB
     1(   FOUND    ,NNODE    ,NEDGEL   ,NODCHK   ,IELPRP(IPOS,IELIDN))
          IF(.NOT.FOUND)CALL ERRPRT('ED0012')
C
C Get the global coordinates of the nodes of the loaded edge
          DO 104 INODEG=1,NNODEG
            IPOIN=NOPRS(INODEG)
            DO 100 IDIME=1,NDIME
              ELCOD(IDIME,INODEG)=COORD(IDIME,IPOIN,1)
  100       CONTINUE
C
            DO 102 II=1,NNODEG
              IF(IABS(LNODS(IELEM,NODCHK(II))).EQ.IPOIN)NODAUX(IPOIN)=II
  102       CONTINUE
  104     CONTINUE
C Extrapolate thickness to nodes (for plane stress only)
          IF(NTYPE.EQ.1)THEN
            IPOS=NGAUSP*NDIME+NGAUSP+1
            CALL EXTNOD
     1(   RELPRP(IPOS,IELIDN),
     2    THKGP(1,IELEM,1)   ,THKN     ,1         ,NGAUSP     ,NNODE   )
          ENDIF
C
C Loop for (boundary) numerical integration over loaded edge
C
C Set number of boundary dimensions (surfaces, for 3-D; curves, for 2-D)
          IF(NTYPE.EQ.4)THEN
            NDIMEB=2
          ELSE
            NDIMEB=1
          ENDIF
C
          DO 150 IGAUSB=1,NGAUSB
C Evaluate the shape functions at the boundary sampling points
            IPPOS=NGAUSP*NDIME+NGAUSP+NGAUSP*NNODE+(IGAUSB-1)*NDIMEB+1
            DO IDIMEB=1,NDIMEB
              EISCRD(IDIMEB)=RELPRP(IPPOS,IELIDN)
              IPPOS=IPPOS+1
            ENDDO
            IPWEI=NGAUSP*NDIME+NGAUSP+NGAUSP*NNODE+NGAUSB*NDIMEB+1
            WEIGPB=RELPRP(IPWEI-1+IGAUSB,IELIDN)
C
            CALL SHPFUN
     1(   DERIV      ,EISCRD     ,1          ,IELTYP     ,
     2    MDIME      ,SHAPE      )
C Calculate components of the equivalent nodal loads
            DO 114 IDOFN=1,NDOFN
              PGASH(IDOFN)=R0
              DGASH(IDOFN)=R0
              EGASH(IDOFN)=R0
C
              DO 110 INODEG=1,NNODEG
C - NOPRS(INODEG) has global node number where pressure is applied
C   (local node INODEG)
C - NODAUX(IPOIN) has local number (on element boundary) of global node 
C   IPOIN
                II=NODAUX(NOPRS(INODEG))
C Pressure components interpolated at current Gauss point using element
C boundary shape functions
                PGASH(IDOFN)=PGASH(IDOFN)+
     1                       PRESS(INODEG,IDOFN)*SHAPE(II)
C First tangent vector: derivative of position with respect to boundary 
C XI coordinate
                DGASH(IDOFN)=DGASH(IDOFN)+
     1                       ELCOD(IDOFN,INODEG)*DERIV(1,II)
C Second tangent vector (only defined in 3-D): derivative of position
C with respect to boundary ETA coordinate
                IF(NTYPE.EQ.4)THEN
                  EGASH(IDOFN)=EGASH(IDOFN)+
     1                          ELCOD(IDOFN,INODEG)*DERIV(2,II)
                ENDIF
  110         CONTINUE
  114       CONTINUE
C
C Calculate components of interpolated Gauss point pressures in global 
C X,Y (and Z) coordinate system
            IF(NTYPE.EQ.4)THEN
C In 3-D, external normal vector (not unit!) is the cross product of the
C two tangent vectors
              EXTNOR(1)=DGASH(2)*EGASH(3)-DGASH(3)*EGASH(2)
              EXTNOR(2)=DGASH(3)*EGASH(1)-DGASH(1)*EGASH(3)
              EXTNOR(3)=DGASH(1)*EGASH(2)-DGASH(2)*EGASH(1)
C DGASH and EGASH are not unit vectors
              RDGASH=NORM2(DGASH)
              REGASH=NORM2(EGASH)
C
              PXCOMP=-PGASH(1)*EXTNOR(1)+PGASH(2)*DGASH(1)*REGASH+
     1                PGASH(3)*EGASH(1)*RDGASH
              PYCOMP=-PGASH(1)*EXTNOR(2)+PGASH(2)*DGASH(2)*REGASH+
     1                PGASH(3)*EGASH(2)*RDGASH
              PZCOMP=-PGASH(1)*EXTNOR(3)+PGASH(2)*DGASH(3)*REGASH+
     1                PGASH(3)*EGASH(3)*RDGASH
            ELSE
              PXCOMP=DGASH(1)*PGASH(2)-DGASH(2)*PGASH(1)
              PYCOMP=DGASH(1)*PGASH(1)+DGASH(2)*PGASH(2)
            ENDIF
C
            DVOLU=WEIGPB
            IF(NTYPE.EQ.1)THEN
C interpolate to find thickness at boundary gauss point (plane stress)
              THICK=R0
              DO 115 INODEG=1,NNODEG
                II=NODAUX(NOPRS(INODEG))
                INODE=NODCHK(II)
                THICK=THICK+THKN(INODE)*SHAPE(II)
  115         CONTINUE
              DVOLU=DVOLU*THICK
            ELSEIF(NTYPE.EQ.3)THEN
C interpolate to find radius at boundary gauss point (axisymmetric case)
              RADUS=R0
              DO 117 INODEG=1,NNODEG
                II=NODAUX(NOPRS(INODEG))
                RADUS=RADUS+SHAPE(II)*ELCOD(NAXIS,INODEG)
  117         CONTINUE
              DVOLU=DVOLU*TWOPI*RADUS
            ENDIF
C
C Add the equivalent nodal loads (consistent loads) to the element force
C vector:
C contribution of current boundary Gauss point to integral of:
C (boundary shape functions) * (nodal prescribed pressures)
            DO 130 INODEG=1,NNODEG
              INODE=NODCHK(INODEG)
              NGASH=(INODE-1)*NDOFN+1
              MGASH=(INODE-1)*NDOFN+2
              RLOAD(NGASH,IELEM)=RLOAD(NGASH,IELEM)+
     1                           SHAPE(INODEG)*PXCOMP*DVOLU
              RLOAD(MGASH,IELEM)=RLOAD(MGASH,IELEM)+ 
     1                           SHAPE(INODEG)*PYCOMP*DVOLU
C Z component needed only in 3-D
              IF(NTYPE.EQ.4)THEN
                LGASH=(INODE-1)*NDOFN+3
                RLOAD(LGASH,IELEM)=RLOAD(LGASH,IELEM)+ 
     1                             SHAPE(INODEG)*PZCOMP*DVOLU
              ENDIF
  130       CONTINUE
  150     CONTINUE
  160   CONTINUE
      ENDIF
C
C Issue error message if I/O error occurs while reading input data file
      GOTO 940
  900 CALL ERRPRT('ED0214')
  910 CALL ERRPRT('ED0215')
  920 CALL ERRPRT('ED0216')
  930 CALL ERRPRT('ED0217')
C
  940 CONTINUE
      RETURN
      END
CDOC END_SUBROUTINE INLOAD
