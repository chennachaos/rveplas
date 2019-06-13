CDOC BEGIN_SUBROUTINE IFSTD
CDOC Internal force vector computation for 2-D and 3-D elements of class
CDOC STDARD
CDOC
CDOC This routine computes the element internal force vector for all
CDOC 2-D (for plane strain, plane stress and axisymmetric analysis) 
CDOC and 3-D elements of the class STDARD: standard isoparametric 
CDOC elements.
CDOC
CDOC BEGIN_PARAMETERS
CDOC INTEGER          IELEM  >  Number of the element whose internal
CDOC C                          force vector is to be computed.
CDOC LOGICAL          INCCUT <  Increment cutting flag. Return value set
CDOC C                          to .FALSE. if internal force vector
CDOC C                          was successfully evaluated and set to
CDOC C                          .TRUE. otherwise. When the return value
CDOC C                          is set to .TRUE., the main program will
CDOC C                          activate increment cutting and the
CDOC C                          current increment will be divided into
CDOC C                          sub-increments.
CDOC INTEGER          MDIME  >  Maximum permissible number of spatial
CDOC C                          dimensions.
CDOC INTEGER          MELEM  >  Maximum permissible number of elements
CDOC C                          in a mesh.
CDOC INTEGER          MPOIN  >  Maximum permissible number of nodal
CDOC C                          points in a mesh.
CDOC INTEGER          MSTRE  >  Maximum permissible number of stress
CDOC C                          (or resultant force) components.
CDOC INTEGER          MTOTV  >  Maximum permissible degrees of freedom
CDOC C                          in a mesh.
CDOC INTEGER          NAXIS  >  Symmetry axis flag. Used only for
CDOC C                          axisymmetric analysis. NAXIS=1 if
CDOC C                          symmetric about the X-axis and
CDOC C                          NAXIS=2 if symmetric about Y-axis.
CDOC INTEGER          NLARGE >  Large deformation flag. Large
CDOC C                          deformation analysis if
CDOC C                          NLARGE=1 and infinitesimal
CDOC C                          deformation analysis otherwise.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC DOUBLE_PRECISION COORD1 >  Array of coordinates of all nodes of
CDOC C                          the mesh in the current configuration.
CDOC DOUBLE_PRECISION DINCR  >  Global array of incremental
CDOC C                          displacements of the entire mesh.
CDOC DOUBLE_PRECISION DTIME  >  Time increment.
CDOC DOUBLE_PRECISION ELOAD  <  Internal force vector of the current
CDOC C                          element.
CDOC INTEGER          IELPRP >  Array of integer element properties
CDOC C                          of the current element group.
CDOC INTEGER          IPROPS >  Array of integer material properties
CDOC C                          of the current element group.
CDOC LOGICAL          LALGVA <  Array of current logical algorithmic
CDOC C                          variables of all Gauss points of the
CDOC C                          current element.
CDOC INTEGER          LNODS  >  Global array of connectivity of the
CDOC C                          entire mesh.
CDOC DOUBLE_PRECISION RALGVA <  Array of current real algorithmic
CDOC C                          variables of all Gauss points of the
CDOC C                          current element.
CDOC DOUBLE_PRECISION RELPRP >  Array of real element properties of
CDOC C                          the current element group.
CDOC DOUBLE_PRECISION RPROPS >  Array of real material properties of
CDOC C                          the current element group.
CDOC DOUBLE_PRECISION RSTAVA <  Array of current real state
CDOC C                          variables of all Gauss points of the
CDOC C                          current element.
CDOC DOUBLE_PRECISION STRSG  <  Array of resultant forces of all
CDOC C                          centroids line Gauss points of the
CDOC C                          current element.
CDOC DOUBLE_PRECISION THKGP  <  Array of current thickness of
CDOC C                          all Gauss points of the current element.
CDOC C                          Used only in plane stress analysis.
CDOC C                          Updated here only in large strain
CDOC C                          analysis.
CDOC DOUBLE_PRECISION TDISP  >  Global array of total displacement
CDOC C                          of all nodes of the mesh.
CDOC INTEGER          NDIME  >  Number of spatial dimensions.
CDOC INTEGER          NDOFN  >  Number of degrees of freedom per node.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, June      1996: Initial coding (as INTFOR)
CHST
CHST E.de Souza Neto, February  1998: Change of program structure
CHST
CHST E.de Souza Neto, May       1998: Common blocks eliminated
CHST
CHST E.de Souza Neto, September 2002: F-bar element removed
CHST
CHST E.de Souza Neto & F.Adziman
CHST                , September 2012: Time increment info added.
CHST
CHST      F. Adziman, September 2013: Passing DVOLU into MATISU
CHST                                  for calc. of volume fraction
CHST                                  in SUMEPC 
CHST
CHST      F. Adziman, October   2013: Passing IELEM and IGAUSP into
CHST                                  MATISU for multivariant calc.
CHST                                  in SUMEPC      
CHST   D. de Bortoli, March     2015: Added 3-D case
CHST                                  Changed name from IFSTD2 to IFSTD
CHST
      SUBROUTINE IFSTD
     1(   IELEM      ,INCCUT     ,MDIME      ,MELEM      ,MPOIN      ,
     2    MSTRE      ,MTOTV      ,NAXIS      ,NLARGE     ,NTYPE      ,
     3    COORD1     ,DINCR      ,DTIME      ,ELOAD      ,IELPRP     ,
     4    IPROPS     ,LALGVA     ,LNODS      ,RALGVA     ,RELPRP     ,
     5    RPROPS     ,RSTAVA     ,STRSG      ,THKGP      ,TDISP      ,
     6    TTIME      ,NDIME      ,NDOFN      ,DVOLU      ,STREPG     )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE '../ELEMENTS.INC'
      INCLUDE '../MATERIAL.INC'
C
      PARAMETER( MGDIM=9 ,MBDIM=6 )
C Arguments
      LOGICAL  INCCUT ,LALGVA, OUTOFP
      DIMENSION
     1    COORD1(MDIME,MPOIN),DINCR(MTOTV)       ,ELOAD(MEVAB)       ,
     2    IELPRP(*)          ,IPROPS(*)         ,LALGVA(MLALGV,MTOTG),
     3    LNODS(MELEM,MEVAB),RALGVA(MRALGV,MTOTG),RELPRP(*)          ,
     4    RPROPS(*)          ,RSTAVA(MRSTAV,MTOTG),STRSG(MSTRE,MTOTG),
     5    THKGP(MTOTG)       ,TDISP(MTOTV)        ,STREPG(MSTRE,MTOTG)
C Local variables and arrays
      LOGICAL  SUFAIL
      DIMENSION
     1    BMATX(MBDIM,MEVAB) ,CARTD(NDIME,MNODE) ,DELDIS(MDOFN,MNODE),
     2    DERIV(NDIME,MNODE) ,EINCR(MBDIM)       ,ELCOD(NDIME,MNODE) ,
     3    FINCIN(3,3)        ,FINCR(3,3)         ,FINV(3,3)          ,
     4    GMATX(MGDIM,MEVAB) ,GPCOD(NDIME)       ,SHAPE(MNODE)       ,
     5    TELDIS(MDOFN,MNODE),FTOTA(3,3)         ,EISCRD(NDIME)
C Local numerical constants
      DATA
     1    R0   ,R1   ,R8   /
     2    0.0D0,1.0D0,8.0D0/
C***********************************************************************
C COMPUTE INTERNAL FORCE VECTOR OF ALL ELEMENTS OF CLASS 'STDARD'
C (STANDARD ISOPARAMETRIC ELEMENTS) IN 2-D (PLANE STRAIN, PLANE STRESS
C AND AXISYMMETRIC PROBLEMS) AND 3-D
C
C REFERENCE: Section 4.1.2
C            Box 4.2, item (viii) 
C***********************************************************************
      IF(NTYPE.EQ.1)THEN
        NBDIM=3
        OUTOFP=.FALSE.
      ELSEIF(NTYPE.EQ.2)THEN
        NBDIM=3
        OUTOFP=.FALSE.
      ELSEIF(NTYPE.EQ.3)THEN
        TWOPI=R8*ATAN(R1)
        NBDIM=4
        OUTOFP=.TRUE.
      ELSEIF(NTYPE.EQ.4)THEN
        NBDIM=6
        OUTOFP=.FALSE.
      ELSE
        CALL ERRPRT('EI0012')
      ENDIF
C Identify element type
      IELTYP=IELPRP(1)
C retrieve some element integer properties
      NNODE =IELPRP(3)
      NGAUSP=IELPRP(4)
      NEVAB =IELPRP(5)
C
C Set element arrays of current nodal coordinates, total and incremental
C displacements
      DO 20 INODE =1,NNODE
        LNODE=IABS(LNODS(IELEM,INODE))
        NPOSN=(LNODE-1)*NDOFN
        DO 10 IDOFN=1,NDOFN
          NPOSN=NPOSN+1
          IF(NLARGE.EQ.1)THEN
            ELCOD(IDOFN,INODE)=COORD1(IDOFN,LNODE)
            TELDIS(IDOFN,INODE)=-TDISP(NPOSN)
            DELDIS(IDOFN,INODE)=-DINCR(NPOSN)
          ELSE
            ELCOD(IDOFN,INODE)=COORD1(IDOFN,LNODE)
            DELDIS(IDOFN,INODE)=DINCR(NPOSN)
          ENDIF
   10   CONTINUE
   20 CONTINUE
C
C Initialise element force vector
C
      CALL RVZERO(ELOAD,NEVAB)
C
C=======================================================================
C                 Begin loop over Gauss points                         |
C=======================================================================
C
      IPWEI=NGAUSP*NDIME
      DO 40 IGAUSP=1,NGAUSP
        IPPOS=NDIME*(IGAUSP-1)
C Set Gauss points positions and weights
        DO IDIME=1,NDIME
          EISCRD(IDIME)=RELPRP(IPPOS+IDIME)
        ENDDO
        WEIGP=RELPRP(IPWEI+IGAUSP)
C Evaluate shape functions and their derivatives (use current
C configuration for large strains)
        CALL SHPFUN
     1(   DERIV      ,EISCRD      ,0          ,IELTYP     ,
     2    NDIME      ,SHAPE      )
        CALL JACDET
     1(   CARTD      ,DERIV      ,DETJAC     ,ELCOD      ,IELEM      ,
     2    NDIME      ,NDIME      ,NNODE      )
        IF(DETJAC.LE.R0)THEN
C ...cut increment if current jacobian is not positive definite
          CALL ERRPRT('WE0009')
          INCCUT=.TRUE.
          GOTO 999
        ENDIF
        IF(NTYPE.EQ.3)CALL GETGCO
     1(   GPCOD      ,ELCOD      ,NDIME      ,NDIME      ,NNODE      ,
     2    SHAPE      )
C Evaluate symmetric gradient operator B
        CALL GETBMX
     1(   BMATX      ,GPCOD      ,CARTD      ,NDIME      ,MBDIM      ,
     2    NAXIS      ,NNODE      ,NTYPE      ,SHAPE      )
C
C Compute basic kinematic variables needed for state update
C =========================================================
C
        IF(NLARGE.EQ.1)THEN
C Large strains: compute incremental deformation gradient
C -------------------------------------------------------
C gradient operator G in current configuration
          CALL GETGMX
     1(   GPCOD      ,CARTD      ,GMATX      ,NDIME      ,MGDIM      ,
     2    NAXIS      ,NNODE      ,NTYPE      ,SHAPE      )
C incremental deformation gradient
          CALL DEFGRA
     1(   DELDIS     ,FINCIN     ,GMATX      ,MDOFN      ,MGDIM      ,
     2    NDOFN      ,NTYPE      ,NNODE      )
          IF(NTYPE.EQ.3.AND.FINCIN(3,3).LE.R0)THEN
C ...cut increment if determinant of incr. def. gradient non positive
            CALL ERRPRT('WE0010')
            INCCUT=.TRUE.
            GOTO 999
          ENDIF
          CALL RVZERO(FINCR,3*3) 
          CALL INVF
     1(   FINCIN     ,FINCR      ,NTYPE      )
C... compute determinant of total deformation gradient
          CALL DEFGRA
     1(   TELDIS     ,FINV       ,GMATX      ,MDOFN      ,MGDIM      ,
     2    NDOFN      ,NTYPE      ,NNODE      )     
C
          DETFIN=DETM23(3, FINV, NDIME, OUTOFP)
C stops program if deformation gradient is not positive definite
          IF((NTYPE.EQ.3).AND.(FINV(3,3).LE.R0))THEN
              CALL ERRPRT('WE0010')
              INCCUT=.TRUE.
              GOTO 999
          ENDIF
C
          DETF=R1/DETFIN
C... compute total deformation gradient
          CALL RVZERO(FTOTA,3*3)
          CALL INVF
     1(   FINV     ,FTOTA      ,NTYPE      )
        ELSE		
C Small strains: compute incremental infinitesimal strain
C -------------------------------------------------------
          CALL LISTRA
     1(   BMATX      ,DELDIS     ,MDOFN      ,MBDIM      ,NDOFN      ,
     2    NNODE      ,NTYPE      ,EINCR      )
        ENDIF
C
C Evaluate elemental volume for plane strain (ntype=2)
C ----------------------------------------------------
        DVOLU=DETJAC*WEIGP 
C        
C Call material interface routine for state update calls: Update stress
C and other state variables
C =====================================================================
C
        CALL MATISU
     1(   DETF       ,NLARGE     ,NTYPE      ,SUFAIL     ,
     2    THKGP(IGAUSP)          ,DTIME      ,EINCR      ,FINCR      ,
     3    IPROPS     ,LALGVA(1,IGAUSP)       ,RALGVA(1,IGAUSP)       ,
     4    RPROPS     ,RSTAVA(1,IGAUSP)       ,STRSG(1,IGAUSP)        ,
     5    TTIME      ,FTOTA      ,DVOLU      ,IELEM      ,IGAUSP     ,    
     6    STREPG(1,IGAUSP)       )
        IF(SUFAIL)THEN
C State updating failed for current Gauss point: break loop over Gauss
C points and exit with increment cutting flag activated
          INCCUT=.TRUE.
          GOTO 999
        ENDIF
C
C Add current Gauss point contribution to the element internal force
C vector
C ==================================================================
C
C evaluate elemental volume
        IF(NTYPE.EQ.1)THEN
          DVOLU=DVOLU*THKGP(IGAUSP)
        ELSEIF(NTYPE.EQ.3)THEN
          DVOLU=DVOLU*TWOPI*GPCOD(NAXIS)
        ENDIF
C                           T
C add current gauss point  B [sigma]
        CALL RTV
     1(   0          ,MBDIM      ,NEVAB      ,NBDIM      ,ELOAD      ,
     2    BMATX      ,STRSG(1,IGAUSP)        ,DVOLU      )
C
   40 CONTINUE
C
C=======================================================================
C                 End of loop over Gauss points                        |
C=======================================================================
C
  999 CONTINUE
      RETURN
      END
CDOC END_SUBROUTINE IFSTD
