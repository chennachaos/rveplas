CDOC BEGIN_SUBROUTINE IFFBA
CDOC Internal force vector computation for 2-D and 3-D F-bar elements
CDOC
CDOC This routine computes the element internal force vector for all
CDOC 2-D and 3-D elements of the class FBAR: F-bar elements for plane 
CDOC strain, axisymmetric and three-dimensional analysis under large 
CDOC deformations.
CDOC
CDOC BEGIN_PARAMETERS
CDOC INTEGER          IELEM  >  Number of the element whose internal
CDOC C                          force vector is to be computed.
CDOC LOGICAL          INCCUT <  Increment cutting flag. Return value set
CDOC C                          to .FALSE. if internal force vector was
CDOC C                          successfully evaluated and set to .TRUE.
CDOC C                          otherwise. When the return value is set
CDOC C                          to .TRUE., the main program will
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
CDOC C                          NAXIS=2 if symmetric about the Y-axis.
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
CHST E.de Souza Neto, September 2002: Initial coding
CHST
CHST E.de Souza Neto & F.Adziman,
CHST                  September 2012: Time increment added
CHST
CHST F. Adziman,      September 2013: Passing DVOLU into MATISU
CHST                                  for calc. of volume fraction
CHST                                  in SUMEPC 
CHST
CHST F. Adziman,      October 2013  : Passing IELEM and IGAUSP into
CHST                                  MATISU for multivariant calc.
CHST                                  in SUMEPC 
CHST   D. de Bortoli, March     2015: Added 3-D case
CHST                                  Changed name from IFFBA2 to IFFBA
CHST
      SUBROUTINE IFFBA
     1(   IELEM      ,INCCUT     ,MDIME      ,MELEM      ,MPOIN      ,
     2    MSTRE      ,MTOTV      ,NAXIS      ,NTYPE      ,
     3    COORD1     ,DINCR      ,DTIME      ,ELOAD      ,IELPRP     ,
     4    IPROPS     ,LALGVA     ,LNODS      ,RALGVA     ,RELPRP     ,
     5    RPROPS     ,RSTAVA     ,STRSG      ,THKGP      ,TDISP      ,
     6    TTIME      ,NDIME      ,NDOFN      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE '../ELEMENTS.INC'
      INCLUDE '../MATERIAL.INC'
C
      PARAMETER( MGDIM=9 ,MBDIM=6 )
C Arguments
      LOGICAL  INCCUT, LALGVA, OUTOFP
      DIMENSION
     1    COORD1(MDIME,MPOIN),DINCR(MTOTV)       ,ELOAD(MEVAB)       ,
     2    IELPRP(*)          ,IPROPS(*)         ,LALGVA(MLALGV,MTOTG),
     3    LNODS(MELEM,MEVAB),RALGVA(MRALGV,MTOTG),RELPRP(*)          ,
     4    RPROPS(*)          ,RSTAVA(MRSTAV,MTOTG),STRSG(MSTRE,MTOTG),
     5    THKGP(MTOTG)       ,TDISP(MTOTV)
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
     1    R0   ,RP5  ,R1   ,R3   ,R8   /
     2    0.0D0,0.5D0,1.0D0,3.0D0,8.0D0/
C***********************************************************************
C COMPUTE INTERNAL FORCE VECTOR OF ALL ELEMENTS OF CLASS 'FBAR'
C (F-Bar ELEMENTS) IN 2-D (PLANE STRAIN AND AXISYMMETRIC) AND 3-D
C
C REFERENCE: Box 15.1
C***********************************************************************
      R1D3=R1/R3
      IF(NTYPE.EQ.2)THEN
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
C F-bar implementation valid only for plane strain, axisymmetric and
C three-dimensional cases
        CALL ERRPRT('EI0033')
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
          ELCOD(IDOFN,INODE)=COORD1(IDOFN,LNODE)
          TELDIS(IDOFN,INODE)=-TDISP(NPOSN)
          DELDIS(IDOFN,INODE)=-DINCR(NPOSN)
   10   CONTINUE
   20 CONTINUE
C
C Initialise element force vector
C
      CALL RVZERO(ELOAD,NEVAB)
C
C-----------------------------------------------------------------------
C Calculation of the F-bar deformation gradient determinat
C-----------------------------------------------------------------------
C Evaluate inverse of the incremental deformation gradient at the
C centroid of the F-bar element
      NGAUSB=IELPRP(8)
C Number of spatial dimensions of element boundary (only used to locate
C centroid isoparametric coordinates in array RELPRP)
      IF(NTYPE.EQ.4)THEN
        NDIMEB=2
      ELSE
        NDIMEB=1
      ENDIF
      IPOS=NGAUSP*NDIME+NGAUSP+NGAUSP*NNODE+NGAUSB*NDIMEB+NGAUSB+1
      DO IDIME=1,NDIME
        EISCRD(IDIME)=RELPRP(IPOS)
        IPOS=IPOS+1
      ENDDO
C
      CALL SHPFUN
     1(   DERIV      ,EISCRD      ,0          ,IELTYP     ,
     2    NDIME      ,SHAPE      )
      CALL JACDET
     1(   CARTD      ,DERIV      ,DETJAC     ,ELCOD      ,IELEM      ,
     2    NDIME      ,NDIME      ,NNODE      )
      IF(DETJAC.LE.R0)THEN
C cut increment if element jacobian is not positive definite
        CALL ERRPRT('WE0021')
        INCCUT=.TRUE.
        GOTO 999
      ENDIF
      IF(NTYPE.EQ.3)CALL GETGCO
     1(   GPCOD      ,ELCOD      ,NDIME      ,NDIME      ,NNODE      ,
     2    SHAPE      )
C          
      CALL GETGMX
     1(   GPCOD      ,CARTD      ,GMATX      ,NDIME      ,MGDIM      ,
     2    NAXIS      ,NNODE      ,NTYPE      ,SHAPE      )
C Determinant of the incremental deformation gradient at the centroid
      CALL DEFGRA
     1(   DELDIS     ,FINCIN     ,GMATX      ,MDOFN      ,MGDIM      ,
     2    NDOFN      ,NTYPE      ,NNODE      )
C     
      IF(NTYPE.EQ.2)THEN
        AFACT=RP5
      ELSEIF((NTYPE.EQ.3).OR.(NTYPE.EQ.4))THEN
        AFACT=R1D3
      ENDIF
C... cut increment if incr. def. gradient is not positive definite
      IF((NTYPE.EQ.3).AND.(FINCIN(3,3).LE.R0))THEN
        CALL ERRPRT('WE0022')
        INCCUT=.TRUE.
        GOTO 999
      ENDIF
C
      DET=DETM23(3, FINCIN, NDIME, OUTOFP)
C      
      DET0=R1/DET
C Determinant of the total deformation gradient at the centroid
      CALL DEFGRA
     1(   TELDIS     ,FINV       ,GMATX      ,MDOFN      ,MGDIM      ,
     2    NDOFN      ,NTYPE      ,NNODE      )
C... cut increment if deformation gradient is not positive definite
      IF((NTYPE.EQ.3).AND.(FINV(3,3).LE.R0))THEN
        CALL ERRPRT('WE0023')
        INCCUT=.TRUE.
        GOTO 999
      ENDIF
C
      DET=DETM23(3, FINV, NDIME, OUTOFP)
C      
      DETF0=R1/DET
C-----------------------------------------------------------------------
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
C        
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
          CALL ERRPRT('WE0024')
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
          CALL ERRPRT('WE0025')
          INCCUT=.TRUE.
          GOTO 999
        ENDIF
        CALL RVZERO(FINCR,3*3)
        CALL INVF
     1(   FINCIN     ,FINCR      ,NTYPE      )
C Modified incremental deformation gradient for F-bar element
        DET=DETM23(3, FINCR, NDIME, OUTOFP)
C
        FACTOR=(DET0/DET)**AFACT
        FINCR(1:NDIME,1:NDIME)=FACTOR*FINCR(1:NDIME,1:NDIME)
        IF(NTYPE.EQ.3)FINCR(3,3)=FACTOR*FINCR(3,3)
        
        
        
C... determinant of total deformation gradient
        DETF=DETF0
C... compute total deformation gradient
        CALL RVZERO(FTOTA,3*3)
        CALL INVF
     1(   FINV     ,FTOTA      ,NTYPE      )
C
C Evaluate elemental volume for plane strain (ntype=2)
C ----------------------------------------------------
        DVOLU=DETJAC*WEIGP
C
C Call material interface routine for state update calls: Update stress
C and other state variables
C =====================================================================
C
        NLARGE=1
        CALL MATISU
     1(   DETF       ,NLARGE     ,NTYPE      ,SUFAIL     ,
     2    THKGP(IGAUSP)          ,DTIME      ,EINCR      ,FINCR      ,
     3    IPROPS     ,LALGVA(1,IGAUSP)       ,RALGVA(1,IGAUSP)       ,
     4    RPROPS     ,RSTAVA(1,IGAUSP)       ,STRSG(1,IGAUSP)        ,
     5    TTIME      ,FTOTA      ,DVOLU      ,IELEM      ,IGAUSP     )
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
        IF(NTYPE.EQ.3)THEN
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
CDOC END_SUBROUTINE IFFBA
