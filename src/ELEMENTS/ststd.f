CDOC BEGIN_SUBROUTINE STSTD
CDOC Stiffness matrix computation for elements of class STDARD
CDOC
CDOC This routine computes the tangent stiffness matrix for all
CDOC elements of the class STDARD: standard isoparametric 3-D elements
CDOC and 2-D elements (plane strain, plane stress and axisymmetric 
CDOC analysis).
CDOC
CDOC BEGIN_PARAMETERS
CDOC INTEGER          IELEM  >  Number of the element whose stiffness
CDOC C                          matrix is to be computed.
CDOC INTEGER          KUNLD  >  Unloading flag. KUNLD is set to 1 if the
CDOC C                          loading programme is currently
CDOC C                          unloading.
CDOC INTEGER          MDIME  >  Maximum permissible number of spatial
CDOC C                          dimensions.
CDOC INTEGER          MELEM  >  Maximum permissible number of elements.
CDOC INTEGER          MPOIN  >  Maximum permissible number of nodal
CDOC C                          points in a mesh.
CDOC INTEGER          MSTRE  >  Maximum permissible number of stress
CDOC C                          (or resultant force) components.
CDOC INTEGER          MTOTV  >  Maximum permissible number degrees of
CDOC C                          freedom in a mesh.
CDOC INTEGER          NAXIS  >  Symmetry axis flag. Used only for
CDOC C                          axisymmetric analysis. NAXIS=1 if
CDOC C                          symmetric about the X-axis and
CDOC C                          NAXIS=2 if symmetric about the Y-axis.
CDOC INTEGER          NLARGE >  Large deformation flag. Large
CDOC C                          deformation analysis if NLARGE=1 and
CDOC C                          infinitesimal deformation analysis
CDOC C                          otherwise.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC LOGICAL          UNSYM  >  Tangent stiffness symmetry flag.
CDOC C                          .TRUE. if tangent stiffness is
CDOC C                          unsymmetric, .FALSE. otherwise.
CDOC DOUBLE_PRECISION COORD1 >  Array of coordinates of all nodes of
CDOC C                          the mesh in the current configuration.
CDOC DOUBLE_PRECISION DINCR  >  Global array of incremental
CDOC C                          displacements of the entire mesh.
CDOC DOUBLE_PRECISION DTIME  >  Time increment.
CDOC DOUBLE_PRECISION ESTIF  <  Tangent stiffness matrix of the current
CDOC C                          element.
CDOC INTEGER          IELPRP >  Array of integer element properties
CDOC C                          of the current element group.
CDOC INTEGER          IPROPS >  Array of integer material properties
CDOC C                          of the current element group.
CDOC LOGICAL          LALGVA >  Array of current logical algorithmic
CDOC C                          variables of all Gauss points of the
CDOC C                          current element.
CDOC INTEGER          LNODS  >  Global array of connectivity of the
CDOC C                          entire mesh.
CDOC DOUBLE_PRECISION RALGVA >  Array of current real algorithmic
CDOC C                          variables of all Gauss points of the
CDOC C                          current element.
CDOC DOUBLE_PRECISION RELPRP >  Array of real element properties of
CDOC C                          the current element group.
CDOC DOUBLE_PRECISION RPROPS >  Array of real material properties of
CDOC C                          the current element group.
CDOC DOUBLE_PRECISION RSTAVA >  Array of current real state
CDOC C                          variables of all Gauss points of the
CDOC C                          current element.
CDOC DOUBLE_PRECISION RSTAV2 >  Array of real state
CDOC C                          variables of all Gauss points of the
CDOC C                          current element at the previous
CDOC C                          equilibrium (converged) configuration.
CDOC DOUBLE_PRECISION STRSG  >  Array of current stress components of
CDOC C                          all Gauss points of the current element.
CDOC DOUBLE_PRECISION THKGP  >  Array of current thickness of
CDOC C                          all Gauss points of the current element.
CDOC C                          Used only in plane stress analysis.
CDOC DOUBLE_PRECISION TDISP  >  Global array of total displacement
CDOC C                          of all nodes of the mesh.
CDOC INTEGER          NDIME  >  Number of spatial dimensions.
CDOC INTEGER          NDOFN  >  Number of degrees of freedom per node.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, June      1996: Initial coding (as STIFF)
CHST
CHST E.de Souza Neto, February  1998: Change of program structure
CHST
CHST E.de Souza Neto, May       1998: Common blocks eliminated
CHST
CHST E.de Souza Neto & F.Adziman,
CHST                  September 2012: Time increment added
CHST
CHST   D. de Bortoli, March     2015: Added 3-D case
CHST                                  Changed name from STSTD2 to STSTD
CHST
      SUBROUTINE STSTD
     1(   IELEM      ,KUNLD      ,MDIME      ,MELEM      ,
     2    MPOIN      ,MSTRE      ,MTOTV      ,NAXIS      ,NLARGE     ,
     3    NTYPE      ,UNSYM      ,
     4    COORD1     ,DINCR      ,DTIME      ,ESTIF      ,IELPRP     ,
     5    IPROPS     ,LALGVA     ,LNODS      ,RALGVA     ,RELPRP     ,
     6    RPROPS     ,RSTAVA     ,RSTAV2     ,STRSG      ,THKGP      ,
     7    TDISP      ,NDIME      ,NDOFN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C Hyplas database
      INCLUDE '../MATERIAL.INC'
      INCLUDE '../ELEMENTS.INC'
C
      PARAMETER( MGDIM=9 ,MBDIM=6 )
C Arguments
      LOGICAL  LALGVA ,UNSYM ,OUTOFP
      DIMENSION
     1    COORD1(MDIME,MPOIN),DINCR(MTOTV)       ,ESTIF(MEVAB,MEVAB) ,
     2    IELPRP(*)          ,IPROPS(*)         ,LALGVA(MLALGV,MTOTG),
     3    LNODS(MELEM,MEVAB) ,RALGVA(MRALGV,MTOTG),RELPRP(*)         ,
     4    RPROPS(*)        ,RSTAVA(MRSTAV,MTOTG),RSTAV2(MRSTAV,MTOTG),
     5    STRSG(MSTRE,MTOTG) ,THKGP(MTOTG)       ,TDISP(MTOTV)
C Local arrays and variables
      DIMENSION
     1    AUXM(MEVAB,MGDIM)  ,AMATX(MGDIM,MGDIM) ,BMATX(MBDIM,MEVAB) ,
     2    CARTD(NDIME,MNODE) ,DELDIS(MDOFN,MNODE),DERIV(NDIME,MNODE) ,
     3    DMATX(MBDIM,MBDIM) ,EINCR(MBDIM)       ,ELCOD(NDIME,MNODE) ,
     4    FINCIN(3,3)        ,FINCR(3,3)         ,FINV(3,3)          ,
     5    GMATX(MGDIM,MEVAB) ,GPCOD(NDIME)       ,SHAPE(MNODE)       ,
     6    TELDIS(MDOFN,MNODE),EISCRD(NDIME)
      DATA R0   ,R1   ,R8   /
     1     0.0D0,1.0D0,8.0D0/
C***********************************************************************
C EVALUATES THE ELEMENT TANGENT STIFFNESS MATRIX FOR ELEMENTS OF CLASS
C 'STDARD' (STANDARD ISOPARAMETRIC ELEMENTS) IN 2-D (PLANE STRAIN,
C PLANE STRESS AND AXISYMMETRIC PROBLEMS) AND 3-D
C
C REFERENCE: Sections 4.2.4, 4.3.4
C            Box 4.2, item (iii)
C            Box 4.3
C***********************************************************************
      IF(NTYPE.EQ.1)THEN
          NBDIM=3
          NGDIM=4
          OUTOFP=.FALSE.
      ELSEIF(NTYPE.EQ.2)THEN
          NBDIM=3
          NGDIM=4
          OUTOFP=.FALSE.
      ELSEIF(NTYPE.EQ.3)THEN
          NBDIM=4
          NGDIM=5
          OUTOFP=.TRUE.
          TWOPI=R8*ATAN(R1)
      ELSEIF(NTYPE.EQ.4)THEN
          NBDIM=6
          NGDIM=9
          OUTOFP=.FALSE.
      ELSE
        CALL ERRPRT('EI0079')
      ENDIF
C Identify element type
      IELTYP=IELPRP(1)
C Recover element information
      NNODE =IELPRP(3)
      NGAUSP=IELPRP(4)
      NEVAB =IELPRP(5)
C
C Set element nodal coordinates, total and incremental displacements
C vectors
C
      DO 20 INODE =1,NNODE
        LNODE=IABS(LNODS(IELEM,INODE))
        NPOSN=(LNODE-1)*NDOFN
        DO 10 IDOFN=1,NDOFN
          NPOSN=NPOSN+1
          ELCOD(IDOFN,INODE)=COORD1(IDOFN,LNODE)
          IF(NLARGE.EQ.1)THEN
            TELDIS(IDOFN,INODE)=-TDISP(NPOSN)
            DELDIS(IDOFN,INODE)=-DINCR(NPOSN)
          ELSE
            DELDIS(IDOFN,INODE)=DINCR(NPOSN)
          ENDIF
   10   CONTINUE
   20 CONTINUE
C
C Initialize the element stiffness matrix
C
      DO 40 IEVAB=1,NEVAB
        DO 30 JEVAB=1,NEVAB
          ESTIF(IEVAB,JEVAB)=R0
   30   CONTINUE
   40 CONTINUE
C
C=======================================================================
C                 Begin loop over Gauss points                         |
C=======================================================================
      IPWEI=NGAUSP*NDIME
      DO 70 IGAUSP=1,NGAUSP
        IPPOS=NDIME*(IGAUSP-1)
        DO IDIME=1,NDIME
          EISCRD(IDIME)=RELPRP(IPPOS+IDIME)
        ENDDO
        WEIGP=RELPRP(IPWEI+IGAUSP)
C Evaluate the shape functions and derivatives
        CALL SHPFUN
     1(   DERIV      ,EISCRD      ,0          ,IELTYP     ,
     2    NDIME      ,SHAPE      )
        CALL JACDET
     1(   CARTD      ,DERIV      ,DETJAC     ,ELCOD      ,IELEM      ,
     2    NDIME      ,NDIME      ,NNODE      )
        IF(DETJAC.LE.R0)THEN
C... stops program if element jacobian is not positive definite
          CALL ERRPRT('EE0003')
        ENDIF
        IF(NTYPE.EQ.3)CALL GETGCO
     1(   GPCOD      ,ELCOD      ,NDIME      ,NDIME      ,NNODE      ,
     2    SHAPE      )
C
        IF(NLARGE.EQ.1)THEN
C Large strains: compute incremental deformation gradient
C -------------------------------------------------------
C gradient operator G in the current configuration
          CALL GETGMX
     1(   GPCOD      ,CARTD      ,GMATX      ,NDIME      ,MGDIM      ,
     2    NAXIS      ,NNODE      ,NTYPE      ,SHAPE      )
C inverse of incremental deformation gradient
          CALL DEFGRA
     1(   DELDIS     ,FINCIN     ,GMATX      ,MDOFN      ,MGDIM      ,
     2    NDOFN      ,NTYPE      ,NNODE      )
C stops program if deformation gradient determinant is non-positive
          IF(NTYPE.EQ.3.AND.FINCIN(3,3).LE.R0)CALL ERRPRT('EE0004')
C incremental deformation gradient
          CALL INVF
     1(   FINCIN     ,FINCR      ,NTYPE      )
C... and the determinant of the total deformation gradient
          CALL DEFGRA
     1(   TELDIS     ,FINV       ,GMATX      ,MDOFN      ,MGDIM      ,
     2    NDOFN      ,NTYPE      ,NNODE      )
C
          DETFIN=DETM23(3, FINV, NDIME, OUTOFP)
C stops program if deformation gradient is not positive definite
          IF((NTYPE.EQ.3).AND.(FINV(3,3).LE.R0))CALL ERRPRT('EE0004')
C
          DETF=R1/DETFIN
        ELSE
C Small strains: compute incremental infinitesimal strain
C -------------------------------------------------------
C compute the symmetric gradient operator B
          CALL GETBMX
     1(   BMATX      ,GPCOD      ,CARTD      ,NDIME      ,MBDIM      ,
     2    NAXIS      ,NNODE      ,NTYPE      ,SHAPE      )
C and the incremental infinitesimal strain
          CALL LISTRA
     1(   BMATX      ,DELDIS     ,MDOFN      ,MBDIM      ,NDOFN      ,
     2    NNODE      ,NTYPE      ,EINCR      )
        ENDIF
C
C Call material interface routine for consistent tangent computation
C calls: Compute either the standard consistent tangent matrix DMATX
C (small strains) or the spatial tangent modulus AMATX (large strains)
C ====================================================================
C
        CALL MATICT
     1(   DETF       ,KUNLD      ,MBDIM      ,MGDIM      ,MSTRE      ,
     2    NLARGE     ,NTYPE      ,
     3    AMATX      ,DMATX      ,DTIME      ,EINCR      ,FINCR      ,
     4    IPROPS     ,LALGVA(1,IGAUSP)       ,RALGVA(1,IGAUSP)       ,
     5    RPROPS     ,RSTAVA(1,IGAUSP)       ,RSTAV2(1,IGAUSP)       ,
     6    STRSG(1,IGAUSP)        ,FINV       )      
C
C Add current Gauss point contribution to element stiffness
C =========================================================
C
C Compute elemental volume
        DVOLU=DETJAC*WEIGP
        IF(NTYPE.EQ.1)THEN
          DVOLU=DVOLU*THKGP(IGAUSP)
        ELSEIF(NTYPE.EQ.3)THEN
          DVOLU=DVOLU*TWOPI*GPCOD(NAXIS)
        ENDIF
        IF(NLARGE.EQ.1)THEN
C                                                         T
C Large strains:  assemble the element stiffness as  K:= G [a] G
C --------------------------------------------------------------
          CALL RTSR
     1(   AUXM       ,0          ,MEVAB      ,MGDIM      ,NEVAB      ,
     2    NGDIM      ,ESTIF      ,GMATX      ,AMATX      ,DVOLU      ,
     3    UNSYM      )
        ELSE
C                                                        T
C Small strains: assemble the element stiffness as  K:= B D B
C -----------------------------------------------------------
          CALL RTSR
     1(   AUXM       ,0          ,MEVAB      ,MBDIM      ,NEVAB      ,
     2    NBDIM      ,ESTIF      ,BMATX      ,DMATX      ,DVOLU      ,
     3    UNSYM      )
        ENDIF
C
   70 CONTINUE
C=======================================================================
C                 End of loop over Gauss points                        |
C=======================================================================
      RETURN
      END
CDOC END_SUBROUTINE STSTD
