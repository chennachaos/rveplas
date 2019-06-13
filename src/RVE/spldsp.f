CDOC BEGIN_SUBROUTINE SPLDSP
CDOC Split into uniform displacements to accomodate different RVE
CDOC kinematical constraint choices.
CDOC
CDOC BEGIN_PARAMETERS
CDOC INTEGER          NINCS  >  Number of total load increment.
CDOC DOUBLE_PRECISION DUISTAR<  Incremental uniform strain displacement
CDOC                            for 'interior' nodes.
CDOC DOUBLE_PRECISION DUPSTAR<  Incremental uniform strain displacement
CDOC                            for 'plus' nodes.
CDOC DOUBLE_PRECISION DUMSTAR<  Incremental uniform strain displacement
CDOC                            for 'minus' nodes.
CDOC DOUBLE_PRECISION DUCSTAR<  Incremental uniform strain displacement
CDOC                            for 'corner' nodes.
CDOC END_PARAMETERS
CHST
CHST M.F. Adziman    ,  July 2013: Initial coding by referring to
CHST                               MICROPLAST (MATLAB)
CHST
      SUBROUTINE SPLDSP 
     1(    DUTYLR     ,NINCS    )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C Hyplas global database
      INCLUDE '../MAXDIM.INC'
      INCLUDE '../MATERIAL.INC'
      INCLUDE '../ELEMENTS.INC'
      INCLUDE '../GLBDBASE.INC'
      INCLUDE '../RVE.INC'
C Local arrays and variables
C      DIMENSION 
C     1   UISTAR(MTOTV) ,UCSTAR(MTOTV)  ,UPSTAR(MTOTV)  ,UMSTAR(MTOTV) ,
C     2   DUISTAR(MTOTV),DUCSTAR(MTOTV) ,DUPSTAR(MTOTV) ,DUMSTAR(MTOTV),
C     3   DUTYLR(MTOTV) ,Y(2)
      DIMENSION 
     1   DUTYLR(MTOTV)
      DATA
     1    R0   ,R1   ,R2   /
     2    0.0D0,1.0D0,2.0D0/
C   
C **********************************************************************
C SET RVE RELATED SUBMAIN VARIABLES FROM SPLMES AND SET PRESCRIBED 
C UNIFORM STRAIN COMPONENT OF ALL NODAL DISPLACEMENTS OF THE RVE MESH  
C **********************************************************************
C
C Set prescribed uniform strain component of all nodal displacements
C ==================================================================
C
      DUTYLR=R0
      DO IPOIN=1,NPOIN
        IF(NTYPE==4)THEN
          IXDOF=3*IPOIN-2
          IYDOF=IXDOF+1
          IZDOF=IYDOF+1
          XCRD=COORD(1,IPOIN,0)
          YCRD=COORD(2,IPOIN,0)
          ZCRD=COORD(3,IPOIN,0)
          IF(NLARGE==0)THEN
C                        DRV_STRAIN(1): epsilon_xx
C                        DRV_STRAIN(2): epsilon_yy
C                        DRV_STRAIN(3): epsilon_zz
C                        DRV_STRAIN(4): epsilon_xy
C                        DRV_STRAIN(5): epsilon_yz
C                        DRV_STRAIN(6): epsilon_xz
            DUTYLR(IXDOF)=DRV_STRAIN(1)*XCRD
     1                   +DRV_STRAIN(4)*YCRD
     2                   +DRV_STRAIN(6)*ZCRD
            DUTYLR(IYDOF)=DRV_STRAIN(4)*XCRD
     1                   +DRV_STRAIN(2)*YCRD
     2                   +DRV_STRAIN(5)*ZCRD
            DUTYLR(IZDOF)=DRV_STRAIN(6)*XCRD
     1                   +DRV_STRAIN(5)*YCRD
     2                   +DRV_STRAIN(3)*ZCRD
          ELSE
C                        DRV_STRAIN(1): Fxx
C                        DRV_STRAIN(2): Fxy
C                        DRV_STRAIN(3): Fxz
C                        DRV_STRAIN(4): Fyx
C                        DRV_STRAIN(5): Fyy
C                        DRV_STRAIN(6): Fyz
C                        DRV_STRAIN(7): Fzx
C                        DRV_STRAIN(8): Fzy
C                        DRV_STRAIN(9): Fzz
            DUTYLR(IXDOF)=(DRV_STRAIN(1)-R1)*XCRD
     1                   +DRV_STRAIN(2)*YCRD
     2                   +DRV_STRAIN(3)*ZCRD
            DUTYLR(IYDOF)=DRV_STRAIN(4)*XCRD
     1                   +(DRV_STRAIN(5)-R1)*YCRD
     2                   +DRV_STRAIN(6)*ZCRD
            DUTYLR(IZDOF)=DRV_STRAIN(7)*XCRD
     1                   +DRV_STRAIN(8)*YCRD
     2                   +(DRV_STRAIN(9)-R1)*ZCRD
          ENDIF
        ELSE
          IXDOF=2*IPOIN-1
          IYDOF=IXDOF+1
          XCRD=COORD(1,IPOIN,0)
          YCRD=COORD(2,IPOIN,0)
          IF(NLARGE==0)THEN
C                        DRV_STRAIN(1): epsilon_xx
C                        DRV_STRAIN(2): epsilon_yy
C                        DRV_STRAIN(3): epsilon_xy
            DUTYLR(IXDOF)=DRV_STRAIN(1)*XCRD+DRV_STRAIN(3)*YCRD
            DUTYLR(IYDOF)=DRV_STRAIN(3)*XCRD+DRV_STRAIN(2)*YCRD
          ELSE
C                        DRV_STRAIN(1): Fxx
C                        DRV_STRAIN(2): Fxy
C                        DRV_STRAIN(3): Fyx
C                        DRV_STRAIN(4): Fyy
            DUTYLR(IXDOF)=(DRV_STRAIN(1)-R1)*XCRD
     1                   +DRV_STRAIN(2)*YCRD
            DUTYLR(IYDOF)=DRV_STRAIN(3)*XCRD
     1                   +(DRV_STRAIN(4)-R1)*YCRD
          ENDIF
        ENDIF
      ENDDO
      RETURN
      END
CDOC END_SUBROUTINE SPLDSP
