CDOC BEGIN_SUBROUTINE INITIA
CDOC Initialises some arrays and problem control variables.
CDOC
CDOC This routine initialises various arrays and problem control
CDOC variables. Gauss-point-related quantities (such as stresses and
CDOC other state and algorithmic variables) are initialised by
CDOC calling the corresponding material interface routine.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DLAMD  <  Iterative load factor.
CDOC INTEGER          IFNEG  <  Signum (+1/-1) of the stiffness matrix
CDOC C                          determinant.
CDOC INTEGER          KLUND  <  Unloading flag.
CDOC DOUBLE_PRECISION TFACT  <  Total load factor.
CDOC DOUBLE_PRECISION TTIME  <  Time.
CDOC END_PARAMETERS
CDOC
      SUBROUTINE INITIA
     1(   DLAMD      ,IFNEG      ,KUNLD      ,TFACT      ,TTIME      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C Hyplas global database
      INCLUDE '../MAXDIM.INC'
      INCLUDE '../MATERIAL.INC'
      INCLUDE '../ELEMENTS.INC'
      INCLUDE '../GLBDBASE.INC'
C Local variables and numerical constants
      LOGICAL  LDUMMY
      DATA R0   /
     1     0.0D0/
C***********************************************************************
C INITIALISES SOME ARRAYS AND VARIABLES
C
C REFERENCE: Section 5.3.4
C***********************************************************************
      KUNLD=0
      TFACT=R0
      TTIME=R0
      DLAMD=R0
      IFNEG=1
      DO 10 IELEM=1,NELEM
        IGRUP=IGRPID(IELEM)
        IELIDN=IELTID(IGRUP)
        NEVAB=IELPRP(5,IELIDN)
        CALL RVZERO(ELOAD(1,IELEM),NEVAB)
        CALL RVZERO(ELOADO(1,IELEM),NEVAB)
   10 CONTINUE
      CALL RVZERO(DTANG,NTOTV)
      CALL RVZERO(TDISP,NTOTV)
      CALL RVZERO(TDISPO,NTOTV)
      CALL RVZERO(DINCR,NTOTV)
      CALL RVZERO(DINCRO,NTOTV)
      CALL RVZERO(DITER,NTOTV)
C Arrays from common block STATE
      DO 30 IELEM=1,NELEM
        IGRUP=IGRPID(IELEM)
        IELIDN=IELTID(IGRUP)
        NGAUSP=IELPRP(4,IELIDN)
        DO 20 IGAUSP=1,NGAUSP
C Call material interface routine to initialise material-specific Gauss
C point data
          MODE=0
          CALL MATISW
     1(   MODE       ,NLARGE     ,NTYPE      ,
     2    IPROPS(1,MATTID(IGRUP)),LALGVA(1,IGAUSP,IELEM,1)           ,
     3    LDUMMY     ,RALGVA(1,IGAUSP,IELEM,1)           ,DUMMY      ,
     4    RPROPS(1,MATTID(IGRUP))            ,
     5    RSTAVA(1,IGAUSP,IELEM,1)           ,DUMMY      ,
     6    STRSG(1,IGAUSP,IELEM,1)            ,DUMMY      )
   20   CONTINUE
   30 CONTINUE
      RETURN
      END
CDOC END_SUBROUTINE INITIA
