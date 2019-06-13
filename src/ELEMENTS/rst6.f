CDOC BEGIN_SUBROUTINE RST6
CDOC Set properties arrays for element type TRI6
CDOC
CDOC This routine sets the element properties arrays for elements type
CDOC TRI6: Standard isoparametric 6-noded quadratic triangle for 
CDOC plane strain,plane stress and axisymmetric analysis. It also echoes
CDOC the properties to the results file and sets the unsymmetric tangent
CDOC stiffness flag.
CDOC
CDOC BEGIN_PARAMETERS
CDOC INTEGER          IELPRP <  Array of integer element properties.
CDOC INTEGER          NDATF  >  Data file unit identifier.
CDOC INTEGER          NRESF  >  Results file unit identifier.
CDOC DOUBLE_PRECISION RELPRP <  Array of real element properties.
CDOC LOGICAL          UNSYM  <  Unsymmetric tangent stiffness flag.
CDOC END_PARAMETERS
CHST
CHST 
CHST
      SUBROUTINE RST6
     1(   IELPRP     ,NDATF      ,NRESF      ,RELPRP     ,UNSYM      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER
     1(   MGAUSP=7   ,MNODEG=3   ,NDIME=2    ,NDOFEL=12  ,NEDGEL=3   ,
     2    NGAUSB=2   ,NNODE=6    )
      LOGICAL UNSYM
      DIMENSION
     1    IELPRP(*)          ,RELPRP(*)
      DIMENSION
     1    NORDEB(NNODE,NEDGEL),POSGP(2,MGAUSP)   ,POSGPB(NGAUSB)     ,
     2    WEIGP(MGAUSP)      ,WEIGPB(NGAUSB)
C***********************************************************************
C READ INPUT DATA AND SET PROPERTIES FOR ELEMENT TYPE 'TRI_6'
C (STANDARD ISOPARAMETRIC 6-NODED QUADRATIC TRIANGLE)
C***********************************************************************
 1000 FORMAT(' TRI_6 (standard 6-noded quadratic triangle)'/
     1       ' Integration rule: ',I2,' gauss points')
C
C Read number of gauss points for domain integration
C --------------------------------------------------
      READ(NDATF,*,ERR=100,END=100)NGAUSP
      WRITE(NRESF,1000)NGAUSP
      IF(NGAUSP.NE.1.AND.NGAUSP.NE.3.AND.NGAUSP.NE.6.AND.NGAUSP.NE.7)
     1  CALL ERRPRT('ED0240')
C Set element integer properties (stored in vector IELPRP)
C --------------------------------------------------------
C total number of nodes and gauss points for domain integration
      IELPRP(3)=NNODE
      IELPRP(4)=NGAUSP
C number of degrees of freedom of the element
      IELPRP(5)=NDOFEL
C number of edges of the element
      IELPRP(6)=NEDGEL
C maximum number of nodes per edge
      IELPRP(7)=MNODEG
C number of gauss points for boundary integration
      IELPRP(8)=NGAUSB
C node numbering order on boundaries (set correspondence between local
C element node numbers and "edge" node numbering for boundary
C integration)
      NORDEB(1,1)=1
      NORDEB(2,1)=2
      NORDEB(3,1)=3
      NORDEB(4,1)=0
      NORDEB(5,1)=0
      NORDEB(6,1)=0
C
      NORDEB(1,2)=0
      NORDEB(2,2)=0
      NORDEB(3,2)=1
      NORDEB(4,2)=2
      NORDEB(5,2)=3
      NORDEB(6,2)=0
C
      NORDEB(1,3)=3
      NORDEB(2,3)=0
      NORDEB(3,3)=0
      NORDEB(4,3)=0
      NORDEB(5,3)=1
      NORDEB(6,3)=2
      IPOS=9
      DO 20 IEDGEL=1,NEDGEL
        DO 10 INODE=1,NNODE
          IELPRP(IPOS)=NORDEB(INODE,IEDGEL)
          IPOS=IPOS+1
   10   CONTINUE
   20 CONTINUE
C Set element real properties (stored in vector RELPRP)
C -----------------------------------------------------
C gaussian constants for domain integration
      CALL GAUS2D
     1(  'TRI'       ,NGAUSP     ,POSGP      ,WEIGP      )
      IPOS=1
      DO 30 IGAUSP=1,NGAUSP
        RELPRP(IPOS)=POSGP(1,IGAUSP)
        RELPRP(IPOS+1)=POSGP(2,IGAUSP)
        IPOS=IPOS+NDIME
   30 CONTINUE
      IPOS=NGAUSP*NDIME+1
      DO 40 IGAUSP=1,NGAUSP
        RELPRP(IPOS)=WEIGP(IGAUSP)
        IPOS=IPOS+1
   40 CONTINUE
C set matrix of coefficients for extrapolation from gauss points to
C nodes
      IPOS=NGAUSP*NDIME+NGAUSP+1
      CALL EXT6
     1(   NGAUSP     ,RELPRP(IPOS)   )
C gaussian constants for boundary integration (integration over edges)
      CALL GAUS1D
     1(   NGAUSB     ,POSGPB     ,WEIGPB     )
      IPOS=NGAUSP*NDIME+NGAUSP+NGAUSP*NNODE+1
      DO 50 IGAUSB=1,NGAUSB
        RELPRP(IPOS)=POSGPB(IGAUSB)
        IPOS=IPOS+1
   50 CONTINUE
      IPOS=NGAUSP*NDIME+NGAUSP+NGAUSP*NNODE+NGAUSB+1
      DO 60 IGAUSB=1,NGAUSB
        RELPRP(IPOS)=WEIGPB(IGAUSB)
        IPOS=IPOS+1
   60 CONTINUE
C Set unsymmetric solver flag
C ---------------------------
      UNSYM=.FALSE.
C
      GOTO 200
C Issue error message and abort program execution in case of I/O error
  100 CALL ERRPRT('ED0241')
C
  200 CONTINUE
      RETURN
      END
CDOC END_SUBROUTINE RST6
