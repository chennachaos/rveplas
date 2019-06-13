CDOC BEGIN_SUBROUTINE RSH8
CDOC Read input and set properties for element type HEXA_8
CDOC
CDOC This routine reads data from the input data file and sets the
CDOC element properties arrays for elements type HEXA_8: standard
CDOC isoparametric 8-noded hexahedron for three-dimensional analysis. It 
CDOC also echoes the properties to the results file and sets the 
CDOC unsymmetric tangent stiffness flag.
CDOC
CDOC BEGIN_PARAMETERS
CDOC INTEGER          IELPRP <  Array of integer element properties.
CDOC INTEGER          NDATF  >  Data file unit identifier.
CDOC INTEGER          NRESF  >  Results file unit identifier.
CDOC DOUBLE_PRECISION RELPRP <  Array of real element properties.
CDOC LOGICAL          UNSYM  <  Unsymmetric tangent stiffness flag.
CDOC END_PARAMETERS
CHST
CHST D. de Bortoli, March 2015: Initial coding, based on RSQ4
CHST
      SUBROUTINE RST7
     1(   IELPRP     ,NDATF      ,NRESF      ,RELPRP     ,UNSYM      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER
     1(   MGAUSP=4   ,MNODEG=3    ,NDIME=3    ,NDOFEL=30   ,NEDGEL=4   ,
     2    NGAUSB=3   ,NNODE=10    ,NDIMEB=2)
      LOGICAL UNSYM
      DIMENSION
     1    IELPRP(*)          ,RELPRP(*)
      DIMENSION
     1    NORDEB(NNODE,NEDGEL)      ,POSGP(NDIME,MGAUSP),
     2    POSGPB(NDIMEB,NGAUSB)     ,WEIGP(MGAUSP)      ,WEIGPB(NGAUSB)
C***********************************************************************
C READ INPUT DATA AND SET PROPERTIES FOR ELEMENT TYPE 'HEXA_8'
C (STANDARD ISOPARAMETRIC 8-NODED TRI-LINEAR HEXAHEDRON)
C
C
C                        4-------3
C     ETA               /|      /|     STANDARD ISOPARAMETRIC
C     |                / |     / |     TRI-LINEAR 8-NODE HEXAHEDRON
C     |               /  |    /  |
C     o--- XI        /   1---/---2 
C    /              8-------7   /
C   /               |  /    |  /
C ZETA              | /     | /
C                   |/      |/
C                   5-------6
C
C
C REFERENCE: Figure (6.14), Zienkiewicz, O.C., Taylor, R.L., 
C            Zhu, J.Z., The Finite Element Method: Its Basis and 
C            Fundamentals. 7th ed. Butterworth-Heinemann, 2013.
C***********************************************************************
 1000 FORMAT(' HEXA_8 (standard 8-noded hexahedron)'/
     1       ' Integration rule: ',I2,' gauss points')
C
C Read number of gauss points for domain integration
C --------------------------------------------------
      READ(NDATF,*,ERR=100,END=100)NGAUSP
      WRITE(NRESF,1000)NGAUSP
      IF((NGAUSP.NE.1).AND.(NGAUSP.NE.4))THEN
        CALL ERRPRT('ED0132')
      ENDIF
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
      NORDEB(2,1)=3
      NORDEB(3,1)=5
      NORDEB(4,1)=0
      NORDEB(5,1)=2
      NORDEB(6,1)=4
      NORDEB(7,1)=6
      NORDEB(8,1)=0
      NORDEB(9,1)=0
      NORDEB(10,1)=0
C
      NORDEB(1,2)=0
      NORDEB(2,2)=1
      NORDEB(3,2)=5
      NORDEB(4,2)=3
      NORDEB(5,2)=0
      NORDEB(6,2)=6
      NORDEB(7,2)=0
      NORDEB(8,2)=0
      NORDEB(9,2)=2
      NORDEB(10,2)=4
C
      NORDEB(1,3)=1
      NORDEB(2,3)=0
      NORDEB(3,3)=5
      NORDEB(4,3)=3
      NORDEB(5,3)=0
      NORDEB(6,3)=0
      NORDEB(7,3)=6
      NORDEB(8,3)=2
      NORDEB(9,3)=0
      NORDEB(10,3)=4
C
      NORDEB(1,4)=1
      NORDEB(2,4)=3
      NORDEB(3,4)=0
      NORDEB(4,4)=5
      NORDEB(5,4)=2
      NORDEB(6,4)=0
      NORDEB(7,4)=0
      NORDEB(8,4)=6
      NORDEB(9,4)=4
      NORDEB(10,4)=0
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
      CALL GAUS3D
     1(  'TET'       ,NGAUSP     ,POSGP      ,WEIGP      )
      IPOS=0
      DO IGAUSP=1,NGAUSP
        DO IDIME=1,NDIME
          IPOS=IPOS+1
          RELPRP(IPOS)=POSGP(IDIME,IGAUSP)
        ENDDO
      ENDDO
      IPOS=NGAUSP*NDIME+1
      DO 40 IGAUSP=1,NGAUSP
        RELPRP(IPOS)=WEIGP(IGAUSP)
        IPOS=IPOS+1
   40 CONTINUE
C set matrix of coefficients for extrapolation from gauss points to
C nodes
      IPOS=NGAUSP*NDIME+NGAUSP+1
      CALL EXT7
     1(   NGAUSP     ,RELPRP(IPOS))
C gaussian constants for boundary integration (integration over faces - 
C quadrilaterals)
      CALL GAUS2D
     1(   'TRI',     NGAUSB     ,POSGPB     ,WEIGPB     )
      IPOS=NGAUSP*NDIME+NGAUSP+NGAUSP*NNODE+1
      DO 50 IGAUSB=1,NGAUSB
        DO 55 IDIMEB=1,NDIMEB
          RELPRP(IPOS)=POSGPB(IDIMEB,IGAUSB)
          IPOS=IPOS+1
   55   CONTINUE
   50 CONTINUE
      IPOS=NGAUSP*NDIME+NGAUSP+NGAUSP*NNODE+NGAUSB*NDIMEB+1
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
  100 CALL ERRPRT('ED0237')
C
  200 CONTINUE
      RETURN
      END
CDOC END_SUBROUTINE RST7
