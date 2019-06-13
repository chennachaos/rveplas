CDOC BEGIN_SUBROUTINE SPLMES
CDOC Split RVE mesh nodes/dofs for multi-scale analysis
CDOC
CDOC Splits the nodes and degrees of freedom of a given micro-cell mesh
CDOC into different regions: interior, boundary and different parts of
CDOC boundary, according to the chosen kinematic assumption. Also 
CDOC creates the global dependency matrix required by the uniform 
CDOC boundary traction assumption.
CDOC WARNING:
CDOC   For the periodic boundary displacement assumption, this routine
CDOC   only works for a polygonal RVE with one-to-one correspondence 
CDOC   between nodes on opposing sides of the RVE (i.e. the "plus" and 
CDOC   "minus" sides)
CDOC
CDOC BEGIN_PARAMETERS
CDOC
CDOC END_PARAMETERS
CHST
CHST D. de Bortoli, April 2015: 3-D version coded, only for hexahedral
CHST                            RVEs with faces parallel to the
CHST                            coordinate planes (X-Y, X-Z and Y-Z).
CHST                            Also works for 2-D rectangular RVEs
CHST                            under analogous conditions.
CHST
      SUBROUTINE SPLMES
C     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C Hyplas global database
      INCLUDE '../MAXDIM.INC'
      INCLUDE '../MATERIAL.INC'
      INCLUDE '../ELEMENTS.INC'
      INCLUDE '../GLBDBASE.INC'
C Hyplas multi-scale database
      INCLUDE '../RVE.INC'
C Constants:
      PARAMETER
     1(   R0=0.0D0   ,  R1=1.0D0, TWOPI=8.0D0*ATAN(1.0D0))
C
C Auxiliary variables (described in relevant code sections)
      INTEGER 
     1        EGNODS(MNODE)   ,
     5        NODAUX(MPOIN)   , EXTNDS(MPOIN)   ,
     6        NODCHK(MNODE)
C
      DOUBLE PRECISION 
     3        ELCOD(MDIME,MNODE)        , THKN(MNODE)                 ,
     5        SHAPE(MNODE)    , DERIV(MDIME,MNODE) , DGASH(MDOFN)     ,
     6        EGASH(MDOFN)    , PGASH(MDOFN)       , EXTNOR(MDOFN)    ,
     7        EISCRD(MDIME)
C
      LOGICAL FOUND, SINGUL
      
      
      
      DOUBLE PRECISION
     1   DMIN(NDIME), DMAX(NDIME), DLENGT(NDIME), TOLER(NDIME)
      LOGICAL, DIMENSION(NPOIN) ::
     1   ISLEFT, ! Node belongs to...
     2   ISRIGT,
     3   ISBACK, 
     4   ISFRON,
     5   ISBOTM, 
     6   ISTOP , 
     7   ISCOR ,
     8   ISBOUN, 
     9   ISINT ,
     9   ISPLU ,
     9   ISMIN
      INTEGER
     1   NODES(NPOIN)
      
      INTEGER, ALLOCATABLE, DIMENSION(:)  :: BDNODS, INNODS
      INTEGER, ALLOCATABLE, DIMENSION(:)  :: PLNODS, MINODS, CRNODS
      INTEGER, ALLOCATABLE, DIMENSION(:)  :: CRDOFS, DEPDFS, PREDFS
      LOGICAL, ALLOCATABLE, DIMENSION(:)  :: ISDEP, ISPRES
      
      LOGICAL, ALLOCATABLE, DIMENSION(:)  :: ISMATC
      
C Logical and other arrays used to build IRV_IFILT in periodic 3-D case          
      LOGICAL, DIMENSION(12,NPOIN) :: ISEDGE
      INTEGER, DIMENSION(12) :: NNDEDG
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: NODEDG
C
      LOGICAL, DIMENSION(NTOTV) ::
     1   ISEDDF
      
      
      
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) ::
     1      CGLOB, CDEP, CDEPIN, CFREE
      
      DOUBLE PRECISION, DIMENSION(NDIME) :: CRDMIN, CRDJND
      
      DOUBLE PRECISION, PARAMETER :: GEOTOL=1.0D-8
      
      LOGICAL IS3D
C
C***********************************************************************
C SPLITTING OF NODES AND DEGREES OF FREEDOM INTO SEPARATE ARRAYS 
C (BOUNDARY, CORNERS, "PLUS"/"MINUS" SIDES, "FREE" AND DEPENDENT DOFS)
C
C Simplified algorithm:
C ----------------------------------------------------------------------
C      i. Find boundary edges (edges that appear in only one element - 
C         both internal and external boundaries are identified)
C     ii. Find an external boundary node
C    iii. Find list of external boundary edges and associated data:
C         starting at external boundary node, walk along boundary 
C         edges until a closed loop is formed, producing list of all 
C         external edges ordered in counter-clockwise sense
C   iv.a. For the periodic assumption (RVE must be a polygon):
C        1. Neighbouring edges that have same external unit normal 
C           vector are grouped together (they belong to the same 
C           geometry side of the RVE's external boundary)
C        2. For every geometry side, find the opposing side (i.e., with 
C           same length and opposing normal vector) and create a "plus" 
C           and "minus" pairing between them. Nodes associated to the
C           element edges on those sides are associated to the
C           corresponding part of the node split ("plus" or "minus")
C        3. Corner nodes are defined at intersections of adjacent 
C           geometry sides
C   iv.b. For the uniform traction assumption:
C        1. For every edge in the external boundary loop, compute the 
C           integrals necessary to build the dependency matrix
C        2. Split boundary degrees of freedom in terms of dependent and
C           free and compute final form of the dependency matrix
C
C   Assumptions:
C ----------------------------------------------------------------------
C     - Two-dimensional RVE
C     - Periodic assumption: RVE must be a polygon
C     - Periodic assumption: one-to-one correspondence between nodes at 
C       opposing sides of the RVE boundary
C     - RVE boundary is one continuous closed curve
C     - Elements' connectivities defined in counter-clockwise order
C     - Nodes numbered consecutively starting at 1:
C             1,2,3,...,NPOIN-1,NPOIN
C   
C***********************************************************************
C
C
C Check if kinematic assumption chosen is supported
      IF((IRV_RVEOPT.NE.1).AND.(IRV_RVEOPT.NE.2)
     1                    .AND.(IRV_RVEOPT.NE.3))THEN
        CALL ERRPRT('ED0402')
      ENDIF
C Check dimension of problem
      IF((NTYPE==1).OR.(NTYPE==2).OR.(NTYPE==3))THEN
        IS3D=.FALSE.
      ELSEIF(NTYPE==4)THEN
        IS3D=.TRUE.
      ELSE
        STOP 'Invalid NTYPE in SPLMES.f'
      ENDIF
C
C Initialise RVE COMMON block variables
      IRV_NDSP=0
      IRV_NDSPPT=0
      IRV_DFSP=0
      IRV_DFSPPT=0
      DRV_BNDR=R0
      IRV_IFILT=0
C
C Find mesh maximum and minimum coordinates (in X, Y and Z directions)
      DO IDIME=1,NDIME
        DMIN(IDIME)=MINVAL(COORD(IDIME,1:NPOIN,0))
        DMAX(IDIME)=MAXVAL(COORD(IDIME,1:NPOIN,0))
      ENDDO
C
C Find RVE characteristic dimension in all directions and RVE volume
      DLENGT=ABS(DMAX-DMIN)
      DRV_CELVOL=PRODUCT(DLENGT)
C 
      TOLER=GEOTOL*DLENGT     
C     
C Logical array pointing nodes that belong to each of the RVE faces
      ISLEFT=ABS(COORD(1,1:NPOIN,0)-DMIN(1)) < TOLER(1)
      ISRIGT=ABS(COORD(1,1:NPOIN,0)-DMAX(1)) < TOLER(1)
      ISBACK=ABS(COORD(2,1:NPOIN,0)-DMIN(2)) < TOLER(2)
      ISFRON=ABS(COORD(2,1:NPOIN,0)-DMAX(2)) < TOLER(2)
      IF(IS3D)THEN
        ISBOTM=ABS(COORD(3,1:NPOIN,0)-DMIN(3)) < TOLER(3)
        ISTOP =ABS(COORD(3,1:NPOIN,0)-DMAX(3)) < TOLER(3)
      ENDIF
C
C Logical arrays pointing nodes that belong to RVE boundary and RVE 
C interior
      IF(IS3D)THEN
        ISBOUN=ISLEFT.OR.ISRIGT.OR.ISBACK.OR.ISFRON.OR.ISBOTM.OR.ISTOP
      ELSE
        ISBOUN=ISLEFT.OR.ISRIGT.OR.ISBACK.OR.ISFRON
      ENDIF
      ISINT =.NOT.ISBOUN
C Corner nodes
      IF(IS3D)THEN
        ISCOR=(ISLEFT.AND.ISBACK.AND.ISBOTM)
     1    .OR.(ISLEFT.AND.ISBACK.AND.ISTOP )
     2    .OR.(ISLEFT.AND.ISFRON.AND.ISBOTM)
     3    .OR.(ISLEFT.AND.ISFRON.AND.ISTOP )
     4    .OR.(ISRIGT.AND.ISBACK.AND.ISBOTM)
     5    .OR.(ISRIGT.AND.ISBACK.AND.ISTOP )
     6    .OR.(ISRIGT.AND.ISFRON.AND.ISBOTM)
     7    .OR.(ISRIGT.AND.ISFRON.AND.ISTOP )
      ELSE
        ISCOR=(ISLEFT.AND.ISBACK)
     1    .OR.(ISLEFT.AND.ISFRON)
     2    .OR.(ISRIGT.AND.ISBACK)
     3    .OR.(ISRIGT.AND.ISFRON)
      ENDIF
C Remove corner nodes from boundary nodes
      ISBOUN=ISBOUN.AND.(.NOT.ISCOR)
C
      NCORND=COUNT(ISCOR)
      IF(IS3D.AND.(NCORND /= 8))THEN
        STOP 'Could not find 8 corner nodes'
      ELSEIF((.NOT.IS3D).AND.(NCORND /= 4))THEN
        STOP 'Could not find 4 corner nodes'
      ENDIF
      ALLOCATE(CRNODS(NCORND))
C Make corner nodes ordered
C     CRNODS=PACK(NODES,ISCOR)
C Can use nodal coordinates for that?
C CRNODS=PACK(NODES,ISCOR)
C N1= % some function of MINVAL, MAXVAL
C CRNODS=[N1,N2,N3...,N8]
C
C Node numbers need to start at 1 and there can be no gaps...!
      NODES=[(I,I=1,NPOIN)]
      IF(IS3D)THEN
        
        IF(COUNT(ISLEFT.AND.ISBACK.AND.ISBOTM)/=1)THEN
          STOP 'E1'
        ENDIF
        IF(COUNT(ISRIGT.AND.ISBACK.AND.ISBOTM)/=1)THEN
          STOP 'E2'
        ENDIF
        IF(COUNT(ISRIGT.AND.ISFRON.AND.ISBOTM)/=1)THEN
          STOP 'E3'
        ENDIF
        IF(COUNT(ISLEFT.AND.ISBACK.AND.ISBOTM)/=1)THEN
          STOP 'E4'
        ENDIF
        IF(COUNT(ISLEFT.AND.ISBACK.AND.ISTOP)/=1)THEN
          STOP 'E5'
        ENDIF
        IF(COUNT(ISRIGT.AND.ISBACK.AND.ISTOP)/=1)THEN
          STOP 'E6'
        ENDIF
        IF(COUNT(ISRIGT.AND.ISFRON.AND.ISTOP)/=1)THEN
          STOP 'E7'
        ENDIF
        IF(COUNT(ISLEFT.AND.ISBACK.AND.ISTOP)/=1)THEN
          STOP 'E8'
        ENDIF
        CRNODS(1:1)=PACK(NODES,ISLEFT.AND.ISBACK.AND.ISBOTM) ! Bottom left, bottom face
        CRNODS(2:2)=PACK(NODES,ISRIGT.AND.ISBACK.AND.ISBOTM) ! Bottom right, bottom face
        CRNODS(3:3)=PACK(NODES,ISRIGT.AND.ISFRON.AND.ISBOTM) ! Top right, bottom face
        CRNODS(4:4)=PACK(NODES,ISLEFT.AND.ISFRON.AND.ISBOTM) ! Top left, bottom face
        CRNODS(5:5)=PACK(NODES,ISLEFT.AND.ISBACK.AND.ISTOP)  ! Bottom left, top face
        CRNODS(6:6)=PACK(NODES,ISRIGT.AND.ISBACK.AND.ISTOP)  ! Bottom right, top face
        CRNODS(7:7)=PACK(NODES,ISRIGT.AND.ISFRON.AND.ISTOP)  ! Top right, top face
        CRNODS(8:8)=PACK(NODES,ISLEFT.AND.ISFRON.AND.ISTOP)  ! Top left, top face
        
      ELSE
        IF(COUNT(ISLEFT.AND.ISBACK)/=1)THEN
          STOP 'E1'
        ENDIF
        IF(COUNT(ISRIGT.AND.ISBACK)/=1)THEN
          STOP 'E2'
        ENDIF
        IF(COUNT(ISRIGT.AND.ISFRON)/=1)THEN
          STOP 'E3'
        ENDIF
        IF(COUNT(ISLEFT.AND.ISBACK)/=1)THEN
          STOP 'E4'
        ENDIF
        CRNODS(1:1)=PACK(NODES,ISLEFT.AND.ISBACK) ! Bottom left
        CRNODS(2:2)=PACK(NODES,ISRIGT.AND.ISBACK) ! Bottom right
        CRNODS(3:3)=PACK(NODES,ISRIGT.AND.ISFRON) ! Top right
        CRNODS(4:4)=PACK(NODES,ISLEFT.AND.ISFRON) ! Top left
      
      ENDIF

        
C Boundary minus corner nodes
      NBOUND=COUNT(ISBOUN)
      ALLOCATE(BDNODS(NBOUND))
      BDNODS=PACK(NODES,ISBOUN)
C
      NINTND=COUNT(ISINT)
      ALLOCATE(INNODS(NINTND))
      INNODS=PACK(NODES,ISINT)
C 
C Linear and uniform traction cases
C =================================
      IF((IRV_RVEOPT==1).OR.(IRV_RVEOPT==3))THEN
C Plus nodes: first half of boundary nodes
C Minus nodes: second half of boundary nodes
        NPLUND=NBOUND/2
        NMINND=NBOUND-NPLUND
        IRV_NDSP(1:NPOIN)=[INNODS, BDNODS, CRNODS]
C Periodic case
C =================================
      ELSEIF(IRV_RVEOPT==2)THEN
C Plus and minus node allocation
        IF(MOD(NBOUND,2)/=0)THEN
          STOP 'Uneven number of boundary nodes'
        ENDIF
        ALLOCATE(ISMATC(NBOUND))
        ISMATC=.FALSE.
C
        NPLUND=NBOUND/2
        NMINND=NPLUND     
        ALLOCATE(PLNODS(NPLUND))
        ALLOCATE(MINODS(NMINND))
C        
        IPLUND=0
        DO IBDND=1,NBOUND
          IPOIN=BDNODS(IBDND)
C Skip already matched nodes
          IF(ISMATC(IBDND))THEN
            CYCLE
          ENDIF
C Match left and right faces
          IF(ISRIGT(IPOIN))THEN
            IDIME=1
          ELSEIF(ISFRON(IPOIN))THEN
            IDIME=2
          ELSEIF(IS3D.AND.(ISTOP(IPOIN)))THEN
            IDIME=3
          ELSE
            CYCLE
          ENDIF
          CRDMIN=COORD(1:NDIME,IPOIN,0)
C Coordinates of corresponding left node
          CRDMIN(IDIME)=CRDMIN(IDIME)-DLENGT(IDIME)
C
          DO JBDND=1,NBOUND
            IF(ISMATC(JBDND))THEN
              CYCLE
            ENDIF
            JPOIN=BDNODS(JBDND)
C
            CRDJND=COORD(1:NDIME,JPOIN,0)
            IF(ALL(ABS(CRDJND-CRDMIN) <= TOLER))THEN
              ISMATC(IBDND)=.TRUE.
              ISMATC(JBDND)=.TRUE.
              IPLUND=IPLUND+1
              PLNODS(IPLUND)=IPOIN
              MINODS(IPLUND)=JPOIN
              EXIT
            ENDIF
          ENDDO
        ENDDO

C Check if all have been assigned
        IF(ANY(.NOT.ISMATC))THEN
          WRITE(*,*)PACK(BDNODS,.NOT.ISMATC)
          STOP 'Not all boundary was able to be matched'
        ENDIF
C
        IRV_NDSP(1:NPOIN)=[INNODS, PLNODS, MINODS, CRNODS]
      ENDIF
      
C Internal nodes pointers
      IRV_NDSPPT(1)=1
      IRV_NDSPPT(2)=NINTND
C Plus nodes pointers (first half of boundary nodes)
      IRV_NDSPPT(3)=IRV_NDSPPT(2)+1
      IRV_NDSPPT(4)=IRV_NDSPPT(2)+NPLUND
C Minus nodes pointers (second half of boundary nodes)
      IRV_NDSPPT(5)=IRV_NDSPPT(4)+1
      IRV_NDSPPT(6)=IRV_NDSPPT(4)+NMINND
C Corner nodes pointers
      IRV_NDSPPT(7)=IRV_NDSPPT(6)+1
      IRV_NDSPPT(8)=IRV_NDSPPT(6)+NCORND
      
C
C Perform some checks on node lists built
C ----------------------------------------------------------------------
C
C Check if total number of external boundary nodes match with the sum
C of the number of nodes in "plus", "minus" and "corner" groups
      IF((NMINND+NPLUND) /= NBOUND)THEN
        STOP 'RVE boundary not correctly divided into "plus", "minus"
     1 and corner regions'
      ENDIF
C
C Check if total number of nodes in mesh matches sum of number of
C external nodes and number of internal nodes
      IF((NINTND+NBOUND+NCORND) /= NPOIN)THEN
        STOP 'Could not find a group for every node in RVE'
      ENDIF
C
C Check if all external boundary nodes were assigned to one and only one
C region of the boundary
C
      ISPLU=.FALSE.
      IFNDSP=IRV_NDSPPT(3)
      ILNDSP=IRV_NDSPPT(4)
      DO INODSP=IFNDSP,ILNDSP
        INOD=IRV_NDSP(INODSP)
        ISPLU(INOD)=.TRUE.
      ENDDO
C
      ISMIN=.FALSE.
      IFNDSP=IRV_NDSPPT(5)
      ILNDSP=IRV_NDSPPT(6)
      DO INODSP=IFNDSP,ILNDSP
        INOD=IRV_NDSP(INODSP)
        ISMIN(INOD)=.TRUE.
      ENDDO
      
C Check for combinations of categories
C Check for absence of all categories
C ISINT
C ISPLU
C ISMIN
C ISCOR
      IF(ANY(ISINT.AND.ISPLU))THEN
        STOP '1'
      ENDIF
      IF(ANY(ISINT.AND.ISMIN))THEN
        STOP '2'
      ENDIF
      IF(ANY(ISINT.AND.ISCOR))THEN
        STOP '3'
      ENDIF
      IF(ANY(ISPLU.AND.ISMIN))THEN
        STOP '4'
      ENDIF
      IF(ANY(ISPLU.AND.ISCOR))THEN
        STOP '5'
      ENDIF
      IF(ANY(ISMIN.AND.ISCOR))THEN
        STOP '6'
      ENDIF
      IF(ANY((.NOT.ISINT).AND.(.NOT.ISPLU).AND.(.NOT.ISMIN)
     1                   .AND.(.NOT.ISCOR)))THEN
        STOP '7'
      ENDIF      
      
C From node split vectors IRV_NDSP and IRV_NDSPPT, create degree of 
C freedom split vectors IRV_DFSP and IRV_DFSPPT
C ----------------------------------------------------------------------
C
C NDFSPL: current number of degrees of freedom in the split; is used to 
C         keep track of the current position in IRV_DFSP, the vector 
C         that contains the dof numbers of the different categories 
C         ("plus", "minus", corner, ...)
C
C Internal degrees of freedom list
      NDFSPL=0
      IFNDSP=IRV_NDSPPT(1)
      ILNDSP=IRV_NDSPPT(2)
      DO INODSP=IFNDSP,ILNDSP
        INOD=IRV_NDSP(INODSP)
        DO IDIME=1,NDIME
          NDFSPL=NDFSPL+1
          IRV_DFSP(NDFSPL)=NDIME*(INOD-1)+IDIME
        ENDDO
      ENDDO
      IRV_DFSPPT(1)=1
      IRV_DFSPPT(2)=NDFSPL
C Plus degrees of freedom list
      IFNDSP=IRV_NDSPPT(3)
      ILNDSP=IRV_NDSPPT(4)
      DO INODSP=IFNDSP,ILNDSP
        INOD=IRV_NDSP(INODSP)
        DO IDIME=1,NDIME
          NDFSPL=NDFSPL+1
          IRV_DFSP(NDFSPL)=NDIME*(INOD-1)+IDIME
        ENDDO
      ENDDO
      IRV_DFSPPT(3)=IRV_DFSPPT(2)+1
      IRV_DFSPPT(4)=NDFSPL
C Minus degrees of freedom list
      IFNDSP=IRV_NDSPPT(5)
      ILNDSP=IRV_NDSPPT(6)
      DO INODSP=IFNDSP,ILNDSP
        INOD=IRV_NDSP(INODSP)
        DO IDIME=1,NDIME
          NDFSPL=NDFSPL+1
          IRV_DFSP(NDFSPL)=NDIME*(INOD-1)+IDIME
        ENDDO
      ENDDO
      IRV_DFSPPT(5)=IRV_DFSPPT(4)+1
      IRV_DFSPPT(6)=NDFSPL
C Corner degrees of freedom list
      IFNDSP=IRV_NDSPPT(7)
      ILNDSP=IRV_NDSPPT(8)
      DO INODSP=IFNDSP,ILNDSP
        INOD=IRV_NDSP(INODSP)
        DO IDIME=1,NDIME
          NDFSPL=NDFSPL+1
          IRV_DFSP(NDFSPL)=NDIME*(INOD-1)+IDIME
        ENDDO
      ENDDO
      IRV_DFSPPT(7)=IRV_DFSPPT(6)+1
      IRV_DFSPPT(8)=NDFSPL      
C
C
C ===========================================
C Boundary integration: uniform traction only
      dep_matrix: IF(IRV_RVEOPT==3)THEN
      
      
C List corner dofs
        IFCRDF=IRV_DFSPPT(7)
        ILCRDF=IRV_DFSPPT(8)
        NCRDOF=ILCRDF-IFCRDF+1
        ALLOCATE(CRDOFS(NCRDOF))
        CRDOFS=IRV_DFSP(IFCRDF:ILCRDF)
        
        ALLOCATE(ISPRES(NCRDOF))
        ALLOCATE( ISDEP(NCRDOF))
        ISPRES=.FALSE.
        ISDEP =.FALSE.
C First do IRV_NDSP and IRV_NDSPPT, then do DOF ones
C Linear/Uniform traction
C 3-D: need 6 prescribed DOFS
C 8 corner nodes -> prescribed dofs can come from corner nodes
C that leaves 24-6=18 other dofs free, that can be used to define
C dependent dofs
C           9 dependent DOFS 
C Prescribed
C 2-D: C1, X and Y -> 1,2
C      C2, Y       -> 4
C 3-D: C1, X,Y,Z   -> 1,2,3
C      C2, Y       -> 5
C      C5, X       -> 13
C      C6, Z       -> 18
        IF(IS3D)THEN
          IF(NLARGE == 1)THEN
C 9 dependent dofs: 3X (4,7,10), 3Y (8,11,14), 3Z (6,9,12)
C [1,2,3,8,22], [13,14,15,8,22]
C [1,2,3,8], [13,14,15,8]
C ISDEP ([4,5,6,7,9,10,11,12,20]), [4,5,6,10,11,12,13,14,15]
            ISPRES([1,2,3])=.TRUE.
            ISDEP ([4,5,6,7,9,10,11,12,20])=.TRUE.
          ELSE
C 6 prescribed dofs: [1,2,3,5,13,17],[1,2,3,8,22,23], [1,2,3,8,22,23], [10,11,15,16,17,18]
C 6 dependent dofs:  [4,6,7,8,9,10], [4,6,7,9,10,11], [4,6,7,9,16,14], [3,4,6,7,8,9]
C            ISPRES([10,11,15,16,17,18])=.TRUE.
C            ISDEP ([3,4,6,7,8,9])=.TRUE.
            ISPRES([1,2,3,8,22,23])=.TRUE.
            ISDEP ([4,6,7,9,16,14])=.TRUE.
          ENDIF
        ELSE
          IF(NLARGE == 1)THEN
C 4 dependent dofs:[1,2,7,8]; [1,2,3,4]; [1,6,7,8]; [2,3,4,5]; [3,4,5,6]; [3,5,6,8]
C 2-D: 3 prescribed dofs: [1,2,4]; [1,2,7]; [3,4,5]; [3,5,6]; [6,7,8]
            ISPRES([1,2])=.TRUE.
            ISDEP ([3,4,5,6])=.TRUE.
          ELSE
C 3 dependent dofs: [1,2,3]; [1,7,8]; [3,4,5]; [3,6,8]; [6,7,8]
            ISPRES([1,2,4])=.TRUE.
            ISDEP ([3,6,8])=.TRUE.
          ENDIF
        ENDIF
        
C Check for no overlaps in dependent and prescribed dofs
        IF(ANY(ISDEP.AND.ISPRES))THEN
          STOP 'Error 123'
        ENDIF
C Number of prescribed and dependent dofs
        NDEPDF=COUNT(ISDEP)
        NPREDF=COUNT(ISPRES)      
        
        ALLOCATE(DEPDFS(NDEPDF))
        ALLOCATE(PREDFS(NPREDF))
        DEPDFS=PACK(CRDOFS,ISDEP)
        PREDFS=PACK(CRDOFS,ISPRES)
        
      
C Number of boundary dofs
        NBOUDF=NDOFN*(NBOUND+NCORND)
C Number of free dofs
        NFREDF=NBOUDF-NDEPDF-NPREDF
C
        ALLOCATE(CGLOB(NDEPDF,NTOTV))
        CGLOB=R0
C
C CDEP: part of the global dependency matrix associated to the dependent
C       degrees of freedom
C CDEPIN: inverse of CDEP
C CFREE: part of the global dependency matrix associated to the free
C        degrees of freedom
C
        ALLOCATE(CDEP(NDEPDF,NDEPDF))
        ALLOCATE(CDEPIN(NDEPDF,NDEPDF))
        ALLOCATE(CFREE(NDEPDF,NFREDF))
        CDEP=R0
        CDEPIN=R0
        CFREE=R0        
        

        
C Loop through elements, listing nodes of edges
        elem_loop: DO IELEM=1,NELEM
C Set properties of the current element
          IGRUP=IGRPID(IELEM)
          IELIDN=IELTID(IGRUP)
          NNODE=IELPRP(3,IELIDN)
          NEDGEL=IELPRP(6,IELIDN)
          MNODEG=IELPRP(7,IELIDN)
          edge_loop: DO IEDGEL=1,NEDGEL
C Create list of edge nodes by going through element connectivity
C
C NNODEG: number of nodes of current edge
C EGNODS: temporary storage for ordered list of current edge node 
C         numbers
            NNODEG=0
            EGNODS=0
            DO INODEG=1,MNODEG
C IPOS: pointer to where current edge node ordering is stored in element 
C properties
              IPOS=9+(IEDGEL-1)*NNODE
              DO INODE=1,NNODE
                NODPOS=IELPRP(IPOS,IELIDN)
                IPOS=IPOS+1
                IF(NODPOS.EQ.INODEG)THEN
                  NNODEG=NNODEG+1
                  EGNODS(NNODEG)=IABS(LNODS(IELEM,INODE))
                  EXIT
                ENDIF
              ENDDO
            ENDDO
C Integrate on element edges/faces that have all nodes on external 
C boundary
            bnd_edge: IF(ALL(ISBOUN(EGNODS(1:NNODEG))
     1                   .OR.ISCOR(EGNODS(1:NNODEG))))THEN
C Set extra properties of the current element
              IELTYP=IELPRP(1,IELIDN)
              NGAUSP=IELPRP(4,IELIDN)
              NGAUSB=IELPRP(8,IELIDN)
              IPOS=9
C Check that global node numbers supplied correspond exactly to an edge
C of the current element
              NODCHK(1:NNODE)=0
              DO INODEG=1,NNODEG
                DO INODE=1,NNODE
                  IPOIN=IABS(LNODS(IELEM,INODE))
                  IF(IPOIN.EQ.EGNODS(INODEG))NODCHK(INODE)=1
                ENDDO
              ENDDO
C
              CALL CHKNDB
     1(   FOUND    ,NNODE    ,NEDGEL   ,NODCHK   ,IELPRP(IPOS,IELIDN))
              IF(.NOT.FOUND)CALL ERRPRT('ED0012')
C
C Get the global coordinates of the nodes of the loaded edge
              DO INODEG=1,NNODEG
                IPOIN=EGNODS(INODEG)
                ELCOD(1:NDIME,INODEG)=COORD(1:NDIME,IPOIN,0)
                DO II=1,NNODEG
                  IF(IABS(LNODS(IELEM,NODCHK(II))).EQ.IPOIN)THEN
                    NODAUX(IPOIN)=II
                  ENDIF
                ENDDO
              ENDDO
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
              IF(IS3D)THEN
                NDIMEB=2
              ELSE
                NDIMEB=1
              ENDIF
C
              bnd_gauss: DO IGAUSB=1,NGAUSB
C Evaluate the shape functions at the boundary sampling points
                IPPOS=NGAUSP*NDIME+NGAUSP+
     1                NGAUSP*NNODE+(IGAUSB-1)*NDIMEB+1
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
                DO IDOFN=1,NDOFN
                  PGASH(IDOFN)=R0
                  DGASH(IDOFN)=R0
                  EGASH(IDOFN)=R0
C
C Apply negative unit normal edge load and zero tangential edge load
C (we are interested in integrating the element external unit normal
C vector on the boundary)
                  IF(IDOFN.EQ.1)THEN
                    PRESS=-R1
                  ELSE
                    PRESS=R0
                  ENDIF
                  DO INODEG=1,NNODEG
C - EGNODS(INODEG) has global node number where pressure is applied
C   (local node INODEG)
C - NODAUX(IPOIN) has local number (on element boundary) of global node 
C   IPOIN
                    II=NODAUX(EGNODS(INODEG))
C Pressure components interpolated at current Gauss point using element
C boundary shape functions
                    PGASH(IDOFN)=PGASH(IDOFN)+PRESS*SHAPE(II)
C First tangent vector: derivative of position with respect to boundary 
C XI coordinate
                    DGASH(IDOFN)=DGASH(IDOFN)+
     1                           ELCOD(IDOFN,INODEG)*DERIV(1,II)
C Second tangent vector (only defined in 3-D): derivative of position
C with respect to boundary ETA coordinate
                    IF(IS3D)THEN
                      EGASH(IDOFN)=EGASH(IDOFN)+
     1                             ELCOD(IDOFN,INODEG)*DERIV(2,II)
                    ENDIF
                  ENDDO
                ENDDO
C
C Calculate components of interpolated Gauss point pressures in global 
C X,Y (and Z) coordinate system
                IF(IS3D)THEN
C In 3-D, external normal vector (not unit!) is the cross product of the
C two tangent vectors
                  EXTNOR(1)=DGASH(2)*EGASH(3)-DGASH(3)*EGASH(2)
                  EXTNOR(2)=DGASH(3)*EGASH(1)-DGASH(1)*EGASH(3)
                  EXTNOR(3)=DGASH(1)*EGASH(2)-DGASH(2)*EGASH(1)
C Norm of tangent vectors (DGASH and EGASH are not unit vectors)
                  RDGASH=NORM2(DGASH)
                  REGASH=NORM2(EGASH)
C
                  PXCOMP=-PGASH(1)*EXTNOR(1)+PGASH(2)*DGASH(1)*REGASH+
     1                    PGASH(3)*EGASH(1)*RDGASH
                  PYCOMP=-PGASH(1)*EXTNOR(2)+PGASH(2)*DGASH(2)*REGASH+
     1                    PGASH(3)*EGASH(2)*RDGASH
                  PZCOMP=-PGASH(1)*EXTNOR(3)+PGASH(2)*DGASH(3)*REGASH+
     1                    PGASH(3)*EGASH(3)*RDGASH
                ELSE
                  PXCOMP=DGASH(1)*PGASH(2)-DGASH(2)*PGASH(1)
                  PYCOMP=DGASH(1)*PGASH(1)+DGASH(2)*PGASH(2)
                ENDIF
C
                DVOLU=WEIGPB
                IF(NTYPE.EQ.1)THEN
C interpolate to find thickness at boundary gauss point (plane stress)
                  THICK=R0
                  DO INODEG=1,NNODEG
                    II=NODAUX(EGNODS(INODEG))
                    INODE=NODCHK(II)
                    THICK=THICK+THKN(INODE)*SHAPE(II)
                  ENDDO
                  DVOLU=DVOLU*THICK
                ELSEIF(NTYPE.EQ.3)THEN
C interpolate to find radius at boundary gauss point (axisymmetric case)
                  RADUS=R0
                  DO INODEG=1,NNODEG
                    II=NODAUX(EGNODS(INODEG))
                    RADUS=RADUS+SHAPE(II)*ELCOD(NAXIS,INODEG)
                  ENDDO
                  DVOLU=DVOLU*TWOPI*RADUS
                ENDIF
C
C Add the equivalent nodal loads (consistent loads) to the element force
C vector:
C contribution of current boundary Gauss point to integral of:
C (boundary shape functions) * (nodal prescribed pressures)
                bnd_node: DO INODEG=1,NNODEG
                  INODE=EGNODS(INODEG)
                  NGASH=(INODE-1)*NDOFN+1
                  MGASH=(INODE-1)*NDOFN+2
C Store integration results in global dependency matrix
                  DVOLSP=SHAPE(INODEG)*DVOLU
                  ATERM=DVOLSP*PXCOMP
                  BTERM=DVOLSP*PYCOMP
                  IF(.NOT.IS3D)THEN
C Small strains, 2-D
                    IF(NLARGE==0)THEN
                      CGLOB(1,NGASH)=CGLOB(1,NGASH)+ATERM
                      CGLOB(2,MGASH)=CGLOB(2,MGASH)+BTERM
                      CGLOB(3,NGASH)=CGLOB(3,NGASH)+BTERM
                      CGLOB(3,MGASH)=CGLOB(3,MGASH)+ATERM
C Large strains, 2-D
                    ELSE
                      CGLOB(1,NGASH)=CGLOB(1,NGASH)+ATERM
                      CGLOB(2,MGASH)=CGLOB(2,MGASH)+BTERM
                      CGLOB(3,NGASH)=CGLOB(3,NGASH)+BTERM
                      CGLOB(4,MGASH)=CGLOB(4,MGASH)+ATERM  
                    ENDIF
                  ELSE
C Z component needed only in 3-D
                    LGASH=(INODE-1)*NDOFN+3
                    CTERM=DVOLSP*PZCOMP
C Small strains, 3-D
                    IF(NLARGE==0)THEN
                      CGLOB(1,NGASH)=CGLOB(1,NGASH)+ATERM
                      CGLOB(2,MGASH)=CGLOB(2,MGASH)+BTERM
                      CGLOB(3,LGASH)=CGLOB(3,LGASH)+CTERM
                      CGLOB(4,NGASH)=CGLOB(4,NGASH)+BTERM
                      CGLOB(4,MGASH)=CGLOB(4,MGASH)+ATERM
                      CGLOB(5,MGASH)=CGLOB(5,MGASH)+CTERM
                      CGLOB(5,LGASH)=CGLOB(5,LGASH)+BTERM
                      CGLOB(6,NGASH)=CGLOB(6,NGASH)+CTERM
                      CGLOB(6,LGASH)=CGLOB(6,LGASH)+ATERM
C Large strains, 3-D
                    ELSE
                      CGLOB(1,NGASH)=CGLOB(1,NGASH)+ATERM
                      CGLOB(2,MGASH)=CGLOB(2,MGASH)+ATERM
                      CGLOB(3,LGASH)=CGLOB(3,LGASH)+ATERM
                      CGLOB(4,NGASH)=CGLOB(4,NGASH)+BTERM
                      CGLOB(5,MGASH)=CGLOB(5,MGASH)+BTERM
                      CGLOB(6,LGASH)=CGLOB(6,LGASH)+BTERM
                      CGLOB(7,NGASH)=CGLOB(7,NGASH)+CTERM
                      CGLOB(8,MGASH)=CGLOB(8,MGASH)+CTERM
                      CGLOB(9,LGASH)=CGLOB(9,LGASH)+CTERM
                    ENDIF
                  ENDIF
                ENDDO bnd_node
              ENDDO bnd_gauss
            ENDIF bnd_edge
          ENDDO edge_loop
        ENDDO elem_loop

        CDEP(:,1:NDEPDF)=CGLOB(:,[DEPDFS])
        
        SINGUL=.FALSE.
        
        CALL GAUSIN( CDEP, CDEPIN, NDEPDF, NDEPDF, SINGUL )
        
        IF(SINGUL)THEN
          STOP 'Set of dependent dofs is not suitable'
        ENDIF       
C Add all boundary dofs (not corner) to free dofs list in IRV_DFSP
        IRV_DFSPPT(9)=IRV_DFSPPT(3)
        IDFSPL=IRV_DFSPPT(9)-1
        DO IBOUND=1,NBOUND
          INODE=BDNODS(IBOUND)
          DO IDIME=1,NDIME
            IDFSPL=IDFSPL+1
            IRV_DFSP(IDFSPL)=NDIME*(INODE-1)+IDIME
          ENDDO
        ENDDO
        
C Add corner dofs which are not dependent nor prescribed        
        DO I=1,NCRDOF
          IF((.NOT.ISDEP(I)).AND.(.NOT.ISPRES(I)))THEN
            IDFSPL=IDFSPL+1
            IRV_DFSP(IDFSPL)=CRDOFS(I)
          ENDIF
        ENDDO
        IRV_DFSPPT(10)=IDFSPL
        
C Put all dependent dofs in IRV_DFSP
        IRV_DFSPPT(11)=IRV_DFSPPT(10)+1
        IRV_DFSPPT(12)=IRV_DFSPPT(10)+NDEPDF
        IRV_DFSP(IRV_DFSPPT(11):IRV_DFSPPT(12))=DEPDFS
        
C Add prescribed dofs to end of IRV_DFSP
C ADD CHECKING TO SEE IF NUMBERS CORRESPOND!
        IRV_DFSPPT(13)=IRV_DFSPPT(12)+1
        IRV_DFSPPT(14)=IRV_DFSPPT(12)+NPREDF
        IRV_DFSP(IRV_DFSPPT(13):IRV_DFSPPT(14))=PREDFS
               
      
C Final form of the dependency matrix
C ----------------------------------------------------------------------
C
C Free part of dependency matrix
        IFREMT=0
        IFREDS=IRV_DFSPPT(9)
        IFREDE=IRV_DFSPPT(10)
        DO IFREDF=IFREDS,IFREDE
          IDOF=IRV_DFSP(IFREDF)
          IFREMT=IFREMT+1
          CFREE(:,IFREMT)=CGLOB(:,IDOF)
        ENDDO
C     
C Multiply the inverse of the dependent part of the matrix with its free
C part to obtain the final dependency matrix
        DRV_BNDR(1:NDEPDF,1:NFREDF)=-MATMUL(CDEPIN,CFREE)
      

        
      ENDIF dep_matrix
      
      
      
C
C Define IRV_IFILT, the array that relates a degree of freedom's number
C in the global stiffness matrix to its number in the reduced stiffness
C matrix (see RVE.INC for more details on the definition of IRV_IFILT)
C ======================================================================
C
C 1) Linear kinematical assumption
C =================================
      IF(IRV_RVEOPT==1)THEN
C
C IRV_IFILT includes only interior dofs (boundary dofs have IRV_IFILT=0)
C
C
C Interior dofs have IRV_IFILT values starting at 1, 2, ... NDOFI
C (number of interior dofs)
        NDOFRD=0
        IDOFIS=IRV_DFSPPT(1)
        IDOFIE=IRV_DFSPPT(2)
        DO IDOFIN=IDOFIS,IDOFIE
          IDOF=IRV_DFSP(IDOFIN)
          NDOFRD=NDOFRD+1
          IRV_IFILT(IDOF)=NDOFRD
        ENDDO
C
C 2) Periodic kinematical assumption
C =================================
      ELSEIF(IRV_RVEOPT==2)THEN
C
C IRV_IFILT includes interior dofs and plus/minus dofs accordingly
C matched. Corresponding plus/minus dofs have the same value of 
C IRV_IFILT. Thus, the pairs of plus/minus dofs will have the same row/
C column numbers in the reduced stiffness matrix assembly. This 
C corresponds to summing the respective parts of the stiffness matrix 
C (Kpp+Kmm, Kpi+Kmi, etc), since the sparse solver will sum any two 
C entries of the matrix that have the same row and column numbers.
C Corner dofs have IRV_IFILT=0.
C
C
C Interior dofs have IRV_IFILT values starting at 1, 2, ... NDOFI
C (NDOFI: number of interior dofs)
        NDOFRD=0
        IDOFIS=IRV_DFSPPT(1)
        IDOFIE=IRV_DFSPPT(2)
        DO IDOFIN=IDOFIS,IDOFIE
          IDOF=IRV_DFSP(IDOFIN)
          NDOFRD=NDOFRD+1
          IRV_IFILT(IDOF)=NDOFRD
        ENDDO
C
        IF(.NOT.IS3D)THEN
C In 2-D, plus/minus dofs have IRV_IFILT values starting at:
C NDOFI+1, ..., NDOFI+NDOFP
C (NDOFP: number of plus dofs - should be equal to number of minus dofs)
C
          NDOFPL=IRV_DFSPPT(4)-IRV_DFSPPT(3)+1
          DO IDOF=1,NDOFPL
C Find IDOF-th dofs on the plus list and minus list, and assign the same
C IRV_IFILT value to both
            IDOFPL=IRV_DFSP(IRV_DFSPPT(3)+IDOF-1)
            IDOFMI=IRV_DFSP(IRV_DFSPPT(5)+IDOF-1)
            NDOFRD=NDOFRD+1
            IRV_IFILT(IDOFPL)=NDOFRD
            IRV_IFILT(IDOFMI)=NDOFRD
          ENDDO
C
        ELSE
C In 3-D, the treatment of plus and minus nodes is the same, except for
C those that belong to the nodes on the external edges of the hexahedral
C RVE (that are not corner). These have three other corresponding dofs 
C each, instead of a single one as in the 2-D case, so they are treated 
C separately:
C                        ·-------·
C     z                 /|      /|     For instance, the Y-dofs of the 
C     |                / |     / |     four nodes labelled 'X' on the
C     |               X  |    X  |     RVE to the left need to have the
C     o--- y         /   ·---/---·     same value of displacement
C    /              ·-------·   /      fluctuation, so they need to be
C   /               |  /    |  /       assigned the same value of
C  x                | X     | X        IRV_IFILT.
C                   |/      |/
C                   ·-------·
C Find nodes that are on each of the four edges along the x direction
C (edges 1, 2, 3 and 4)
          ISEDGE( 1,:)=.NOT.(ISCOR).AND.(ISBACK.AND.ISBOTM)
          ISEDGE( 2,:)=.NOT.(ISCOR).AND.(ISBACK.AND.ISTOP)
          ISEDGE( 3,:)=.NOT.(ISCOR).AND.(ISFRON.AND.ISBOTM)
          ISEDGE( 4,:)=.NOT.(ISCOR).AND.(ISFRON.AND.ISTOP)
C Find nodes that are on each of the four edges along the y direction
C (edges 5, 6, 7 and 8)
          ISEDGE( 5,:)=.NOT.(ISCOR).AND.(ISLEFT.AND.ISBOTM)
          ISEDGE( 6,:)=.NOT.(ISCOR).AND.(ISLEFT.AND.ISTOP)
          ISEDGE( 7,:)=.NOT.(ISCOR).AND.(ISRIGT.AND.ISBOTM)
          ISEDGE( 8,:)=.NOT.(ISCOR).AND.(ISRIGT.AND.ISTOP)
C Find nodes that are on each of the four edges along the z direction
C (edges 9, 10, 11 and 12)
          ISEDGE( 9,:)=.NOT.(ISCOR).AND.(ISBACK.AND.ISLEFT)
          ISEDGE(10,:)=.NOT.(ISCOR).AND.(ISBACK.AND.ISRIGT)
          ISEDGE(11,:)=.NOT.(ISCOR).AND.(ISFRON.AND.ISLEFT)
          ISEDGE(12,:)=.NOT.(ISCOR).AND.(ISFRON.AND.ISRIGT)
C Find number of nodes along each edge (not including corner nodes)
          DO IEDGE=1,12
            NNDEDG(IEDGE)=COUNT(ISEDGE(IEDGE,:))
          ENDDO
C
C All edges along x direction must have same number of nodes
          IF((NNDEDG(1)/=NNDEDG(2)).OR.(NNDEDG(1)/=NNDEDG(3)).OR.
     1       (NNDEDG(1)/=NNDEDG(4)))THEN
            STOP 'Not all edges along x direction have same no. nodes'
          ENDIF
C Same for edges along y and z directions
          IF((NNDEDG(5)/=NNDEDG(6)).OR.(NNDEDG(5)/=NNDEDG(7)).OR.
     1       (NNDEDG(5)/=NNDEDG(8)))THEN
            STOP 'Not all edges along y direction have same no. nodes'
          ENDIF
          IF((NNDEDG(9)/=NNDEDG(10)).OR.(NNDEDG(9)/=NNDEDG(11)).OR.
     1       (NNDEDG(9)/=NNDEDG(12)))THEN
            STOP 'Not all edges along z direction have same no. nodes'
          ENDIF
C
C Build arrays containing node numbers for each edge
          ALLOCATE(NODEDG(12,MAXVAL(NNDEDG)))
          NODEDG=0
          DO IEDGE=1,12
            NODEDG(IEDGE,1:NNDEDG(IEDGE))=PACK(NODES,ISEDGE(IEDGE,:))
          ENDDO
C
C Now find corresponding nodes in each edge (i.e. find which nodes
C correspond at edges x1, x2, x3 and x4) and assign their dofs the same
C value of IRV_IFILT. Also flag those dofs with ISEDDF=.TRUE., so that
C plus/minus dofs that don't belong to edges can be treated later.
          ISEDDF=.FALSE.
C
          DO IDIME=1,3
C First and last edge numbers in each dimension
            IF(IDIME==1)THEN
              IEDGES=1
              IEDGEE=4
            ELSEIF(IDIME==2)THEN
              IEDGES=5
              IEDGEE=8
            ELSEIF(IDIME==3)THEN
              IEDGES=9
              IEDGEE=12
            ENDIF
C Loop through nodes of the first edge of current dimension
            DO INODE=1,NNDEDG(IEDGES)
C IGNODE: global node number of current node of first edge
              IGNODE=NODEDG(IEDGES,INODE)
C CRD: relevant coordinate for node matching (e.g. x-coordinate if 
C matching nodes from x edges)
              CRD=COORD(IDIME,IGNODE,0)
C Assign new values of IRV_IFILT for each dof of the node
              DO IDOFN=1,NDIME
                IGDOF=NDIME*(IGNODE-1)+IDOFN
                IRV_IFILT(IGDOF)=NDOFRD+IDOFN
                ISEDDF(IGDOF)=.TRUE.
              ENDDO
C Now find corresponding nodes in each of the other three edges, by
C looping through each edge and comparing relevant coordinate of the
C nodes to that of the node on the first edge
              DO JEDGE=IEDGES+1,IEDGEE
C FOUND=.TRUE. when a matching node was found in each edge
                FOUND=.FALSE.
                DO JNODE=1,NNDEDG(JEDGE)
                  JGNODE=NODEDG(JEDGE,INODE)
C If the relevant coordinate is the same (within tolerance), the matched
C node's dofs receive the same IRV_IFILT values
                  IF(ABS(COORD(IDIME,JGNODE,0)-CRD) < TOLER(IDIME))THEN
                    DO JDOFN=1,NDIME
                      JGDOF=NDIME*(JGNODE-1)+JDOFN
                      IRV_IFILT(JGDOF)=NDOFRD+JDOFN
                      ISEDDF(JGDOF)=.TRUE.
                    ENDDO
C Flag node as matched and skip rest of loop for this edge
                    FOUND=.TRUE.
                    EXIT
                  ENDIF
                ENDDO
                IF(.NOT.FOUND)THEN
                  STOP 'Could not find matching node for edge'
                ENDIF
              ENDDO
C
              NDOFRD=NDOFRD+3
            ENDDO
          ENDDO
C         
C Plus/minus dofs that don't belong to edges are treated in the same way
C as in 2-D: they receive the same value of IRV_IFILT
C
C Number of plus dofs
          NDOFPL=IRV_DFSPPT(4)-IRV_DFSPPT(3)+1
          DO IDOF=1,NDOFPL
C IDOFPL and IDOFMI are corresponding plus/minus dofs
            IDOFPL=IRV_DFSP(IRV_DFSPPT(3)+IDOF-1)
            IDOFMI=IRV_DFSP(IRV_DFSPPT(5)+IDOF-1)
C If any of them are edge dofs, they have already receive IRV_IFILT
C values and can be skipped
            IF(ISEDDF(IDOFPL))CYCLE
            IF(ISEDDF(IDOFMI))CYCLE
            NDOFRD=NDOFRD+1
            IRV_IFILT(IDOFPL)=NDOFRD
            IRV_IFILT(IDOFMI)=NDOFRD
          ENDDO
C
        ENDIF
C
C 3) Uniform traction kinematical assumption
C =================================
      ELSEIF(IRV_RVEOPT==3)THEN
C Define IRV_IFILT so that interior dofs correspond to first linear system
C variables, with free dofs following. Dependent dofs are treated
C differently in the assembly of the linear system, being labelled by
C the array ISDEP.
C
C
        NDOFRD=0
C Pointers for first and last interior dofs
        IDOFIS=IRV_DFSPPT(1)
        IDOFIE=IRV_DFSPPT(2)
C Assign first variables of reduced system to interior dofs
        DO IDOFIN=IDOFIS,IDOFIE
          IDOF=IRV_DFSP(IDOFIN)
          NDOFRD=NDOFRD+1
          IRV_IFILT(IDOF)=NDOFRD
        ENDDO
C
C Pointers for first and last free dofs
        IDOFFS=IRV_DFSPPT(9)
        IDOFFE=IRV_DFSPPT(10)
C Assign next variables of reduced system to free dofs
        DO IDOFFR=IDOFFS,IDOFFE
          IDOF=IRV_DFSP(IDOFFR)
          NDOFRD=NDOFRD+1
          IRV_IFILT(IDOF)=NDOFRD
        ENDDO
C
C Pointers for first and last dependent dofs
        IDOFDS=IRV_DFSPPT(11)
        IDOFDE=IRV_DFSPPT(12)
C Flag dependent dofs in array ISDEP
        NDOFD=0
        DO IDOFDP=IDOFDS,IDOFDE
          IDOF=IRV_DFSP(IDOFDP)
          NDOFD=NDOFD+1
          IRV_IFILT(IDOF)=NDOFD
        ENDDO
C
C
      ENDIF

         

C Write output information (standard output and result file)
C ======================================================================
C      WRITE( *,10000)
C      WRITE(16,10000)
C10000 FORMAT(//
C     1 '  Summary from mesh splitting routine: '/
C     2 '=========================================')
C
C      WRITE( *,10010)DRV_GEOTOL,DRV_CELVOL,XLENGT,YLENGT
C      WRITE(16,10010)DRV_GEOTOL,DRV_CELVOL,XLENGT,YLENGT
C10010 FORMAT('Geometric tolerance used: ', G15.6/
C     1       '   (relative to RVE size)'//
C     2       'RVE volume (including holes): ', G15.6/
C     3       '   Characteristic dimensions: '/
C     4       '                           X: ', G15.6/
C     5       '                           Y: ', G15.6/)
CC
C      WRITE( *,10020)NTSGUN,NSEGUN,NBDEDG
C      WRITE(16,10020)NTSGUN,NSEGUN,NBDEDG
C10020 FORMAT('             Segments/edges'/
C     1       '----------------------------------------'/
C     2       '        Total boundary segments:', I8/
C     3       '     External boundary segments:', I8/
C     4       'External boundary element edges:', I8/)    
CC
CC Periodic case
C      IF(IRV_RVEOPT.EQ.2)THEN
C        WRITE( *,10030)NEXTSD,NSDPAR
C        WRITE(16,10030)NEXTSD,NSDPAR
C        DO IEXTSD=1,NEXTSD
C          WRITE( *,10031)IEXTSD
C          WRITE(16,10031)IEXTSD
C        ENDDO
C      ENDIF
C10030   FORMAT('            Periodic geometry'/
C     1         '---------------------------------------------'/
C     2         'Number of external sides identified: ', I8/
C     3         '     Number of side pairings formed: ', I8/)
C10031   FORMAT('    Side number: ', I6/
C     1         '         length: ', G15.6/
C     2         '  normal vector: ', 2G15.6/)
CC
C      WRITE( *,10040)
C      WRITE(16,10040)
CC Linear or uniform traction cases
C      IF((IRV_RVEOPT.EQ.1).OR.(IRV_RVEOPT.EQ.3))THEN
C        WRITE( *,10041)NINTND,NEXTND
C        WRITE(16,10041)NINTND,NEXTND
CC Periodic case
C      ELSEIF(IRV_RVEOPT.EQ.2)THEN
C        WRITE( *,10042)NINTND,NEXTND,NPLUND,NMINND,NCORND
C        WRITE(16,10042)NINTND,NEXTND,NPLUND,NMINND,NCORND
C      ENDIF
C10040 FORMAT('            Node split'/
C     1       '---------------------------------')
C10041 FORMAT('         Internal nodes: ', I8/
C     2       'External boundary nodes: ', I8/)
C10042 FORMAT('         Internal nodes: ', I8/
C     2       'External boundary nodes: ', I8/
C     3       '           "plus" nodes: ', I8/
C     4       '          "minus" nodes: ', I8/
C     5       '         "corner" nodes: ', I8/)
CC
C      WRITE( *,10050)
C      WRITE(16,10050)
C      NDOFIN=IRV_DFSPPT( 2)-IRV_DFSPPT( 1)+1
C      NDOFPL=IRV_DFSPPT( 4)-IRV_DFSPPT( 3)+1
C      NDOFMI=IRV_DFSPPT( 6)-IRV_DFSPPT( 5)+1
C      NDOFCO=IRV_DFSPPT( 8)-IRV_DFSPPT( 7)+1
C      NDOFFR=IRV_DFSPPT(10)-IRV_DFSPPT( 9)+1
C      NDOFDE=IRV_DFSPPT(12)-IRV_DFSPPT(11)+1
CC Linear case
C      IF(IRV_RVEOPT.EQ.1)THEN
C        WRITE( *,10051)NDOFIN
C        WRITE(16,10051)NDOFIN
CC Periodic case
C      ELSEIF(IRV_RVEOPT.EQ.2)THEN
C        WRITE( *,10052)NDOFIN,NDOFPL,NDOFMI,NDOFCO
C        WRITE(16,10052)NDOFIN,NDOFPL,NDOFMI,NDOFCO
CC Uniform traction cases
C      ELSEIF(IRV_RVEOPT.EQ.3)THEN
C        WRITE( *,10053)NDOFIN,NDOFFR,NDOFDE
C        WRITE(16,10053)NDOFIN,NDOFFR,NDOFDE
C      ENDIF
C10050 FORMAT('             Dof split'/
C     1       '---------------------------------')
C10051 FORMAT('          Internal dofs: ', I8/)
C10052 FORMAT('          Internal dofs: ', I8/
C     2       '            "plus" dofs: ', I8/
C     3       '           "minus" dofs: ', I8/
C     4       '          "corner" dofs: ', I8/)
C10053 FORMAT('          Internal dofs: ', I8/
C     2       '            "free" dofs: ', I8/
C     3       '       "dependent" dofs: ', I8/)
CC      
      END SUBROUTINE SPLMES
C
CDOC END_SUBROUTINE SPLMES
