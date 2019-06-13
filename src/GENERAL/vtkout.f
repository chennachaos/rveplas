C Create VTK file (for Paraview)
      SUBROUTINE VTKOUT ( IINCS, TFACT, DATFIL, DUTYLR )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C Hyplas database
      INCLUDE '../MAXDIM.INC'
      INCLUDE '../MATERIAL.INC'
      INCLUDE '../ELEMENTS.INC'
      INCLUDE '../GLBDBASE.INC'
      INCLUDE '../RVE.INC'
C Local arrays and variables
      PARAMETER( MGDIM=9 )
C
      DOUBLE PRECISION 
     1        ELCOD(NDIME,MNODE)  ,SHAPEF(MNODE)   ,DERIV(NDIME,MNODE),
     2        EISCRD(NDIME)       ,CARTD(NDIME,MNODE)  , PSTRS(3)
      
      DOUBLE PRECISION DUTYLR(MTOTV)


      
      CHARACTER*256     DATFIL, VTKFIL, TEMP
      

      
C Internal representation of elements in VTK file format
      INTEGER, PARAMETER :: NELTYP=10 !Number of element types currently in HYPLAS
      INTEGER VTKCOD(NELTYP)
      INTEGER NODORD(MNODE,NELTYP)
C
      CHARACTER*2       STRCMP(6)
      CHARACTER*3       PSTRCM(3)
C
      DOUBLE PRECISION  STRSN(MSTRE,MNODE), RSTAN(MRSTAV,MNODE),
     1                  RALGN(MRALGV,MNODE)
C Needed for weighted variable averaging (based on element area/volume)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) ::  ELVOL
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) ::  TELVOL
C Variable averaging (not by group)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) ::  STRSA, PSTRSA
C      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) ::  RSTAA, RALGA
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) ::  PRESS, EFFSTR
C           
C Variable averaging (by group)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: 
     1      STRSAG, RSTAAG, RALGAG
      
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) ::  PSTRAG
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: PRESSG, EFFSTG
C      
      INTEGER, PARAMETER :: NMTTYP=14 !Number of material types currently in HYPLAS
      CHARACTER*20      :: RSTACM(MRSTAV,NMTTYP), RALGCM(MRALGV,NMTTYP)
      DOUBLE PRECISION  :: DVLSTA(MRSTAV,NMTTYP), DVLALG(MRALGV,NMTTYP)
      CHARACTER*20      :: MATNAM(NMTTYP)
      
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: VARCOD
      CHARACTER*20      :: ALLVAR(MRSTAV*NMTTYP)
      
      INTEGER           :: NRSTVM(NMTTYP), NRALVM(NMTTYP)
      LOGICAL           :: ISMATP(NMTTYP)
C
      LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: ISGRND
      INTEGER, ALLOCATABLE, DIMENSION(:) :: OLDNOD, NWNDGR
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: NEWNOD
      
C Needed for discontinuous variable output
      INTEGER, ALLOCATABLE, DIMENSION(:,:) ::  NLNODS
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: 
     1                           STRSD, RSTAD, RALGD, PSTRD
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: PRESSD, EFFSTD
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: MYTMPVAR

      DATA R0    ,R1    ,R2    ,R3    /
     1     0.0D0 ,1.0D0 ,2.0D0 ,3.0D0 /

C MATERIAL VARIABLES
C ------------------
      MATNAM=''
      NRSTVM=0
      RSTACM=''
      DVLSTA=R0
      NRALVM=0
      RALGCM=''
      DVLALG=R0
      ALLOCATE(MYTMPVAR(10))
C
C ELASTIC
C =======
C
      MATNAM(ELASTC)='ELASTIC'
C Real state variables (RSTAVA)
C -----------------------------
C 1-4 (2D) / 1-6 (3D): infinitesimal eng. strains (small strains)
C                      logarithmic eng. strains (large strains)
      IF((NTYPE==1).OR.(NTYPE==2).OR.(NTYPE==3))THEN
        NRSTVM(ELASTC)=4
        RSTACM(1:4,ELASTC)=['E_xx','E_yy','E_xy','E_zz']
      ELSEIF(NTYPE==4)THEN
        NRSTVM(ELASTC)=6
        RSTACM(1:6,ELASTC)=['E_xx','E_yy','E_zz','E_xy','E_yz','E_xz']
      ENDIF
C
C OGDEN
C =======
C
      MATNAM(OGDEN)='OGDEN'
C Real state variables (RSTAVA)
C -----------------------------
C 1-4: left Cauchy-Green tensor
      NRSTVM(OGDEN)=4
      RSTACM(1:4,OGDEN)=['b_xx','b_yy','b_xy','b_zz']
      DVLSTA([1,2,4],OGDEN)=R1
C
C VON MISES i changed it to drucker
C =========
C
      MATNAM(DRUPRA)='DRUCKER_PRAGER'
C Real state variables (RSTAVA)
C -----------------------------
C 1-4: infinitesimal eng. strains (small strains)
C      logarithmic eng. strains (large strains)
C   5: accumulated plastic strain
      IF(NTYPE==2)THEN
        NRSTVM(DRUPRA)=5
        RSTACM(1:4,DRUPRA)=['Ee_xx','Ee_yy','Ee_zz','Ee_xy']
        RSTACM(5,DRUPRA)='Eps'
      ELSE
      NRSTVM(DRUPRA)=7
      RSTACM(1:6,DRUPRA)=['Ee_xx','Ee_yy','Ee_zz','Ee_xy',
     1                    'Ee_yz','Ee_xz']
      RSTACM(7,DRUPRA)='Eps'
C Real algorithmic variables (RALGVA)
      ENDIF
C -----------------------------------
C   1: incremental plastic multiplier
      NRALVM(DRUPRA)=1
      RALGCM(1:1,DRUPRA)=['dgama']
C
C PLANAR DOUBLE SLIP ELASTO-PLASTIC SINGLE CRYSTAL
C ================================================
C
      MATNAM(PDSCRY)='ELASTO_PLANAR_CRYSTAL'
C Real state variables (RSTAVA)
C -----------------------------
C 1-4: elastic deformation gradient (in plane components)
      NRSTVM(PDSCRY)=5
      RSTACM(1:4,PDSCRY)=['Fe-xx','Fe-yx','Fe-xy','Fe-yy']
      DVLSTA([1,4],PDSCRY)=R1
C 5: accumulated plastic slip (hardening variable)
      RSTACM(5,PDSCRY)='hrvar'
C Real algorithmic variables (RALGVA)
C -----------------------------------
C 1-4: incremental plastic multipliers for each 4 slip systems
      NRALVM(PDSCRY)=4
      RALGCM(1:4,PDSCRY)=['dgama1','dgama2','dgama3','dgama4']
C
C VISCO-PLASTIC SINGLE CRYSTAL
C ============================
C
      MATNAM(VSCTWD)='VISCO_CRYSTAL'
C Real state variables (RSTAVA)
C -----------------------------
      NRSTVM(VSCTWD)=19
C 1-9: elastic deformation gradient
      RSTACM(1:9,VSCTWD)=['Fe-xx','Fe-yx','Fe-zx','Fe-xy','Fe-yy',
     1                    'Fe-zy','Fe-xz','Fe-yz','Fe-zz']
      DVLSTA([1,5,9],VSCTWD)=R1
C 10-18: plastic deformation gradient
      RSTACM(10:18,VSCTWD)=['Fp-xx','Fp-yx','Fp-zx','Fp-xy','Fp-yy',
     1                      'Fp-zy','Fp-xz','Fp-yz','Fp-zz']
      DVLSTA([10,14,18],VSCTWD)=R1
C 19: accumulated plastic strain (hardening variable)
      RSTACM(19,VSCTWD)='hvar'
C
C VISCO-PLASTIC MARTENSITIC TRANSFORMATION SINGLE CRYSTAL
C ======================================================
C
      MATNAM(MTEPTD)='MARTEN_VISCO'
C Real state variables (RSTAVA)
C -----------------------------
      NRSTVM(MTEPTD)=52
C 1-9: elastic deformation gradient
      RSTACM(1:9,MTEPTD)=['Fe-xx','Fe-yx','Fe-zx','Fe-xy','Fe-yy',
     1                    'Fe-zy','Fe-xz','Fe-yz','Fe-zz']
      DVLSTA([1,5,9],MTEPTD)=R1
C 10-18: transformation deformation gradient
      RSTACM(10:18,MTEPTD)=['Ftr-xx','Ftr-yx','Ftr-zx','Ftr-xy','Ftr-yy'
     1                     ,'Ftr-zy','Ftr-xz','Ftr-yz','Ftr-zz']
      DVLSTA([10,14,18],MTEPTD)=R1
C 19: transformation multiplier - gamma
      RSTACM(19,MTEPTD)='gamma'
C 20-28: work conjugate of transformation deformation gradient
      RSTACM(20:28,MTEPTD)=['Tn-xx','Tn-yx','Tn-zx','Tn-xy','Tn-yy',
     1                      'Tn-zy','Tn-xz','Tn-yz','Tn-zz']
C 29-37: total deformation gradient
      RSTACM(29:37,MTEPTD)=['Fto-xx','Fto-yx','Fto-zx','Fto-xy','Fto-yy'
     1                     ,'Fto-zy','Fto-xz','Fto-yz','Fto-zz']
      DVLSTA([29,33,37],MTEPTD)=R1
C 38: differential volume (for volume fraction calculation) - dvolu
      RSTACM(38,MTEPTD)='dvolu'
C 39: active variant of the transformation system (integer) - iactrs
      RSTACM(39,MTEPTD)='iactrs'
C 40: whether further variant selection is required
      RSTACM(40,MTEPTD)='ifvars'
C 41: hardening internal variable (accumulated plastic slip?) - hrvar
      RSTACM(41,MTEPTD)='hrvar'
C 42-50: austenite plastic deformation gradient
      RSTACM(42:50,MTEPTD)=['Fpa-xx','Fpa-yx','Fpa-zx','Fpa-xy','Fpa-yy'
     1                     ,'Fpa-zy','Fpa-xz','Fpa-yz','Fpa-zz']
      DVLSTA([42,46,50],MTEPTD)=R1
C 51: state of transformation in previous converged solution (yes/no)
      RSTACM(51,MTEPTD)='ifprtr'
C 52: state of transformation (yes/no)
      RSTACM(52,MTEPTD)='iftran'
C
C
C ======================================================================
C
C VISCO-PLASTIC MARTENSITIC TRANSFORMATION SINGLE CRYSTAL (MODIFIED)
C ======================================================
C
      MATNAM(MTCTWD)='MARTEN_VISCO2'
C Real state variables (RSTAVA)
C -----------------------------
      NRSTVM(MTCTWD)=44
C 1-9: elastic deformation gradient
      RSTACM(1:9,MTCTWD)=['Fe-xx','Fe-yx','Fe-zx','Fe-xy','Fe-yy',
     1                    'Fe-zy','Fe-xz','Fe-yz','Fe-zz']
      DVLSTA([1,5,9],MTCTWD)=R1
C 10-18: transformation deformation gradient
      RSTACM(10:18,MTCTWD)=['Fpa-xx','Fpa-yx','Fpa-zx','Fpa-xy','Fpa-yy'
     1                     ,'Fpa-zy','Fpa-xz','Fpa-yz','Fpa-zz']
      DVLSTA([10,14,18],MTCTWD)=R1
C 19: hardening internal variable (accumulated plastic slip) - hrvar
      RSTACM(19,MTCTWD)='hrvar'
C 20-43: volume fractions martensite
      RSTACM(20:28,MTCTWD)=['mvf1','mvf2','mvf3','mvf4','mvf5','mvf6',
     1  'mvf7','mvf8','mvf9']
      RSTACM(29:43,MTCTWD)=['mvf10','mvf11','mvf12','mvf13','mvf14',
     2  'mvf15','mvf16','mvf17','mvf18','mvf19','mvf20','mvf21','mvf22',
     3  'mvf23','mvf24']
C 44: total martensite volume fraction
      RSTACM(44,MTCTWD)='gamma'
C
C
C ======================================================================
C
C Element types (VTK codes) and node ordering
      NODORD=0
C Linear triangles - VTK code 5 (VTK_TRIANGLE)
      VTKCOD(TRI3)=5
      NODORD(1:3,TRI3)=[(I,I=1,3)]
C Quadratic triangles - VTK code 22 (VTK_QUADRATIC_TRIANGLE)
C TRI6 elements in VTK have a different node order than in HYPLAS
      VTKCOD(TRI6)=22
      NODORD(1:6,TRI6)=[1,3,5,2,4,6]
      VTKCOD(TETA7)=24
      NODORD(1:10, TETA7)=[1,2,3,4,5,6,7,8,9,10]
C Linear quadrilaterals - VTK code 9 (VTK_QUAD)
      VTKCOD( QUAD4)=9
      VTKCOD(QUA4FB)=9
      NODORD(1:4, QUAD4)=[(I,I=1,4)]
      NODORD(1:4,QUA4FB)=[(I,I=1,4)]
C Quadratic quadrilaterals - VTK code 23 (VTK_QUADRATIC_QUAD)
C QUAD8 elements in VTK have a different node order than in HYPLAS
      VTKCOD( QUAD8)=23 
      NODORD(1:8,QUAD8)=[1,3,5,7,2,4,6,8]
C Linear hexahedra - VTK code 12 (VTK_HEXAHEDRON)
      VTKCOD( HEXA8)=12
      VTKCOD(HEX8FB)=12
      NODORD(1:8, HEXA8)=[(I,I=1,8)]
      NODORD(1:8,HEX8FB)=[(I,I=1,8)]
      
      
      IF(NTYPE==1)THEN
        NSTRE=3
        STRCMP(1:NSTRE)=['xx','yy','xy']
      ELSEIF(NTYPE==2)THEN
        NSTRE=4
        STRCMP(1:NSTRE)=['xx','yy','xy','zz']
      ELSEIF(NTYPE==4)THEN
        NSTRE=6
        STRCMP(1:NSTRE)=['xx','yy','zz','xy','yz','xz']
      ENDIF
      IF((NTYPE==1).OR.(NTYPE==2).OR.(NTYPE==3))THEN
        PSTRCM=['max','min','ang']
      ELSEIF(NTYPE==4)THEN
        PSTRCM=['max','mid','min']
      ENDIF
      
C Loop through elements, calculating their (current) volume
      write(*,*) "NELEM=", NELEM
      ALLOCATE(ELVOL(NELEM))
      ELVOL=R0
      DO IELEM=1,NELEM
        IGRUP =IGRPID(IELEM)
        IELIDN=IELTID(IGRUP)
        IELTYP=IELPRP(1,IELIDN)
        NNODE =IELPRP(3,IELIDN)
        NGAUSP=IELPRP(4,IELIDN)
C Compute coordinates of the element nodal points
        CALL RVZERO(ELCOD,NDIME*MNODE)
        DO INODE=1,NNODE
          LNODE=IABS(LNODS(IELEM,INODE))
          DO IDIME=1,NDIME
            ELCOD(IDIME,INODE)=COORD(IDIME,LNODE,1)
          ENDDO
        ENDDO
C Loop for numerical integration over element domain
C
        IPWEI=NGAUSP*NDIME
        DO IGAUSP=1,NGAUSP
          IPPOS=NDIME*(IGAUSP-1)
          DO IDIME=1,NDIME
            EISCRD(IDIME)=RELPRP(IPPOS+IDIME,IELIDN)
          ENDDO
          WEIGP=RELPRP(IPWEI+IGAUSP,IELIDN)
C Compute the shape functions at the sampling points and elemental
C volume
          CALL SHPFUN
     1(   DERIV      ,EISCRD      ,0          ,IELTYP     ,
     2    NDIME      ,SHAPEF      )
C     
          CALL JACDET
     1(   CARTD      ,DERIV      ,DETJAC     ,ELCOD     ,IELEM      ,
     2    NDIME      ,NDIME      ,NNODE      )
      
          DVOLU=DETJAC*WEIGP
          IF(NTYPE.EQ.1)THEN
            DVOLU=DVOLU*THKGP(IGAUSP,IELEM,1)
          ENDIF
          
          ELVOL(IELEM)=ELVOL(IELEM)+DVOLU
        ENDDO
      ENDDO
C
C TELVOL(IPOIN, IGRUP): total volume of elements in group IGRUP that 
C share node IPOIN
      ALLOCATE(TELVOL(NPOIN,NGRUP))
      TELVOL=R0
      DO IELEM=1,NELEM
        IGRUP=IGRPID(IELEM)
        IELIDN=IELTID(IGRUP)
        NNODE=IELPRP(3,IELIDN)
        DO INODE=1,NNODE
          LNODE=IABS(LNODS(IELEM,INODE))
          TELVOL(LNODE,IGRUP)=TELVOL(LNODE,IGRUP)+ELVOL(IELEM)
        ENDDO
      ENDDO
C
C Generate .vtk file based on input file name with results from current
C load increment
C ----------------------------------------------------------------------
      IUNIT=19
      VTKFIL=''
      I=INDEX(DATFIL,' ')-1
      VTKFIL(1:I-3)=DATFIL(1:I-3)
      WRITE(TEMP,'(I4.4)')IINCS
      VTKFIL(I-3:I+6)='__'//TRIM(TEMP)//'.vtk'
C AVERAGING BY GROUP
C ======================================================================
C Averaging by group: 
C Find nodes that are at interfaces of groups (they might be at more than one interface)
C They need to be copied as many times as they appear in groups
C        
C NEWNOD: for node IPOIN in an element belonging to group IGRUP, stores its new node number
C Loop through elements, writing down instead of connectivity:
C OLD: LNODS(IELEM,1:NNODE))
C NEW: NEWNOD(LNODS(IELEM,1:NNODE)),IGRUP)
      OPEN(UNIT=IUNIT,FILE=VTKFIL,STATUS='UNKNOWN')
C File header
      WRITE(IUNIT,'(A)')'# vtk DataFile Version 3.1'
      WRITE(IUNIT,'(A)')DATFIL(1:I)//' - results from increment '
     1                          //TRIM(TEMP)
      WRITE(IUNIT,'(A)')'ASCII'
      WRITE(IUNIT,'(A)')'DATASET UNSTRUCTURED_GRID'
C
      ALLOCATE(ISGRND(NPOIN,NGRUP))
      ISGRND=.FALSE.
      DO IGRUP=1,NGRUP
        DO IELEM=1,NELEM
          LGRUP=IGRPID(IELEM)
          IF(LGRUP.NE.IGRUP)CYCLE
          IELIDN=IELTID(IGRUP)
          NNODE =IELPRP(3,IELIDN)
          DO INODE=1,NNODE
            IPOIN=IABS(LNODS(IELEM,INODE))
            ISGRND(IPOIN,IGRUP)=.TRUE.
          ENDDO
        ENDDO
      ENDDO
C New number of nodes
      NEWNPO=COUNT(ISGRND)
      
C Each new node will only belong to elements of a particular group
C Need to know to which group a node belongs
      ALLOCATE(NWNDGR(NEWNPO))
      ALLOCATE(OLDNOD(NEWNPO))
      ALLOCATE(NEWNOD(NPOIN,NGRUP))
      NWNDGR=0
      OLDNOD=0
      NEWNOD=0
C
      INEWND=0
      DO IPOIN=1,NPOIN
        DO IGRUP=1,NGRUP
C If node belongs to a particular group
          IF(ISGRND(IPOIN,IGRUP))THEN
            INEWND=INEWND+1
            OLDNOD(INEWND)=IPOIN
            NWNDGR(INEWND)=IGRUP
            NEWNOD(IPOIN,IGRUP)=INEWND
          ENDIF
        ENDDO
      ENDDO
C 
C Nodal coordinates
C -----------------
      WRITE(IUNIT,'(/A,I0,A)')'POINTS ',NEWNPO,' double'
      WRITE(IUNIT,'(3E24.15e3)')COORD(:,OLDNOD(1:NEWNPO),0)
      
C Element connectivities
C ----------------------
C Calculate total number of values to write (NVAL)
      NVAL=0
      DO IELEM=1,NELEM
        NNODE=IELPRP(3,IELTID(IGRPID(IELEM)))
        NVAL=NVAL+NNODE+1
      ENDDO
      WRITE(IUNIT,'(/A,I0,A,I0)')'CELLS ',NELEM,' ',NVAL
      DO IELEM=1,NELEM
        IGRUP=IGRPID(IELEM)
        IELIDN=IELTID(IGRUP)
        NNODE=IELPRP(3,IELIDN)
        IELTYP=IELPRP(1,IELIDN)
        WRITE(IUNIT,'(I0,*(1X, I0))')NNODE,
     1          NEWNOD(ABS(LNODS(IELEM,NODORD(1:NNODE,IELTYP))),IGRUP)-1
      ENDDO
C Element types
      WRITE(IUNIT,'(/A,I0)')'CELL_TYPES ',NELEM
      WRITE(IUNIT,'(25(I0,1X))')
     1               VTKCOD(IELPRP(1,IELTID(IGRPID(1:NELEM))))
      
C Point data
C ==========
      WRITE(IUNIT,'(/A,I0)')'POINT_DATA ',NEWNPO
C Displacements
C -------------
      WRITE(IUNIT,'(/A)')'VECTORS Displacement double'
      DO INEWNP=1,NEWNPO
        IPOIN=OLDNOD(INEWNP)
        ISTART=1+NDOFN*(IPOIN-1)
        IEND  =NDOFN*IPOIN
C 2-D analysis
        IF((NTYPE==1).OR.(NTYPE==2).OR.(NTYPE==3))THEN
          WRITE(IUNIT,'(3E24.15e3)')TDISP(ISTART:IEND),R0
C 3-D analysis
        ELSEIF(NTYPE==4)THEN
          WRITE(IUNIT,'(3E24.15e3)')TDISP(ISTART:IEND)
        ENDIF
      ENDDO
C Taylor displacements
C -------------------------
      IF(NMULTI==2)THEN
      WRITE(IUNIT,'(/A)')'VECTORS DisplacementTaylor double'
      DO INEWNP=1,NEWNPO
        IPOIN=OLDNOD(INEWNP)
        ISTART=1+NDOFN*(IPOIN-1)
        IEND  =NDOFN*IPOIN
C 2-D analysis
        IF((NTYPE==1).OR.(NTYPE==2).OR.(NTYPE==3))THEN
          WRITE(IUNIT,'(3E24.15e3)')TFACT*DUTYLR(ISTART:IEND),R0
C 3-D analysis
        ELSEIF(NTYPE==4)THEN
          WRITE(IUNIT,'(3E24.15e3)')TFACT*DUTYLR(ISTART:IEND)
        ENDIF
      ENDDO
C Displacement fluctuations
C -------------------------
      WRITE(IUNIT,'(/A)')'VECTORS DisplacementFluctuation double'
      DO INEWNP=1,NEWNPO
        IPOIN=OLDNOD(INEWNP)
        ISTART=1+NDOFN*(IPOIN-1)
        IEND  =NDOFN*IPOIN
C 2-D analysis
        IF((NTYPE==1).OR.(NTYPE==2).OR.(NTYPE==3))THEN
          WRITE(IUNIT,'(3E24.15e3)')TDISP(ISTART:IEND)-
     1                              TFACT*DUTYLR(ISTART:IEND),R0
C 3-D analysis
        ELSEIF(NTYPE==4)THEN
          WRITE(IUNIT,'(3E24.15e3)')TDISP(ISTART:IEND)-
     1                              TFACT*DUTYLR(ISTART:IEND)
        ENDIF
      ENDDO
      ENDIF
C
C Stress averaging (by element group)
      ALLOCATE(STRSAG(MSTRE,NEWNPO))
      STRSAG=R0
      ALLOCATE(RSTAAG(MRSTAV,NEWNPO))
      RSTAAG=R0
      ALLOCATE(RALGAG(MRALGV,NEWNPO))
      RALGAG=R0
      DO IGRUP=1,NGRUP
        DO IELEM=1,NELEM
          LGRUP=IGRPID(IELEM)
          IF(LGRUP.NE.IGRUP)CYCLE
          IELIDN=IELTID(IGRUP)
          NNODE =IELPRP(3,IELIDN)
          NGAUSP=IELPRP(4,IELIDN)
C Extrapolate stresses and other state and algorithmic variables from
C gauss points to nodes
          IPOS=NGAUSP*NDIME+NGAUSP+1
          CALL EXTNOD
     1(   RELPRP(IPOS,IELIDN),
     2    STRSG(1,1,IELEM,1) ,STRSN    ,MSTRE     ,NGAUSP     ,NNODE   )
          CALL EXTNOD
     1(   RELPRP(IPOS,IELIDN),
     2    RSTAVA(1,1,IELEM,1),RSTAN    ,MRSTAV    ,NGAUSP     ,NNODE   )
          CALL EXTNOD
     1(   RELPRP(IPOS,IELIDN),
     2    RALGVA(1,1,IELEM,1),RALGN    ,MRALGV    ,NGAUSP     ,NNODE   )
C Nodal averaging
          DO INODE=1,NNODE
            IPOIN=IABS(LNODS(IELEM,INODE))
            NEWIPO=NEWNOD(IPOIN,IGRUP)
C Regular averaging (unweighted)
C            R1DVAL=R1/DBLE(NVALEN(IPOIN,IGRUP))
C Weighted averaging (by current element area/volume)
            R1DVAL=ELVOL(IELEM)/TELVOL(IPOIN,IGRUP)
C
            STRSAG(:,NEWIPO)=STRSAG(:,NEWIPO)+STRSN(:,INODE)*R1DVAL
            RSTAAG(:,NEWIPO)=RSTAAG(:,NEWIPO)+RSTAN(:,INODE)*R1DVAL
            RALGAG(:,NEWIPO)=RALGAG(:,NEWIPO)+RALGN(:,INODE)*R1DVAL
          ENDDO
        ENDDO
      ENDDO
      
C Write all stress components to VTK file
      DO ISTRE=1,NSTRE
        WRITE(IUNIT,'(/A,/A)')'SCALARS Stress_'//STRCMP(ISTRE)//
     1                   ' double','LOOKUP_TABLE default'
        WRITE(IUNIT,'(25E24.15e3)')STRSAG(ISTRE,:)
      ENDDO      
C
C Principal stresses
C ---------------------------------
C
      ALLOCATE(PSTRAG(3,NEWNPO))
      PSTRAG=R0
      DO INEWNP=1,NEWNPO
C 2-D analysis
        IF((NTYPE==1).OR.(NTYPE==2).OR.(NTYPE==3))THEN
          CALL PRINC2(PSTRS,STRSAG(1:3,INEWNP))
C 3-D analysis
        ELSEIF(NTYPE==4)THEN
          CALL PRINC3(PSTRS,STRSAG(1:6,INEWNP))
        ENDIF
        PSTRAG(:,INEWNP)=PSTRS
      ENDDO
C
C Write all stress components to VTK file
      DO I=1,3
        WRITE(IUNIT,'(/A,/A)')'SCALARS PrincipalStress_'//PSTRCM(I)//
     1                     ' double','LOOKUP_TABLE default'
        WRITE(IUNIT,'(25E24.15e3)')PSTRAG(I,:)
      ENDDO
C
C Von Mises effective stress and pressure (computed from nodal averaged
C stresses)
      ALLOCATE(PRESSG(NEWNPO))
      ALLOCATE(EFFSTG(NEWNPO))
      PRESSG=R0
      EFFSTG=R0
      
      IF(NTYPE==1)THEN
        PRESSG(:)=(STRSAG(1,:)+STRSAG(2,:))/R3
        EFFSTG(:)=SQRT(R3/R2*((STRSAG(1,:)-PRESSG(:))**2
     1                       +(STRSAG(2,:)-PRESSG(:))**2
     2                       +R2*STRSAG(3,:)**2+PRESSG(:)**2))
      ELSEIF((NTYPE==2).OR.(NTYPE==3))THEN
        PRESSG(:)=(STRSAG(1,:)+STRSAG(2,:)+STRSAG(4,:))/R3
        EFFSTG(:)=SQRT(R3/R2*((STRSAG(1,:)-PRESSG(:))**2
     1                       +(STRSAG(2,:)-PRESSG(:))**2
     2                       +R2*STRSAG(3,:)**2
     3                       +(STRSAG(4,:)-PRESSG(:))**2))
      ELSEIF(NTYPE==4)THEN
        PRESSG(:)=(STRSAG(1,:)+STRSAG(2,:)+STRSAG(3,:))/R3
        EFFSTG(:)=SQRT(R3/R2*((STRSAG(1,:)-PRESSG(:))**2
     1                       +(STRSAG(2,:)-PRESSG(:))**2
     2                       +(STRSAG(3,:)-PRESSG(:))**2
     3                       +R2*STRSAG(4,:)**2
     4                       +R2*STRSAG(5,:)**2
     5                       +R2*STRSAG(6,:)**2))
      ENDIF      
C
      WRITE(IUNIT,'(/A,/A)')'SCALARS Pressure double',
     1                      'LOOKUP_TABLE default'
      WRITE(IUNIT,'(25E24.15e3)')PRESSG(:)
      WRITE(IUNIT,'(/A,/A)')'SCALARS EffectiveStress double',
     1                      'LOOKUP_TABLE default'
      WRITE(IUNIT,'(25E24.15e3)')EFFSTG(:)
C
C Loop through material groups, find which ones are present in analysis
      ISMATP=.FALSE.
      DO IGRUP=1,NGRUP
        IMTTYP=IPROPS(1,MATTID(IGRUP))
        ISMATP(IMTTYP)=.TRUE.
      ENDDO
C
      write(*,*) SIZE(MYTMPVAR), NEWNPO
      !write(*,*) MYTMPVAR
      IF(ALLOCATED(MYTMPVAR)) DEALLOCATE(MYTMPVAR)
      !IF( .NOT. ALLOCATED( MYTMPVAR ) ) ALLOCATE(MYTMPVAR(NEWNPO))
      ALLOCATE(MYTMPVAR(NEWNPO))
C Loop through all material types
      DO IMTTYP=1,NMTTYP
C Skip those not present
        IF(.NOT.ISMATP(IMTTYP))CYCLE
C Loop through state variables of current material considered
        NVAR=NRSTVM(IMTTYP)
        DO IVAR=1,NVAR
C Skip variables with empty name
          IF(RSTACM(IVAR,IMTTYP)=='')CYCLE
C Copy values of nodal averaged variable
          MYTMPVAR=RSTAAG(IVAR,:)
          RVAL=DVLSTA(IVAR,IMTTYP)
C Loop through elements, finding those that do not 
          !DO IELEM=1,NELEM
          !  IGRUP=IGRPID(IELEM)
          !  JMTTYP=IPROPS(1,MATTID(IGRUP))
          !  IF(JMTTYP/=IMTTYP)THEN
          !    IELIDN=IELTID(IGRUP)
          !    NNODE=IELPRP(3,IELIDN)
          !    DO INODE=1,NNODE
          !      IPOIN=IABS(LNODS(IELEM,INODE))
          !      NEWIPO=NEWNOD(IPOIN,IGRUP)
          !      MYTMPVAR(NEWIPO)=RVAL
          !    ENDDO
          !  ENDIF
          !ENDDO
          WRITE(IUNIT,'(/A,/A)')'SCALARS '//TRIM(MATNAM(IMTTYP))//
     1                          '_RSTAV_'//TRIM(RSTACM(IVAR,IMTTYP))//
     2                          ' double','LOOKUP_TABLE default'
          WRITE(IUNIT,'(25E24.15e3)')MYTMPVAR
        ENDDO
C Loop through alg. variables of current material considered
        NVAR=NRALVM(IMTTYP)
        DO IVAR=1,NVAR
C Skip variables with empty name
          IF(RALGCM(IVAR,IMTTYP)=='')CYCLE
C Copy values of nodal averaged variable
          MYTMPVAR=RALGAG(IVAR,:)
          RVAL=DVLALG(IVAR,IMTTYP)
C Loop through elements, finding those that do not 
          DO IELEM=1,NELEM
            IGRUP=IGRPID(IELEM)
            JMTTYP=IPROPS(1,MATTID(IGRUP))
            IF(JMTTYP/=IMTTYP)THEN
              IELIDN=IELTID(IGRUP)
              NNODE=IELPRP(3,IELIDN)
              DO INODE=1,NNODE
                IPOIN=IABS(LNODS(IELEM,INODE))
                NEWIPO=NEWNOD(IPOIN,IGRUP)
                MYTMPVAR(NEWIPO)=RVAL
              ENDDO
            ENDIF
          ENDDO
          WRITE(IUNIT,'(/A,/A)')'SCALARS '//TRIM(MATNAM(IMTTYP))//
     1                          '_RALGV_'//TRIM(RALGCM(IVAR,IMTTYP))//
     2                          ' double','LOOKUP_TABLE default'
          WRITE(IUNIT,'(25E24.15e3)')MYTMPVAR
        ENDDO
      ENDDO
C
C Real state variables output
C      DO IRSTAV=1,NRSTAV
C        WRITE(IUNIT,'(/A,/A)')'SCALARS RSTAV_'//TRIM(RSTACM(IRSTAV))//
C     1                      ' double','LOOKUP_TABLE default'
C        WRITE(IUNIT,'(25E24.15e3)')RSTAAG(IRSTAV,:)
C      ENDDO
C
C Real algorithmic variables output
C      DO IRALGV=1,NRALGV
C        WRITE(IUNIT,'(/A,/A)')'SCALARS RALGV_'//TRIM(RALGCM(IRALGV))//
C     1                      ' double','LOOKUP_TABLE default'
C        WRITE(IUNIT,'(25E24.15e3)')RALGAG(IRALGV,:)
C      ENDDO
C
C 
C Hyplas node number
C ------------------
      WRITE(IUNIT,'(/A,/A)')'SCALARS NodeNo_HYPLAS integer',
     1                        'LOOKUP_TABLE default'
      WRITE(IUNIT,'(25(I0,1X))')OLDNOD(1:NEWNPO)   
C Cell data
C =========
      WRITE(IUNIT,'(/A,I0)')'CELL_DATA ',NELEM
C Material types
C --------------
      WRITE(IUNIT,'(/A,/A)')'SCALARS MaterialID int',
     1                      'LOOKUP_TABLE default'
      WRITE(IUNIT,'(25(I0,1X))')MATTID(IGRPID(1:NELEM))
C Element types
C -------------
      WRITE(IUNIT,'(/A,/A)')'SCALARS ElementID int',
     1                      'LOOKUP_TABLE default'
      WRITE(IUNIT,'(25(I0,1X))')IELTID(IGRPID(1:NELEM))
  
      CLOSE(UNIT=IUNIT,STATUS='KEEP')

      RETURN
      END
CDOC END_SUBROUTINE VTKOUT
