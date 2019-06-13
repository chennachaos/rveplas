C***********************************************************************
C
C 1. RVE data base   ( COMMON blocks )
C
C***********************************************************************
C
C   RVE related main variables
C   --------------------------------------------
C
C
      COMMON/RVE/
     1   IRV_NDSP(MPOIN)     ,IRV_DFSP(MTOTV)         ,IRV_NDSPPT(8)  ,
     2   IRV_DFSPPT(14)      ,DRV_BNDR(9,MTOTV)       ,DRV_STRAIN(9)  ,
     3   DRV_CELVOL          ,IRV_RVEOPT              ,IRV_IFILT(MTOTV)
C            
C
C       - IRV_NDSP   = An array of node numbers resulting from the mesh 
C                      splitting routine. The order of this array
C                      depends on the kinematical constraint chosen. 
C                      It is acessed using the pointer indices stored
C                      in IRV_NDSPPT (refer to its description below for 
C                      more information). 
C
C       - IRV_DFSP   = An array of d.o.f. numbers resulting from the 
C                      mesh splitting routine. The order of this array
C                      depends on the kinematical constraint chosen. 
C                      It is acessed using the pointer indices stored
C                      in IRV_DFSPPT (refer to its description below for 
C                      more information). 
C
C       - IRV_NDSPPT = Pointers for the array IRV_NDSP. Sequence:
C
C                       Node group   |  First index   |   Last index   |
C                  ––––––––––––––––––+––––––––––––––––+––––––––––––––––+
C                    Interior nodes  | IRV_NDSPPT(1)  | IRV_NDSPPT(2)  |
C                        Plus nodes  | IRV_NDSPPT(3)  | IRV_NDSPPT(4)  |
C                       Minus nodes  | IRV_NDSPPT(5)  | IRV_NDSPPT(6)  |
C                      Corner nodes  | IRV_NDSPPT(7)  | IRV_NDSPPT(8)  |
C
C       - IRV_DFSPPT = Pointers for the array IRV_DFSP. Sequence:
C                      
C                        Dof group   |  First index   |   Last index   |
C                  ––––––––––––––––––+––––––––––––––––+––––––––––––––––+
C                     Interior dofs  | IRV_DFSPPT(1)  | IRV_DFSPPT(2)  |
C                         Plus dofs  | IRV_DFSPPT(3)  | IRV_DFSPPT(4)  |
C                        Minus dofs  | IRV_DFSPPT(5)  | IRV_DFSPPT(6)  |
C                       Corner dofs  | IRV_DFSPPT(7)  | IRV_DFSPPT(8)  |
C                         Free dofs  | IRV_DFSPPT(9)  | IRV_DFSPPT(10) |
C                    Dependent dofs  | IRV_DFSPPT(11) | IRV_DFSPPT(12) |
C                   Prescribed dofs  | IRV_DFSPPT(13) | IRV_DFSPPT(14) |
C
C       - DRV_BNDR   = Global dependency matrix for the uniform boundary
C                      traction assumption (integral constraint on
C                      external boundary edges/faces).
C
C
C
C	    - DRV_STRAIN = In small strains, it is the array representation
C                      of the prescribed (macroscopic) infinitesimal
C                      strain tensor:
C
C                      2-D problems:  3 components
C                        DRV_STRAIN(1): epsilon_xx
C                        DRV_STRAIN(2): epsilon_yy
C                        DRV_STRAIN(3): epsilon_xy
C
C                      3-D problems:  6 components
C                        DRV_STRAIN(1): epsilon_xx
C                        DRV_STRAIN(2): epsilon_yy
C                        DRV_STRAIN(3): epsilon_zz
C                        DRV_STRAIN(4): epsilon_xy
C                        DRV_STRAIN(5): epsilon_yz
C                        DRV_STRAIN(6): epsilon_xz
C
C                      In large strains, it is the array representation
C                      of the prescribed (macroscopic) deformation
C                      gradient.
C
C                      2-D problems:  4 components
C                        DRV_STRAIN(1): Fxx
C                        DRV_STRAIN(2): Fxy
C                        DRV_STRAIN(3): Fyx
C                        DRV_STRAIN(4): Fyy
C
C                      3-D problems:  9 components
C                        DRV_STRAIN(1): Fxx
C                        DRV_STRAIN(2): Fxy
C                        DRV_STRAIN(3): Fxz
C                        DRV_STRAIN(4): Fyx
C                        DRV_STRAIN(5): Fyy
C                        DRV_STRAIN(6): Fyz
C                        DRV_STRAIN(7): Fzx
C                        DRV_STRAIN(8): Fzy
C                        DRV_STRAIN(9): Fzz
C
C       - DRV_CELVOL = Volume of micro-cell.
C
C       - IRV_RVEOPT = RVE kinematical constraint option:
C                      1 - Linear boundary displacement fluctuations
C                      2 - Periodic displacement fluctuations
C                      3 - Uniform boundary traction
C
C       - IRV_IFILT  = relates a dof's global number to its number in 
C                      the reduced linear system:
C                      IRD=IRV_IFILT(IGL), for global dof number IGL
C                      and its corresponding reduced system number IRD.
C                      This array is constructed and used differently
C                      for different kinematical assumptions.
C                      Dofs that are not part of the reduced linear
C                      system (e.g. external boundary dofs in the
C                      linear displacement assumption, or corner
C                      dofs in the periodic displacement fluctuations
C                      assumption) have IRV_IFILT=0
C
C-----------------------------------------------------------------------
C
C   RVE related submain variables (for programming convenience)
C   -----------------------------------------------------------
C
C      COMMON/RVE_SUBMAIN/
C     1   IDFSP(MTOTV)  ,IPDFIN(2)     ,IPDFPL(2)     ,IPDFMI(2)        ,
C     2   IPDFCO(2)     ,IPDFFR(2)     ,IPDFDP(2)     ,RBOUND(9,MTOTV)  ,
C     3   INDSP(MPOIN)  ,IPNDIN(2)     ,IPNDPL(2)     ,IPNDMI(2)        ,
C     4   IPNDCO(2)    
C
C       -  IDFSP  = Global dofs pointer vector (equals to IRV_DFSP)
C
C       -  IPDFIN = Interior dofs pointer vector 
C
C       -  IPDFPL = 'plus' boundary dofs pointer vector
C
C       -  IPDFMI = 'minus' boundary dofs pointer vector
C
C       -  IPDFCO = Corner dofs pointer vector
C
C       -  IPDFFR = 'free' boundary dofs pointer vector
C
C       -  IPDFDP = 'dependent' boundary dofs pointer vector 
C
C       -  RBOUND = Integral constraint matrix for 'uniform traction' 
C
C       -  INDSP  = Global nodes pointer vector
C
C       -  IPNDIN = Interior nodes pointer vector
C
C       -  IPNDPL = 'plus' boundary nodes pointer vector
C
C       -  IPNDMI = 'minus' boundary nodes pointer vector
C
C       -  IPNDCO = Corner nodes pointer vector
C
C
C      COMMON/RVE_TEMP/
C     5    CGLOB_TEMP(9,MTOTV)
C
C      PARAMETER(   MDEPDF=4   )
C       - MDEPDF     = Maximum number of dependent dofs. In large 
C                      strains, 4 are needed in 2-D, while 9 are 
C                      necessary in 3-D
C
C      PARAMETER( DRV_GEOTOL  = 1.0D-2 )
C
C       - DRV_GEOTOL = Relative tolerance used for geometric operations
C                      on RVE geometry (splitting of the mesh in
C                      boundary and interior nodes and associated
C                      procedures)