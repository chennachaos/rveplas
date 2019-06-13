C***********************************************************************
C----------------------------------------------------------------------*
C                                                                      *
C                                                                      *
C      H Y P L A S  MAXIMUM PROBLEM SIZE DIMENSIONING PARAMETERS       *
C                                                                      *
C                                                                      *
C----------------------------------------------------------------------*
C***********************************************************************
C
C
C
C  Dimensioning parameters:
C  -----------------------
C
C   * MELEM = Maximum permissible number of elements in the mesh.
C
C   * MFRON = Maximum front width allowed in frontal solver.
C
C   * MGRUP = Maximum permissible number of element groups.
C
C   * MINCS = Maximum permissible number of load increments.
C
C   * MPOIN = Maximum permissible number of nodal points in the mesh.
C
C   * MSUBIN= Dimension of increment cut stack array. Maximum
C             permissible number of consecutive increment cuts is
C             MSUBIN-1.
C
C   * MVFIX = Maximum permissible number of nodes with prescribed
C             displacements. Does not include (slave) nodes with
C             kinematical contraints prescribed by means of Master/Slave
C             constraint specification.
C
C   * MDOFGL= Maximum permissible number of degrees of freedom in
C             reduced global stiffness matrix (only used by solvers that
C             assemble the global stiffness matrix, e.g. MA41)
C
C   * MNNZSP= Maximum non zero entries in the sparse reduced global 
C             stiffness matrix (only used by sparse solvers, e.g. MA41)
C
C   * MA41RW= Maximum size of the real working space used by the solver 
C             MA41
C
C   * MA41IW= Maximum size of the integer working space used by the 
C             solver MA41
C
      PARAMETER( MELEM  = 1000)
      PARAMETER( MFRON  = 500)
      PARAMETER( MGRUP  = 10)
      PARAMETER( MINCS  = 1000)
      PARAMETER( MPOIN  = 1000)
      PARAMETER( MSUBIN = 100  )
      PARAMETER( MVFIX  = 500 )
