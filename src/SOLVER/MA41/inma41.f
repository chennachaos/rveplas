CDOC BEGIN_SUBROUTINE INMA41
CDOC Interface subroutine for MA41 linear system solver
CDOC	
CDOC Assembles the global stiffness matrix in sparse format and calls
CDOC the solver MA41 to obtain the solution to the corresponding linear 
CDOC system(s). More information about the solver, its input, error 
CDOC messaging, etc, can be found on the given manual file, 'ma41.pdf', 
CDOC or on http://www.hsl.rl.ac.uk.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DTIME  >  Time increment.
CDOC INTEGER          IFNEG  <  Signum (+1/-1) of the determinant of the
CDOC C                          stiffness matrix.
CDOC INTEGER          IITER  >  Current equilibrium iteration number.
CDOC LOGICAL          INCCUT <  Load increment cutting flag.
CDOC INTEGER          KRESL  >  Equation resolution index.
CDOC INTEGER          KUNLD  <> Unloading flag.
CDOC LOGICAL          UNSYM  >  Global unsymmetric solver flag.
CDOC
CDOC ACKNOWLEDGMENT
CDOC MA41 - Sparse multifrontal unsymmetric system solver, obtained from 
CDOC http://www.hsl.rl.ac.uk: "HSL(2011). A collection of Fortran codes 
CDOC for large scale scientific computation". The permission to use the 
CDOC code is gratefully acknowledged.
CDOC 
CDOC END_PARAMETERS
CHST
CHST F.M. Andrade Pires, February 2002: Initial coding as LEQSOL
CHST
CHST M.F. Adziman, D. de Bortoli, E.A. de Souza Neto, July 2013: 
CHST      Adjustments to make it compatible with HYPLAS v3.1:
CHST      - Updated ELEIST input parameters: added DTIME, removed IITER	
CHST      - Removed code specific to TRI3FB (F-bar patch elements):
CHST          loops over patch elements and calls to ESTIF using element
CHST          patch number
CHST      - Added master/slave degree of freedom treatment
CHST      - Added unloading flag reinitialization (as in FRONT)
CHST      Other improvements:
CHST      - Added check when using filter IFILT to avoid writing out of
CHST        array bounds
CHST      - MA41 error messages now handled by HYPLAS
CHST
      SUBROUTINE INMA41 ( DTIME, IFNEG, IITER, INCCUT, KRESL, KUNLD,
     1                    UNSYM )
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C 
      !INCLUDE 'mkl_pardiso.fi' ! MKL PARDISO interface headers
C
C Hyplas global database
C
      INCLUDE '../../MAXDIM.INC'
      INCLUDE '../../MATERIAL.INC'
      INCLUDE '../../ELEMENTS.INC'
      INCLUDE '../../GLBDBASE.INC'
C
C Arguments
C
      LOGICAL UNSYM, INCCUT
C
C Local arrays and variables
C
      DIMENSION ESTIF(MEVAB,MEVAB)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IFILT
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: RHS
      LOGICAL IS2D
C
C Sparse matrix representation (in coordinate/COO format)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IRN, JCN
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: ASPK
C
C Variables needed by MA41
      INTEGER :: KEEP(50), ICNTL(20), INFO(20)
      DOUBLE PRECISION :: CNTL(10), RINFO(20)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IW
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: RW, COLSCA, ROWSCA
C 
C Variables needed by MKL PARDISO
      !TYPE(MKL_PARDISO_HANDLE) PT(64) ! internal solver memory pointer
      INTEGER :: PT(64)
      INTEGER :: IPARM(64)
      INTEGER :: PHASE, ERROR
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: X
      INTEGER, ALLOCATABLE, DIMENSION(:) :: PERM
C Sparse matrix representation (in compressed sparse row/CSR format)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ICSR, JCSR
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: ACSR
C Variables needed for conversion from COO to CSR representation
      INTEGER JOBCSR(8)
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: ISDUPL
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ILEN
C     
      DOUBLE PRECISION, PARAMETER :: R0=0.0D0
C Enables printing solver information to the screen and input sparse
C matrix checking for MKL PARDISO
      LOGICAL, PARAMETER :: DEBUG=.FALSE.
      
      
      
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C #DEBUG-BEG#
      LOGICAL, PARAMETER :: USEMKL=.FALSE.
C #DEBUG-END#
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      
      
C
C
 1010 FORMAT(//' *** MA41 Error ***'/' Minimum required size of real' 
     1  ' working space:'/, I15)
 1020 FORMAT(//' *** MA41 Error ***'/' Minimum required size of integer' 
     1  ' working space:'/, I15)
 1030 FORMAT(//' *** MA41 Warning ***'/' Minimum recommended size of ' 
     1  ' real working space:'/, I15)
C
      IF((NTYPE==1).OR.(NTYPE==2).OR.(NTYPE==3))THEN
        IS2D=.TRUE.
      ELSEIF(NTYPE==4)THEN
        IS2D=.FALSE.
      ELSE
        STOP 'Invalid NTYPE in INMA41.f'
      ENDIF
C
C
C ==========================================
C Prepare input for the linear system solver
C ==========================================
C
C IFILT: filter array for the solver. For global dof I, IFILT(I) is the
C        corresponding dof in the reduced system that is solved (i.e.
C        without constrained and slave dofs)
C     N: order of the reduced linear system
C
      ALLOCATE(IFILT(NTOTV))
      IFILT=0
      N=0
      DO I=1,NTOTV
        IF(MASTER(I)==I.AND.IFFIX(I)==0)THEN !is master and not fixed
          N=N+1
          IFILT(I)=N
        ENDIF
      ENDDO
C
      ALLOCATE(RHS(N,2))
      RHS=R0
C Reset unloading flag
      IF(IITER.GT.1)KUNLD=0
C Initialise increment cutting flag
      INCCUT=.FALSE.
C Decide solution required
      IF(NALGO.LT.0.AND.IITER.EQ.1.AND.KRESL.EQ.2)THEN
        NRHS=0
        MODE=0
      ELSE IF(NALGO.LT.0.AND.IITER.EQ.1.AND.KRESL.EQ.1)THEN
        NRHS=1
        MODE=2
      ELSE IF(NALGO.LT.0.AND.IITER.GT.1.AND.KRESL.EQ.1)THEN
        NRHS=2
        MODE=3
      ELSE IF(NALGO.LT.0.AND.IITER.GT.1.AND.KRESL.EQ.2)THEN
        NRHS=1
        MODE=1
      ELSE
        NRHS=1
        MODE=1
      ENDIF
      IF(MODE.EQ.0) GOTO 1000
C
C Estimate number of entries on sparse stiffness matrix:
C each element's matrix has dimension (NEVAB,NEVAB), so it adds at most 
C NEVAB**2 entries (this is a simple approach that leads to lots of
C duplicate entries, i.e. with same line and column number)
      IESTNZ=0
      DO IELEM=1,NELEM
        NEVAB=IELPRP(5,IELTID(IGRPID(IELEM)))
        IESTNZ=IESTNZ+NEVAB**2
      ENDDO
C Allocate memory for COO format stiffness matrix
      ALLOCATE(IRN(IESTNZ))
      ALLOCATE(JCN(IESTNZ))
      ALLOCATE(ASPK(IESTNZ))
      ASPK=R0
      IRN=0
      JCN=0
C
C Loop over all elements, adding their stiffness matrix to the global
C sparse matrix
C ----------------------------------------------------------------------
      INNZSP=0
      DO 500 IELEM=1,NELEM
C Retrieve element properties
        IGRUP=IGRPID(IELEM)
        IELIDN=IELTID(IGRUP)
        NNODE=IELPRP(3,IELIDN)
        NEVAB=IELPRP(5,IELIDN)
C
C Call element interface routine for computation of element stiffness
C 
        CALL ELEIST (DTIME,ESTIF,IELEM,KUNLD,UNSYM )
C
C Transform the stiffness matrix into the local nodal coordinate system
C for prescribed displacements at an angle (for 2-D only)
C 
        IF(IS2D)THEN
          DO 100 INODE=1,NNODE
            LNODE=IABS(LNODS(IELEM,INODE))
            DO 80 IVFIX=1,NVFIX
              IF(NOFIX(IVFIX).EQ.LNODE)THEN
                IF(ANGLE(IVFIX).EQ.R0)GOTO 100
                C=COS(ANGLE(IVFIX))
                S=SIN(ANGLE(IVFIX))
                IEVAB=(INODE-1)*NDOFN+1
                JEVAB=IEVAB+1
                DO 60 KEVAB=1,NEVAB
                  GASHI= C*ESTIF(IEVAB,KEVAB)+S*ESTIF(JEVAB,KEVAB)
                  GASHJ=-S*ESTIF(IEVAB,KEVAB)+C*ESTIF(JEVAB,KEVAB)
                  ESTIF(IEVAB,KEVAB)=GASHI
                  ESTIF(JEVAB,KEVAB)=GASHJ
   60           CONTINUE
                DO 40 KEVAB=1,NEVAB
                  GASHI= ESTIF(KEVAB,IEVAB)*C+ESTIF(KEVAB,JEVAB)*S
                  GASHJ=-ESTIF(KEVAB,IEVAB)*S+ESTIF(KEVAB,JEVAB)*C
                  ESTIF(KEVAB,IEVAB)=GASHI
                  ESTIF(KEVAB,JEVAB)=GASHJ
   40           CONTINUE
                GOTO 100
              ENDIF
   80       CONTINUE
  100     CONTINUE
        ENDIF
C
C Assemble element loads (residual and/or external load) into global
C load vector 
C 
C ======================================================================
        DO 200 INODE=1,NNODE
          IEVAB=(INODE-1)*NDOFN
          LNODE=IABS(LNODS(IELEM,INODE))
          IF(IS2D)THEN
C   Use local nodal coordinate system if node has a prescribed
C   displacements at an angle (for 2-D only)
            DO 210 IVFIX=1,NVFIX
              IF(NOFIX(IVFIX).EQ.LNODE.AND.ANGLE(IVFIX).NE.R0)THEN
                C=COS(ANGLE(IVFIX))
                S=SIN(ANGLE(IVFIX))
                IEVAB=IEVAB+1
                JEVAB=IEVAB+1
C     Position of node's first degree of freedom in global vector
                IPOS=MASTER((LNODE-1)*NDOFN+1)
C     Position of node's first degree of freedom in reduced system (use 
C     filter previously built)
                ICPOS=IFILT(IPOS)
                IF(ICPOS.NE.0)THEN
                  IF(MODE.EQ.1)THEN
                    RHS(ICPOS,1)=RHS(ICPOS,1)+C*ELOAD(IEVAB,IELEM)
     1                                       +S*ELOAD(JEVAB,IELEM)
                  ELSE IF(MODE.EQ.2)THEN
                    RHS(ICPOS,1)=RHS(ICPOS,1)+C*RLOAD(IEVAB,IELEM)
     1                                       +S*RLOAD(JEVAB,IELEM)
                  ELSE IF(MODE.EQ.3)THEN
                    RHS(ICPOS,1)=RHS(ICPOS,1)+C*RLOAD(IEVAB,IELEM)
     1                                       +S*RLOAD(JEVAB,IELEM)
                    RHS(ICPOS,2)=RHS(ICPOS,2)+C*ELOAD(IEVAB,IELEM)
     1                                       +S*ELOAD(JEVAB,IELEM)
                  ENDIF
                ENDIF
C     Position of node's second degree of freedom in global vector
                JPOS=MASTER((LNODE-1)*NDOFN+2)
C     Position of node's second degree of freedom in reduced system (use 
C     filter previously built)
                JCPOS=IFILT(JPOS)
                IF(JCPOS.NE.0)THEN
                  IF(MODE.EQ.1)THEN
                    RHS(JCPOS,1)=RHS(JCPOS,1)+S*ELOAD(IEVAB,IELEM)
     1                                       +C*ELOAD(JEVAB,IELEM)
                  ELSE IF(MODE.EQ.2)THEN
                    RHS(JCPOS,1)=RHS(JCPOS,1)+S*RLOAD(IEVAB,IELEM)
     1                                       +C*RLOAD(JEVAB,IELEM)
                  ELSE IF(MODE.EQ.3)THEN
                    RHS(JCPOS,1)=RHS(JCPOS,1)+S*RLOAD(IEVAB,IELEM)
     1                                       +C*RLOAD(JEVAB,IELEM)
                    RHS(JCPOS,2)=RHS(JCPOS,2)+S*ELOAD(IEVAB,IELEM)
     1                                       +C*ELOAD(JEVAB,IELEM)
                  ENDIF
                ENDIF
                GOTO 200
              ENDIF
  210       CONTINUE
          ENDIF
C
          DO 140 IDOFN=1,NDOFN
            IEVAB=IEVAB+1
            IPOS=MASTER((LNODE-1)*NDOFN+IDOFN)
            ICPOS=IFILT(IPOS)
            IF(ICPOS.NE.0)THEN
              IF(MODE.EQ.1)THEN
                RHS(ICPOS,1)=RHS(ICPOS,1)+ELOAD(IEVAB,IELEM)
              ELSE IF(MODE.EQ.2)THEN
                RHS(ICPOS,1)=RHS(ICPOS,1)+RLOAD(IEVAB,IELEM)
              ELSE IF(MODE.EQ.3)THEN
                RHS(ICPOS,1)=RHS(ICPOS,1)+RLOAD(IEVAB,IELEM)
                RHS(ICPOS,2)=RHS(ICPOS,2)+ELOAD(IEVAB,IELEM)
              ENDIF
            ENDIF
  140     CONTINUE
  200   CONTINUE
C ======================================================================
C
C Assemble element global reduced stiffness matrix (in sparse format)
C 
        DO 300 INODE=1,NNODE
          NODEI=IABS(LNODS(IELEM,INODE))
          DO 220 IDOFN=1,NDOFN
            NROWG=MASTER((NODEI-1)*NDOFN+IDOFN)
            NROWE=(INODE-1)*NDOFN+IDOFN
            DO 240 JNODE=1,NNODE
              NODEJ=IABS(LNODS(IELEM,JNODE))
              DO 260 JDOFN=1,NDOFN
                NCOLG=MASTER((NODEJ-1)*NDOFN+JDOFN)
                NCOLE=(JNODE-1)*NDOFN+JDOFN
C   Add contribution of a prescribed displacement to right hand side
                IF(IFFIX(NCOLG).EQ.1)THEN
                  IF(IFFIX(NROWG).EQ.0)THEN
                    NROWC=IFILT(NROWG)
                    DO 270 IRHS=1,NRHS
                      RHS(NROWC,IRHS)=RHS(NROWC,IRHS)
     1                             -ESTIF(NROWE,NCOLE)*FIXED(NCOLG,IRHS)
  270               CONTINUE
                  ENDIF
                ENDIF
C   Add contribution of element stiffness to reduced global matrix
C   (only non-zero entries)
                IF((IFFIX(NROWG).EQ.0).AND.(IFFIX(NCOLG).EQ.0).AND.
     1             (ESTIF(NROWE,NCOLE).NE.R0))THEN
                  INNZSP=INNZSP+1
                  ASPK(INNZSP)=ESTIF(NROWE,NCOLE)
                  IRN(INNZSP)=IFILT(NROWG)
                  JCN(INNZSP)=IFILT(NCOLG)
                ENDIF
  260         CONTINUE
  240       CONTINUE
  220     CONTINUE
  300   CONTINUE
C
  500 CONTINUE
C
C NZ: number of non-zero entries in the sparse stiffness matrix
C
      NZ=INNZSP
      IF(NZ>IESTNZ)STOP 'Out of bounds for COO sparse matrix - INMA41.f'
      
      
      
      
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C #DEBUG-BEG#
      IF(USEMKL)THEN
C
C =========================================================
C Convert sparse matrix from COO (COOrdinate) format to CSR 
C (Compressed Sparse Row) format (as required by PARDISO)
C =========================================================
C
C Allocate memory for CSR matrix
      ALLOCATE(ACSR(NZ))
      ALLOCATE(JCSR(NZ))
      ALLOCATE(ICSR(N+1))
C
C Function MKL_DCSRCOO DOES NOT sum repeated entries in the COO sparse
C matrix (i.e. two entries with same IRN and JCN)
C Reference: software.intel.com/en-us/node/468628
C
      JOBCSR(1)=2 ! convert COO to CSR and sort column indices in CSR 
                  ! representation in increasing order within each row
      JOBCSR(2)=1 ! use one-based indexing for input COO matrix
      JOBCSR(3)=1 ! use one-based indexing for output CSR matrix
      JOBCSR(5)=NZ
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C #DEBUG-BEG#
      JOBCSR(6)=3 ! fill all output arrays (ACSR, JCSR, ICSR)
C! the correct value should be: JOBCSR(6)=0
C  but on this version of MKL (11.1.3), it doesn't work.
C  -> on current version (11.3), it works as expected
C #DEBUG-END#
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!      CALL MKL_DCSRCOO(JOBCSR, N, ACSR, JCSR, ICSR, NZ, ASPK, IRN, JCN,
!     1                 INFOCSR)
      IF(INFOCSR/=0)THEN
        WRITE(*,*)INFOCSR
        STOP 'Error in function MKL_DCSRCOO'
      ENDIF
C Release memory for COO matrix
      DEALLOCATE(ASPK)
      DEALLOCATE(IRN)
      DEALLOCATE(JCN)
C
C Flag duplicate entries from CSR matrix (they should be adjacent since
C MKL_DCSRCOO sorts the CSR matrix)
      ALLOCATE(ISDUPL(NZ))
      ISDUPL=.FALSE.
      ALLOCATE(ILEN(N))
C Find beginning, end and length of line I in CRS matrix
      DO I=1,N
        IRS=ICSR(I)
        IRE=ICSR(I+1)-1
        ILEN(I)=IRE-IRS+1
C If two adjacent entries in that line refer to same column, flag the 
C first one for removal and add its value to the second one
        DO J=IRS,IRE-1
          IF(JCSR(J)==JCSR(J+1))THEN
            ISDUPL(J)=.TRUE.
            ACSR(J+1)=ACSR(J+1)+ACSR(J)
            ILEN(I)=ILEN(I)-1 ! line length must be adjusted
          ENDIF
        ENDDO
      ENDDO
C
C Now rebuild ACSR and JCSR skipping duplicates
      OLDNZ=NZ
      NZ=0
      DO I=1,OLDNZ
        IF(ISDUPL(I))CYCLE
        NZ=NZ+1
        JCSR(NZ)=JCSR(I)
        ACSR(NZ)=ACSR(I)
      ENDDO     
C Set remaining elements of CSR matrix to zero
      JCSR(NZ+1:OLDNZ)=0
      ACSR(NZ+1:OLDNZ)=R0
C
C Finally, rebuild ICSR taking into account new line lengths
      ICSR(1)=1
      DO I=2,N+1
        ICSR(I)=ICSR(I-1)+ILEN(I-1)
      ENDDO
C
      DEALLOCATE(ISDUPL)
      DEALLOCATE(ILEN)
C
C ===============================================
C Call PARDISO to solve system of equations
C ===============================================
C Reference: software.intel.com/en-us/node/470282
C
      MAXFCT=1 ! Maximum number of factors with same sparsity structure to be kept in memory
      MNUM=1   ! Number of matrix to be factorised (1 ≤ mnum ≤ maxfct)
      ERROR=0  ! initialise error flag
      IF(DEBUG)THEN
        MSGLVL=1 ! 1: print statistical information
      ELSE
        MSGLVL=0 ! 0: no output from PARDISO
      ENDIF
C
      MTYPE=11 ! Matrix type:
               !  1 - real, structurally symmetric
               !  2 - real, symmetric positive definite
               ! -2 - real, symmetric indefinite
               ! 11 - real, nonsymmetric
C
C Set up PARDISO control parameters
C ----------------------------------
C Use default values for chosen MTYPE
      CALL PARDISOINIT(PT, MTYPE, IPARM)
      IPARM(6)=1 ! array RHS stores solution (x is still used internally)
      IPARM(8)=0 ! maximum number of iterative refinement steps
      IF(DEBUG)THEN
        IPARM(27)=1 ! check input sparse matrix for errors (e.g. columns not sorted)
      ELSE
        IPARM(27)=0
      ENDIF
C IPARM(4) - use CGS/CG (iterative methods) based on a factorisation
C previously performed
C IPARM(60) - out of core storage of factors (for problems that don't
C fit in RAM)
C
C PHASE controls the execution of the solver:
C A) Analysis, numerical factorisation, solve, iterative refinement
      ALLOCATE(X(N))
      ALLOCATE(PERM(N))
      PHASE=13
      CALL PARDISO(PT, MAXFCT, MNUM, MTYPE, PHASE, N, ACSR, ICSR, JCSR,
     1             PERM, 1, IPARM, MSGLVL, RHS(:,1), X, ERROR)
C Solve for second right-hand side in arc-length method
C (use same factorisation)
      IF(MODE==3)THEN
        PHASE=33
        CALL PARDISO(PT, MAXFCT, MNUM, MTYPE, PHASE, N, ACSR, ICSR,JCSR,
     1               PERM, 1, IPARM, MSGLVL, RHS(:,2), X, ERROR)
      ENDIF
      IF(ERROR/=0)THEN
        WRITE(*,*)ERROR
        STOP 'PARDISO returned the above error code (PHASE=13)'
      ENDIF
      IF(DEBUG)THEN
        WRITE(*,*)'Number of iterative refinement steps = ',IPARM(7)
      ENDIF
C B) Release internal memory
      PHASE=-1
      CALL PARDISO(PT, MAXFCT, MNUM, MTYPE, PHASE, N, ACSR, ICSR, JCSR,
     1             PERM, 1, IPARM, MSGLVL, RHS(1:1,1), X, ERROR)
      IF(ERROR/=0)THEN
        WRITE(*,*)ERROR
        STOP 'PARDISO returned the above error code (PHASE=-1)'
      ENDIF
      DEALLOCATE(X)
      DEALLOCATE(PERM)
      DEALLOCATE(ACSR)
      DEALLOCATE(JCSR)
      DEALLOCATE(ICSR)
C
      ELSE
C #DEBUG-END#
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
        
        
        
        
        
C
C ===================================
C Call MA41 to solve linear system(s)
C ===================================
C
C Initialization of control parameters for the solver
C
      CALL MA41ID (CNTL, ICNTL, KEEP)
      IF(DEBUG)THEN
        ICNTL(1)=6
        ICNTL(2)=6
        ICNTL(3)=6
        ICNTL(4)=3
      ELSE
C Do not print error, warning and statistics from MA41 execution (the 
C relevant ones will be handled by HYPLAS)
        ICNTL(1)=-1
        ICNTL(2)=-1
        ICNTL(3)=-1
        ICNTL(4)=-1
      ENDIF
C Control whether permutations are performed on the input matrix or not:
C     - if 0, no permutations are performed (default option)
C     - if 1, a maximum transversal algorithm is applied to the matrix
C       (recommended only for very unsymmetric matrices)
      ICNTL(6)=0
C Control the scaling strategy used:
C     - if 0, no scaling is performed (default option)
C     - if 6, use automatic scaling based on MC29 (recommended only for 
C       very badly scaled matrices - see MA41 manual)
      ICNTL(8)=0
C Set maximum number of steps of iterative refinement of solution
C (default is 0 - no refinement performed)
      ICNTL(10)=0
C Compute norm of input matrix and error estimates
      IF(DEBUG)THEN
        ICNTL(11)=1
      ENDIF
C 
      ALLOCATE(COLSCA(N))
      ALLOCATE(ROWSCA(N))
C
C 1) Analysis: symbolic factorisation (MA41 chooses pivot order)
C --------------------------------------------------------------
      JOB=1
C Estimated size of integer working space needed for analysis and
C factorisation (see MA41 manual)
      IIWSIZ=2*NZ+11*N+1
      ALLOCATE(IW(IIWSIZ))
      CALL MA41AD (JOB,N,NZ,IRN,JCN,ASPK,RHS(1,1),COLSCA,ROWSCA,KEEP,
     1             IW,IIWSIZ,DUMMY,IDUMMY,CNTL,ICNTL,INFO,RINFO)
C 2) Numerical factorisation
C --------------------------
      JOB=2
C After analysis, INFO(8) holds the real working space size required for
C the factorisation phase
      IRWSIZ=INFO(8)
      ALLOCATE(RW(IRWSIZ))
      CALL MA41AD (JOB,N,NZ,IRN,JCN,ASPK,RHS(1,1),COLSCA,ROWSCA,KEEP,
     1             IW,IIWSIZ,RW,IRWSIZ,CNTL,ICNTL,INFO,RINFO)
C 3) Solution using previously computed factorisation
C ---------------------------------------------------
      JOB=3
      CALL MA41AD (JOB,N,NZ,IRN,JCN,ASPK,RHS(1,1),COLSCA,ROWSCA,KEEP,
     1             IW,IIWSIZ,RW,IRWSIZ,CNTL,ICNTL,INFO,RINFO)
C Solve for second right-hand side in arc-length method
C (use same factorisation)
      IF(MODE==3)THEN
        CALL MA41AD (JOB,N,NZ,IRN,JCN,ASPK,RHS(1,2),COLSCA,ROWSCA,KEEP,
     1               IW,IIWSIZ,RW,IRWSIZ,CNTL,ICNTL,INFO,RINFO)
      ENDIF
C
C Relay error/warning messages from MA41, in case of execution errors
C -------------------------------------------------------------------
C
      IF(INFO(1).EQ.-5)THEN
        WRITE(*,1010)INFO(2)
        CALL ERRPRT('EE0303')
      ENDIF
C
      IF(INFO(1).EQ.-6)THEN
        CALL ERRPRT('EE0304')
      ENDIF
C
      IF(INFO(1).EQ.-7)THEN
        WRITE(*,1020)INFO(2)
        CALL ERRPRT('EE0305')
      ENDIF
C
      IF(INFO(1).EQ.-8)THEN
        WRITE(*,1020)INFO(2)
        CALL ERRPRT('EE0306')
      ENDIF
C
      IF(INFO(1).EQ.-9)THEN
        WRITE(*,1010)INFO(2)
        CALL ERRPRT('EE0307')
      ENDIF
C
      IF(INFO(1).EQ.-10)THEN
        CALL ERRPRT('EE0308')
      ENDIF
C
      IF(INFO(1).EQ.-11)THEN
        WRITE(*,1010)INFO(2)
        CALL ERRPRT('EE0309')
      ENDIF
C
      IF(INFO(1).EQ.-12)THEN
        WRITE(*,1010)INFO(2)
        CALL ERRPRT('EE0310')
      ENDIF
C
C     MA41 execution warning messages
C
      IF((INFO(1).EQ.1 ).OR.(INFO(1).EQ.3 ).OR.(INFO(1).EQ.5 ).OR.
     1   (INFO(1).EQ.7 ).OR.(INFO(1).EQ.9 ).OR.(INFO(1).EQ.11).OR.
     2   (INFO(1).EQ.13).OR.(INFO(1).EQ.15))THEN
        CALL ERRPRT('WE0301')
      ENDIF
C
      IF((INFO(1).EQ.2 ).OR.(INFO(1).EQ.3 ).OR.(INFO(1).EQ.6 ).OR.
     1   (INFO(1).EQ.7 ).OR.(INFO(1).EQ.10).OR.(INFO(1).EQ.11).OR.
     2   (INFO(1).EQ.14).OR.(INFO(1).EQ.15))THEN
        CALL ERRPRT('WE0302')
      ENDIF
C
      IF((INFO(1).EQ.4 ).OR.(INFO(1).EQ.5 ).OR.(INFO(1).EQ.6 ).OR.
     1   (INFO(1).EQ.7 ).OR.(INFO(1).EQ.12).OR.(INFO(1).EQ.13).OR.
     2   (INFO(1).EQ.14).OR.(INFO(1).EQ.15))THEN
        WRITE(*,1030)INFO(2)
        CALL ERRPRT('WE0303')
      ENDIF
C
      IF((INFO(1).EQ.8 ).OR.(INFO(1).EQ.9 ).OR.(INFO(1).EQ.10).OR.
     1   (INFO(1).EQ.11).OR.(INFO(1).EQ.12).OR.(INFO(1).EQ.13).OR.
     2   (INFO(1).EQ.14).OR.(INFO(1).EQ.15))THEN
        CALL ERRPRT('WE0304')
      ENDIF
      

      
      
      
      
      
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C #DEBUG-BEG#
      ENDIF
C #DEBUG-END#
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
      
      
      
      
      
      
      
C
C ============================================
C Post-processing: update arrays with solution
C ============================================
C
      IF(MODE.EQ.1)THEN
C If MODE=1 (standard system), the solution is already done and needs 
C only to be copied
        DO 600 ITOTV=1,NTOTV
          IF(IFILT(ITOTV).EQ.0)THEN
            DITER(ITOTV)=FIXED(ITOTV,1)
          ELSE 
            DITER(ITOTV)=RHS(IFILT(ITOTV),1)
          ENDIF
  600   CONTINUE
C
      ELSEIF(MODE.EQ.2)THEN
C If MODE=2 (tangential solution of the arc-length method), the solution
C is already done and needs only to be copied
        DO 700 ITOTV=1,NTOTV
          IF(IFILT(ITOTV).EQ.0)THEN
            DTANG(ITOTV)=FIXED(ITOTV,1)
          ELSE 
            DTANG(ITOTV)=RHS(IFILT(ITOTV),1)
          ENDIF
  700   CONTINUE
C
      ELSEIF(MODE.EQ.3)THEN
C For arc-tangent method, copy solutions from both right-hand sides
        DO 800 ITOTV=1,NTOTV
          IF(IFILT(ITOTV).EQ.0)THEN
            DTANG(ITOTV)=FIXED(ITOTV,1)
            DITER(ITOTV)=FIXED(ITOTV,2)
          ELSE 
            DTANG(ITOTV)=RHS(IFILT(ITOTV),1)
            DITER(ITOTV)=RHS(IFILT(ITOTV),2)
          ENDIF
  800   CONTINUE
      ENDIF
C
C Copy displacements of master degrees of freedom to slaves
C
      DO 810 ITOTV=1,NTOTV
        IF(MODE.EQ.1)THEN
          DITER(ITOTV)=DITER(MASTER(ITOTV))
        ELSE IF(MODE.EQ.2)THEN
          DTANG(ITOTV)=DTANG(MASTER(ITOTV))
        ELSE IF(MODE.EQ.3)THEN
          DTANG(ITOTV)=DTANG(MASTER(ITOTV))
          DITER(ITOTV)=DITER(MASTER(ITOTV))
        ENDIF
  810 CONTINUE
C
C Transform local nodal displacements back to global values for
C prescribed displacements at an angle (for 2-D only)
C
      IF(IS2D)THEN
        DO IVFIX=1,NVFIX
          IF(ANGLE(IVFIX)==R0)CYCLE
          IPOIN=NOFIX(IVFIX)
          ISVAB=MASTER((IPOIN-1)*NDOFN+1)
          JSVAB=MASTER((IPOIN-1)*NDOFN+2)
          C=COS(ANGLE(IVFIX))
          S=SIN(ANGLE(IVFIX))
          IF(MODE.EQ.1)THEN
            TEMPI=C*DITER(ISVAB)-S*DITER(JSVAB)
            TEMPJ=S*DITER(ISVAB)+C*DITER(JSVAB)
            DITER(ISVAB)=TEMPI
            DITER(JSVAB)=TEMPJ
          ELSE IF(MODE.EQ.2)THEN
            TEMPI=C*DTANG(ISVAB)-S*DTANG(JSVAB)
            TEMPJ=S*DTANG(ISVAB)+C*DTANG(JSVAB)
            DTANG(ISVAB)=TEMPI
            DTANG(JSVAB)=TEMPJ
          ELSE IF(MODE.EQ.3)THEN
            TEMPI=C*DTANG(ISVAB)-S*DTANG(JSVAB)
            TEMPJ=S*DTANG(ISVAB)+C*DTANG(JSVAB)
            DTANG(ISVAB)=TEMPI
            DTANG(JSVAB)=TEMPJ
            TEMPI=C*DITER(ISVAB)-S*DITER(JSVAB)
            TEMPJ=S*DITER(ISVAB)+C*DITER(JSVAB)
            DITER(ISVAB)=TEMPI
            DITER(JSVAB)=TEMPJ
          ENDIF
        ENDDO
      ENDIF
C
 1000 CONTINUE
C
      RETURN
      END
CDOC END_SUBROUTINE INMA41
