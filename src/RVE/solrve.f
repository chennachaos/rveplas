CDOC BEGIN_SUBROUTINE SOLRVE
CDOC Interface for the solution of the linear system of equations
CDOC associated to the IBVP of an RVE (computation of homogenised
CDOC stresses)
CDOC
CDOC This routine assembles the global stiffness matrix, deals with 
CDOC prescribed kinematical boundary constraints and solves the reduced 
CDOC linear system of equations using the solver MA41 from HSL
CDOC (parallel version of sparse Gaussian elimination)
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DTIME  >  Time increment.
CDOC INTEGER          KUNLD  <> Unloading flag.
CDOC DOUBLE_PRECISION RESVEC >  Reduced vector of global internal force 
CDOC                            specific for the chosen kinematical
CDOC                            constraint
CDOC LOGICAL          UNSYM  >  Stiffness matrix unsymmetry flag.
CDOC END_PARAMETERS
CHST
CHST       M.F. Adziman,July 2013: Initial coding of SOLRVE as coded in 
CHST                               MICROPLAST (Matlab multiscale code). 
CHST                               Uses MA41 solver.
CHST
CHST     D. de Bortoli,March 2014: Dimensioning parameters used for MA41
CHST                               arrays now come from MAXDIM.INC
CHST
CHST D. de Bortoli, December 2014: Stiffness matrix now assembled
CHST                               directly in sparse format for all
CHST                               multi-scale assumptions.
CHST
      SUBROUTINE SOLRVE
     1(   DTIME    ,KUNLD    ,UNSYM     ,INCCUT    ,RESVEC   )
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
!      INCLUDE 'mkl_pardiso.fi' ! MKL PARDISO interface headers
C
C Hyplas global database
      INCLUDE '../MAXDIM.INC'
      INCLUDE '../MATERIAL.INC'
      INCLUDE '../ELEMENTS.INC'
      INCLUDE '../GLBDBASE.INC'
      INCLUDE '../RVE.INC'
C
C Arguments
      LOGICAL UNSYM, INCCUT
      DOUBLE PRECISION RESVEC(MTOTV)
C
C Local arrays and variables
C
      DOUBLE PRECISION ESTIF(MEVAB,MEVAB)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: RHS
C
C Global sparse matrix representation (in coordinate/COO format)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IGRN, JGCN
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: GSPK
C Reduced sparse matrix representation (in coordinate/COO format)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IRRN, JRCN
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: RSPK
C Arrays used in construction of reduced sparse matrix for uniform 
C traction case
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: ISZ, ISDEP
      LOGICAL :: ISZKND, ISZKDN
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: KND, KDN, KDD
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
      
      
      
      DOUBLE PRECISION, PARAMETER :: R0=0.0D0
C
C***********************************************************************
C ASSEMBLE AND SOLVE REDUCED SYSTEM OF RVE IBVP LINEAR EQUATIONS
C***********************************************************************
C
 1010 FORMAT(//' *** MA41 Error ***'/' Minimum required size of real' 
     1  ' working space:'/, I15)
 1020 FORMAT(//' *** MA41 Error ***'/' Minimum required size of integer' 
     1  ' working space:'/, I15)
 1030 FORMAT(//' *** MA41 Warning ***'/' Minimum recommended size of ' 
     1  ' real working space:'/, I15)
C
C Check RVE kinematical choice is valid
      IF(IRV_RVEOPT/=1 .AND. IRV_RVEOPT/=2 .AND. IRV_RVEOPT/=3)THEN
        CALL ERRPRT('ED0402')
      ENDIF
C Stop program if solver is not MA41
      IF(NSOLVE.NE.2)CALL ERRPRT('ED0409')
C
C Initialise increment cutting flag
      INCCUT=.FALSE.
C Decide solution required
      IF(NALGO.LT.0)THEN
        STOP 'ERROR: Arc-length Method is not supported'
      ENDIF
C
C 1. Assemble full stiffness matrix in sparse format
C ==================================================
C
C Estimate number of entries on sparse stiffness matrix:
C each element's matrix has dimension (NEVAB,NEVAB), so it adds at most 
C NEVAB**2 entries (this is a simple approach that leads to lots of
C duplicate entries (same line and column number))
      IESTNZ=0
      DO IELEM=1,NELEM
        NEVAB=IELPRP(5,IELTID(IGRPID(IELEM)))
        IESTNZ=IESTNZ+NEVAB**2
      ENDDO
C Allocate memory for COO format global stiffness matrix
      ALLOCATE(IGRN(IESTNZ))
      ALLOCATE(JGCN(IESTNZ))
      ALLOCATE(GSPK(IESTNZ))
      GSPK=R0
      IGRN=0
      JGCN=0
C
C Add each element's contribution to global stiffness matrix
C
      NZG=0
      DO IELEM=1,NELEM
        NNODE=IELPRP(3,IELTID(IGRPID(IELEM)))
C Get element stiffness matrix
        CALL ELEIST (DTIME, ESTIF, IELEM, KUNLD, UNSYM)
C Add entries to corresponding global dof numbers
        DO INODE=1,NNODE
          NODEI=IABS(LNODS(IELEM,INODE))
          DO IDOFN=1,NDOFN
            IG=(NODEI-1)*NDOFN+IDOFN ! global row number
            IL=(INODE-1)*NDOFN+IDOFN ! local row number
            DO JNODE=1,NNODE
              NODEJ=IABS(LNODS(IELEM,JNODE))
              DO JDOFN=1,NDOFN
                JG=(NODEJ-1)*NDOFN+JDOFN ! global column number
                JL=(JNODE-1)*NDOFN+JDOFN ! local column number
                A=ESTIF(IL,JL)
                IF(A==R0)CYCLE ! zero entries don't need to go on sparse matrix
                NZG=NZG+1
                GSPK(NZG)=A
                IGRN(NZG)=IG
                JGCN(NZG)=JG
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C 2. Assemble reduced stiffness matrix in sparse format for the chosen
C    kinematic assumption
C ===================================================================
C    
C For the uniform traction case, some extra matrices need to be built
C beforehand and the estimated number of non-zero entries revised
C (unlike the other two cases, here the reduced stiffness matrix
C generally has MORE entries than the global one - because of all the
C relationships between dependent and free boundary dofs)
      IF(IRV_RVEOPT==3)THEN
C Number of interior, free and dependent dofs
        NDOFI=IRV_DFSPPT( 2)-IRV_DFSPPT( 1)+1
        NDOFF=IRV_DFSPPT(10)-IRV_DFSPPT( 9)+1
        NDOFD=IRV_DFSPPT(12)-IRV_DFSPPT(11)+1
C Array ISDEP indicates which dofs are dependent
        ALLOCATE(ISDEP(NTOTV))
        ISDEP=.FALSE.
        DO IDP=IRV_DFSPPT(11),IRV_DFSPPT(12)
          IDOF=IRV_DFSP(IDP)
          ISDEP(IDOF)=.TRUE.
        ENDDO
C Allocate dense matrices related to dependent dofs used in assembly of
C reduced sparse matrix:
C KND - relates non-dependent to dependent dofs
C KDN - relates dependent to non-dependent dofs
C KDD - relates dependent to dependent dofs
C 
        ALLOCATE(KND(NDOFI+NDOFF,NDOFD))
        ALLOCATE(KDN(NDOFD,NDOFI+NDOFF))
        ALLOCATE(KDD(NDOFD,NDOFD))
        KND=R0
        KDN=R0
        KDD=R0
C Populate them by scanning global stiffness matrix
        DO IG=1,NZG
          I=IGRN(IG)
          IR=IRV_IFILT(I)
          IF(IR==0)CYCLE ! ignore prescribed dofs
          J=JGCN(IG)
          JR=IRV_IFILT(J)
          IF(JR==0)CYCLE
          A=GSPK(IG)
          IF( (.NOT.ISDEP(I)).AND.ISDEP(J) )THEN
            KND(IR,JR)=KND(IR,JR)+A
          ELSEIF( ISDEP(I).AND.(.NOT.ISDEP(J)) )THEN
            KDN(IR,JR)=KDN(IR,JR)+A
          ELSEIF( ISDEP(I).AND.ISDEP(J) )THEN
            KDD(IR,JR)=KDD(IR,JR)+A
          ENDIF
        ENDDO
C
C Flag rows/columns of KND/KDN that are all zero-valued
        ALLOCATE(ISZ(NDOFI+NDOFF))
        ISZ=.FALSE.
        NNZKND=0
        NNZKDN=0
        DO I=1,NDOFI
          ISZKND=.TRUE.
          ISZKDN=.TRUE.
          IF( ANY( KND(I,:)/=R0 ) )THEN
            NNZKND=NNZKND+1
            ISZKND=.FALSE.
          ENDIF
          IF( ANY( KDN(:,I)/=R0 ) )THEN
            NNZKDN=NNZKDN+1
            ISZKDN=.FALSE.
          ENDIF
          ISZ(I)=ISZKND.AND.ISZKDN
        ENDDO
C Update estimate of number of non-zero entries in reduced stiffness
C matrix
        IESTNZ=IESTNZ+(NDOFI+NDOFF-COUNT(ISZ))*(2*NDOFF)
      ENDIF
C
C Allocate space for reduced stiffness matrix in COO format
      ALLOCATE(RSPK(IESTNZ))
      ALLOCATE(IRRN(IESTNZ))
      ALLOCATE(JRCN(IESTNZ))
      RSPK=R0
      IRRN=0
      JRCN=0
C NZ: number of non-zero entries in it
      NZ=0
C
C Linear and periodic kinematic assumptions
C -----------------------------------------
      IF(IRV_RVEOPT==1 .OR. IRV_RVEOPT==2)THEN
C
C Assemble if according to IRV_IFILT
        DO IG=1,NZG
          I=IGRN(IG)
          IR=IRV_IFILT(I)
          IF(IR==0)CYCLE ! skip dofs with IRV_IFILT=0
          J=JGCN(IG)
          JR=IRV_IFILT(J)
          IF(JR==0)CYCLE
          NZ=NZ+1
          IRRN(NZ)=IR
          JRCN(NZ)=JR
          RSPK(NZ)=GSPK(IG)
        ENDDO
C
C Uniform traction kinematic assumption
C -------------------------------------
      ELSEIF(IRV_RVEOPT==3)THEN
C Terms that do not involve dependent dofs do not need updating and
C are directly added
C K1 = [ Kii  Kif ]
C      [ Kfi  Kff ], where 'i' refers to interior dofs and 'f' to free
        DO IG=1,NZG
          I=IGRN(IG)
          IF(ISDEP(I))CYCLE ! skip dependent dofs (they are treated later)
          J=JGCN(IG)
          IF(ISDEP(J))CYCLE
          AVAL=GSPK(IG)
C
          IR=IRV_IFILT(I)
          IF(IR==0)CYCLE ! skip dofs with IRV_IFILT=0
          JR=IRV_IFILT(J)
          IF(JR==0)CYCLE
C
          NZ=NZ+1
          IRRN(NZ)=IR
          JRCN(NZ)=JR
          RSPK(NZ)=AVAL
        ENDDO 
C
C Now add terms relating dependent dofs to reduced stiffness matrix:
C K = K1 + K2 = [ Kii  Kif ] + [    0               Kid*R          ]
C               [ Kfi  Kff ]   [ R^T*Kdi   R^T*Kdf+Kfd*R+R^T*Kdd*R ],
C
C     where R is the dependency matrix relating dependent and free dofs,
C     'd' denote the dependent dofs and ^T denotes transposition
C
        IINS=IRV_DFSPPT(1)
        IFRS=IRV_DFSPPT(9)
C Loop through all non-dependent dofs (I) and free dofs (J)
        DO IND=1,NDOFI+NDOFF
          IF(ISZ(IND))CYCLE ! skip entries where both KND and KDN are only zero
          DO JFR=1,NDOFF
            AIJ=SUM( KND(IND,:)*DRV_BNDR(1:NDOFD,JFR) ) ! Kid*R, Kfd*R terms
            AJI=SUM( KDN(:,IND)*DRV_BNDR(1:NDOFD,JFR) ) ! R^T*Kdi, R^T*Kdf terms
C Find number of dof I in global system (IG) and reduced system (IR)...
C ...if IND<=NDOFI, dof I is an interior dof
            IF(IND<=NDOFI)THEN
              IG=IRV_DFSP(IINS+IND-1)
C ...else it is a free dof, so add R^T*Kdd*R term
            ELSE
              IFR=IND-NDOFI
              IG=IRV_DFSP(IFRS+IFR-1)
              DO KD=1,NDOFD
                DO LD=1,NDOFD
                  AIJ=AIJ+DRV_BNDR(KD,IFR)*KDD(KD,LD)*DRV_BNDR(LD,JFR)
                ENDDO
              ENDDO
            ENDIF
            IR=IRV_IFILT(IG)
C Find number of dof J in global system (JG) and reduced system (JR)
            JG=IRV_DFSP(IFRS+JFR-1)
            JR=IRV_IFILT(JG)
C
C Add Kid*R, Kfd*R and R^T*Kdd*R terms
            IF(AIJ/=R0)THEN
              NZ=NZ+1
              IRRN(NZ)=IR
              JRCN(NZ)=JR
              RSPK(NZ)=AIJ
            ENDIF
C
C Add R^T*Kdi and R^T*Kdf terms
            IF(AJI/=R0)THEN
              NZ=NZ+1
              IRRN(NZ)=JR ! swapped indices!
              JRCN(NZ)=IR
              RSPK(NZ)=AJI
            ENDIF
          ENDDO
        ENDDO
C       
        DEALLOCATE(ISDEP)
        DEALLOCATE(ISZ)
        DEALLOCATE(KND)
        DEALLOCATE(KDN)
        DEALLOCATE(KDD)
C
      ENDIF
C
C Global stiffness matrix not necessary any more
      DEALLOCATE(GSPK)
      DEALLOCATE(IGRN)
      DEALLOCATE(JGCN)
C
C Checking if allocated memory was enough
      IF(NZ>IESTNZ)STOP 'SOLRVE - RSPK out of bounds!'!CALL ERRPRT('EE0302')
C
C 3. Solve linear system
C ======================
C
C N: number of degrees of freedom in reduced system
      N=MAXVAL(IRV_IFILT(1:NTOTV))
C Assign residual vector to right hand side
      ALLOCATE(RHS(N))
      RHS=RESVEC(1:N)
C
      
      
      WRITE(*,*) "USEMKL", USEMKL
      
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
C matrix (i.e. two entries with same IRRN and JGCN)
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

!      CALL MKL_DCSRCOO(JOBCSR, N, ACSR, JCSR, ICSR, NZ, RSPK, IRRN,JRCN,
!     1                 INFOCSR)
      IF(INFOCSR/=0)THEN
        WRITE(*,*)INFOCSR
        STOP 'Error in function MKL_DCSRCOO'
      ENDIF
C Release COO matrix memory
      DEALLOCATE(RSPK)
      DEALLOCATE(IRRN)
      DEALLOCATE(JRCN)
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
     1             PERM, 1, IPARM, MSGLVL, RHS, X, ERROR)
C      IF(ERROR/=0)THEN
C        WRITE(*,*)ERROR
C        STOP 'PARDISO returned the above error code (PHASE=13)'
C      ENDIF
      IF(DEBUG)THEN
        WRITE(*,*)'Number of iterative refinement steps = ',IPARM(7)
      ENDIF
C B) Release internal memory
      PHASE=-1
      CALL PARDISO(PT, MAXFCT, MNUM, MTYPE, PHASE, N, ACSR, ICSR, JCSR,
     1             PERM, 1, IPARM, MSGLVL, RHS, X, ERROR)
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

        
        
        
        
        
        
C Solve using MA41
C-----------------
C Initialization of control parameters for the solver
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
      CALL MA41AD (JOB,N,NZ,IRRN,JRCN,RSPK,RHS,COLSCA,ROWSCA,KEEP,
     1             IW,IIWSIZ,DUMMY,IDUMMY,CNTL,ICNTL,INFO,RINFO)
C 2) Numerical factorisation
C --------------------------
      JOB=2
C After analysis, INFO(8) holds the real working space size required for
C the factorisation phase
      IRWSIZ=INFO(8)
      ALLOCATE(RW(IRWSIZ))
      CALL MA41AD (JOB,N,NZ,IRRN,JRCN,RSPK,RHS,COLSCA,ROWSCA,KEEP,
     1             IW,IIWSIZ,RW,IRWSIZ,CNTL,ICNTL,INFO,RINFO)
C 3) Solution using previously computed factorisation
C ---------------------------------------------------
      !DO I=1,NZ
        !WRITE(*,*) I, IRRN(I), JRCN(I), RSPK(I)
      !ENDDO

      !DO I=1,N
        !WRITE(*,*) I, RHS(I)
      !ENDDO
      
      
      JOB=3
      CALL MA41AD (JOB,N,NZ,IRRN,JRCN,RSPK,RHS,COLSCA,ROWSCA,KEEP,
     1             IW,IIWSIZ,RW,IRWSIZ,CNTL,ICNTL,INFO,RINFO)
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
C 4. Update iterative global nodal displacements vector
C =====================================================
C
      DITER=R0
C Do update according to IRV_IFILT
      DO I=1,NTOTV
        IPOS=IRV_IFILT(I)
        IF(IPOS==0)CYCLE
        DITER(I)=-RHS(IRV_IFILT(I))
      ENDDO
C      
C For uniform traction case, calculate the displacement of dependent
C dofs (Ud = R*Uf)
      IF(IRV_RVEOPT.EQ.3)THEN
        DO IDP=IRV_DFSPPT(11),IRV_DFSPPT(12)
          IPOS=IRV_DFSP(IDP)
          IR=IRV_IFILT(IPOS)
          DITER(IPOS)=R0
          JTEMP=0
          DO JFR=IRV_DFSPPT(9),IRV_DFSPPT(10)
            JPOS=IRV_DFSP(JFR)
            JR=IRV_IFILT(JPOS)
            JTEMP=JTEMP+1
            DITER(IPOS)=DITER(IPOS) - RHS(JR)*DRV_BNDR(IR,JTEMP)
          ENDDO
        ENDDO
      ENDIF
C
      RETURN
      END
CDOC END_SUBROUTINE SOLRVE
