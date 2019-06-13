CDOC BEGIN_PROGRAM HYPLAS
CDOC This is HYPLAS main program
CDOC
C***********************************************************************
C----------------------------------------------------------------------*
C|                                                                    |*
C|                                                                    |*
C|    H Y P L A S    version 4.0.1                                    |*
C|                                                                    |*
C|                                                                    |*
C|    Program for implicit small and large strain                     |*
C|    finite element analysis of hyperelastic and                     |*
C|    visco-elastoplastic solids.                                     |*
C|    Includes multi-scale capabilities: computation                  |*
C|    of homogenised stresses resulting from a                        |*
C|    Representative Volume Element (RVE)being                        |*
C|    subjected to a given strain path.                               |*
C|                                                                    |*
C|                                                                    |*
C|    Copyright (c) 1996-2015  EA de Souza Neto, D Peric & DRJ Owen   |*
C|                             Civil and Computational Eng. Centre    |*
C|                             School of Engineering                  |*
C|                             Swansea University                     |*
C|                                                                    |*
C|                                                                    |*
C|    This program is a companion to the textbook:                    |*
C|    EA de Souza Neto, D Peric & DRJ Owen. Computational             |*
C|    Methods for Plasticity: Theory and Applications. Wiley,         |*
C|    Chichester, 2008.                                               |*
C|                                                                    |*
C|                                                                    |*
C|    Please send BUG REPORTS to                                      |*
C|                                                                    |*
C|                        hyplas_v2.0@live.co.uk                      |*
C|                                                                    |*
C|    NOTE: Messages sent to the authors' personal email addresses    |*
C|          will NOT be answered.                                     |*
C|                                                                    |*
C----------------------------------------------------------------------*
C***********************************************************************
C----------------------------------------------------------------------*
C                                                                      *
C     COPYRIGHT STATEMENT                                              *
C                                                                      *
C     You may only use this program for your own private purposes.     *
C     You are not allowed, in any circumstances, to distribute this    *
C     program (including its source code, executable and any other     *
C     files related to it, either in their original version or any     *
C     modifications introduced by you, the authors or any other party) *
C     in whole or in part, either freely or otherwise, in any medium,  *
C     without the prior written consent of the copyright holders.      *
C                                                                      *
C     DISCLAIMER                                                       *
C                                                                      *
C     This program (including its source code, executable and any      *
C     other files related to it) is provided "as is" without warranty  *
C     of any kind, either expressed or implied, including, but not     *
C     limited to, any implied warranties of fitness for purpose.       *
C     In particular, THIS PROGRAM IS BY NO MEANS GUARANTEED TO BE FREE *
C     FROM ERRORS.                                                     *
C     The results produced by this program are in no way garanteed to  *
C     be fit for any purpose.                                          *
C     This program (or any modification incorporated to it by you, the *
C     authors or any other party) will be used entirely at your own    *
C     risk.                                                            *
C     Under no circumstances will the authors/copyright holders be     *
C     liable to anyone for damages, including any general, special,    *
C     incidental or consequential damages arising from the use or      *
C     inability to use the program (including, but not limited to,     *
C     loss or corruption of data, failure of the program to operate in *
C     any particular way as well as damages arising from the use of    *
C     any results produced by the program for any purpose).            *
C                                                                      *
C     CONDITIONS OF USE                                                *
C                                                                      *
C     You may only use this program if you fully understand and agree  *
C     with the terms of the above disclaimer. You must not use this    *
C     program if you do not agree with or do not understand (fully or  *
C     in part) these conditions of use.                                *
C                                                                      *
C----------------------------------------------------------------------*
C***********************************************************************
C
C
      PROGRAM HYPLAS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C Hyplas database: Global parameters and common blocks
      INCLUDE 'MAXDIM.INC'
      INCLUDE 'MATERIAL.INC'
      INCLUDE 'ELEMENTS.INC'
      INCLUDE 'GLBDBASE.INC'
      INCLUDE 'RVE.INC'
C Common block of arrays used only by the frontal solver
      COMMON / FRONTA /
     1    EQRHS(MTOTV,2)     ,EQROW(MFRON,MTOTV) ,EQCOL(MFRON,MTOTV) ,
     2    DECAY(MFRON)       ,GLOAD(MFRON,2)     ,VECRV(MFRON,2)     ,
     3    LOCEL(MEVAB,MELEM) ,NACVA(MFRON,MELEM) ,NAMEV(MTOTV)       ,
     4    NDEST(MEVAB,MELEM) ,NPIVO(MTOTV)       ,NFRON(MELEM)
C Logical control flags for main program
      LOGICAL
     1    CONVRG     ,DIVERG     ,INCCUT     ,RSTRT      ,UNSYM
C File names
      CHARACTER*256
     1    DATFIL     ,RESFIL     ,RSTINP     ,RSTOUT     ,PLOTD     ,
     2    DATFIL2
C Increment control arrays for main program
      DIMENSION
     1    DFACTV(MINCS)      ,DFSUB(MSUBIN)      ,DTIMEV(MINCS)     ,
     2    DTSUB(MSUBIN)      ,MITERV(MINCS)      ,NOUTP(5)           ,
     3    NOUTPV(5,MINCS)    ,TOLERV(MINCS)
C Local arrays for micro-to-macro transitions (homogenisation)
      DIMENSION
     1    DUTYLR(MTOTV) ,RESVEC(MTOTV) ,FBD(MTOTV)
C Numerical constants
      PARAMETER
     1(   R0=0.0D0   ,RP5=0.5D0  ,RP7=0.7D0  )
C***********************************************************************
C
 1000 FORMAT(//20X,' H Y P L A S   ANALYSIS  RESULTS'/
     1         20X,'================================='/)
 1010 FORMAT(
     1 4X,' __________________________________________________________'
     1                                                            '__ '/
     2 4X,'|                                                          '
     2                                                            '  |'/
     3 4X,'|   Program compiled with the dimensioning parameters      '
     3                                                            '  |'/
     4 4X,'|__________________________________________________________'
     4                                                            '__|'/
     5 4X,'|                                                 |        '
     5                                                            '  |'/
     6 4X,'| Maximum number of elements              (MELEM) |',I9,' |'/
     7 4X,'| Maximum frontwidth allowed in solution  (MFRON) |',I9,' |'/
     8 4X,'| Maximum number of element groups        (MGRUP) |',I9,' |'/
     9 4X,'| Maximum number of load increments       (MINCS) |',I9,' |'/
     O 4X,'| Maximum number of nodal points          (MPOIN) |',I9,' |'/
     1 4X,'| Size of increment cutting stack array   (MSUBIN)|',I9,' |'/
     2 4X,'| Max. number of nodes with prescr.displ. (MVFIV) |',I9,' |'/
     7 4X,'|_________________________________________________|________'
     7                                                            '__|'/
     8                                                                 )
 1015 FORMAT(/,' Data file name:'/
     1         ' ==============='//,1X,A/)
 1020 FORMAT(////15X,'Reading data...')
 1030 FORMAT(////15X,'H Y P L A S   ANALYSIS starting...'/)
 1050 FORMAT(////
     1' INCREMENT NUMBER',I5,6X,'TIME = ',G12.6,2X,
     2 'LOAD FACTOR = ',G12.6/
     2' --------------------------------------------------------------',
     3'------------'/
     4 4X,'         ',13X,'relative residual',13X,'maximum residual'/
     5 4X,'iteration',13X,'    norm (%)     ',13X,'     norm       '/
     6' --------------------------------------------------------------',
     7'------------')
 1055 FORMAT(////
     1' INCREMENT NUMBER',I5,19X,'       ARC LENGTH =',G15.6/
     2' --------------------------------------------------------------',
     3'------------'/
     4 4X,'         ',6X,'relative residual',4X,'maximum residual',
     5 5X,'   total'/
     6 4X,'iteration',6X,'    norm (%)     ',4X,'     norm       ',
     7 5X,'load factor'/
     8' --------------------------------------------------------------',
     9'------------')
 1060 FORMAT(
     1' --------------------------------------------------------------',
     2'------------')
 1063 FORMAT(5X,'TIME INCREMENT = ',G12.6,3X,
     1           'INCREMENTAL LOAD FACTOR = ',G12.6)
 1065 FORMAT(30X,'CONVERGED TOTAL LOAD FACTOR =',G15.6)
 1067 FORMAT(24X,'CONVERGED INCREMENTAL LOAD FACTOR =',G15.6)
 1070 FORMAT(/////15X,'Program  H Y P L A S  successfully completed.')
 1080 FORMAT(///15X,'Data file name was ---------> ',A)
 1090 FORMAT(// 15X,'Results file name is -------> ',A)
 1092 FORMAT(// 15X,'Kinematical constraint is --> ',A///) 
 1095 FORMAT(   15X,'Re-start file name is ---> ',A//
     1          15X,'Last increment written -->',I5///)
 1040 FORMAT(//' Iterations not converged.')
 1100 FORMAT(//' Iterations diverging.')
 1110 FORMAT(/ ' Re-trying with reduced increment size...'/)
 1120 FORMAT(/ ' Re-trying with reduced arc length...'/)
C
C
C
C Start up. Read data, initialise variables, etc...
C
C REFERENCE: Flowchart of Figure 5.1
C
C *************************************************
C
C Send greeting message to standard output
       CALL GREET
C Read names and open relevant files
      CALL FOPEN(  DATFIL     ,RESFIL     ,RSTOUT     ,NMULTI,PLOTD)
C      OPEN(UNIT=13,FILE=PLOTD.TXT,STATUS='OLD')
C Echo dimensioning parameters defined in file MAXDIM.INC
c       WRITE(16,1000)
c       WRITE(16,1010)MELEM,MFRON,MGRUP,MINCS,MPOIN,MSUBIN,MVFIX
C Echo data file name
       I=INDEX(DATFIL,' ')-1
c       WRITE(16,1015)DATFIL(1:I)
C
C Read relevant data from input data/re-start file
C ------------------------------------------------
C Check if main data is to be read from the input data file or from an
C input re-start file
      CALL RSTCHK(  RSTINP     ,RSTRT      )
c       WRITE(*,1020)
C
      IF(RSTRT)THEN
C Re-start mode: Read main data from input re-start file
         CALL RSTART
     1(  DFOLD      ,DLENG      ,DLENGO     ,DLENM      ,DLAMD      ,
     2   IFNEG      ,IINCS      ,MXFRON     ,NOUTP      ,TFACT      ,
     3   TFACTO     ,TTIME      ,TTIMEO     ,UNSYM      ,RSTINP     ,
     4   RSTOUT     ,0          ,IDUMMY     )
      ELSE
C Not re-start mode: Read main data from input data file
        CALL INDATA(MXFRON ,UNSYM)
C Regular analysis:        
        IF(NMULTI.EQ.1)THEN
C ...read and evaluate the applied external loads
          CALL INLOAD
C Multi-scale analysis:
        ELSEIF(NMULTI.EQ.2)THEN
C ...read prescribed array of macroscopic engineering strain (small 
C    strains) or deformation gradient (large strains)
          CALL RVSTRA
C ...split the nodes and degrees of freedom of the micro-cell mesh into
C different groups (interior, boundary, ...)
          CALL SPLMES
        ENDIF
      ENDIF

C For any mode: Read load incrementation data from input data file
      CALL ININCR
     1(   DFACT      ,DLENP      ,FSTOP      ,ITDES      ,MINCS      ,
     2    MITER      ,NALGO      ,NINCS      ,TOLER      ,
     3    DFACTV     ,DTIMEV     ,MITERV     ,NOUTP      ,NOUTPV     ,
     4    TOLERV     )
C
C Set prescribed component (uniform strain component) of all nodal
C displacements of RVE mesh and set RVE related submain variables
C ----------------------------------------------------------------      
      IF(NMULTI.EQ.2)THEN
        CALL SPLDSP(    DUTYLR      ,NINCS     )
      ENDIF
C
C Initialise some variables and arrays if not in re-start mode
C ------------------------------------------------------------
      INCRST=0
      IF(.NOT.RSTRT)THEN
        CALL INITIA
     1(   DLAMD      ,IFNEG      ,KUNLD      ,TFACT      ,TTIME      )
      ENDIF
C
C
C
C Start incremental finite element analysis...
C ********************************************
C
       WRITE(*,1030)
C
C=======================================================================
C                                                                      |
C                   Start loop over load increments                    |
C                                                                      |
C REFERENCE: Chapter 4 (Boxes 4.1-4) of the companion textbook.        |
C            Section 5.4.                                              |
C            The load incrementation loops carried out here are those  |
C            of the Flowcharts of Figures 5.2-3.                       |
C                                                                      |
C=======================================================================
C
      IF(.NOT.RSTRT)IINCS=0
C
      DO 50 ICOUNT=1,NINCS
C
        IPSUB=1
        IF(NALGO.GT.0)THEN
          DFSUB(1)=DFACTV(ICOUNT)
          DTSUB(1)=DTIMEV(ICOUNT)
          TOLER=TOLERV(ICOUNT)
          MITER=MITERV(ICOUNT)
          NOUTP(1)=NOUTPV(1,ICOUNT)
          NOUTP(2)=NOUTPV(2,ICOUNT)
          NOUTP(3)=NOUTPV(3,ICOUNT)
          NOUTP(4)=NOUTPV(4,ICOUNT)
          NOUTP(5)=NOUTPV(5,ICOUNT)
        ENDIF
C
C Reset converged problem variables
C ---------------------------------
        CALL SWITCH( 1 )
C
   10   CONTINUE
C
C Update increment counter
C ------------------------
        IINCS=IINCS+1
C
C For fixed increments option only: Increment external load according
C to user-prescribed incremental proportional load factor
C -------------------------------------------------------------------
        IF(NALGO.GT.0)THEN
          DFACT=DFSUB(IPSUB)
          DTIME=DTSUB(IPSUB)
          CALL INCREM
     1(   IINCS      ,TFACT      ,TOLER      ,TTIME     ,MITER      ,
     2    NOUTP      ,DFACT      ,DFOLD      ,DTIME     ,KUNLD      )
        ENDIF
C
        IF(NMULTI.EQ.2)THEN
          DO I=1,NTOTV
            DITER(I)=DUTYLR(I)*DFACT
          ENDDO
C Update incremental and total displacements. Also update nodal
C coordinates for large deformation analyses
C -------------------------------------------------------------
          CALL UPCONF
C
C Re-set relevant problem variables to last converged solution
C ------------------------------------------------------------
          CALL SWITCH( 2 )
C
C Update problem variables (stress and other state variables) and
C evaluate internal force vectors of all elements
C ---------------------------------------------------------------
          CALL INTFOR( DTIME, INCCUT, TTIME )
C          
C Update ELOAD (element internal forces)
          CALL CONVRV(CONVRG,DIVERG,IITER,TOLER,TFACT,RESVEC,FBD)
C 
        ENDIF
C______________________________________________________________________
C                                                                      |
C             Start loop over equilibrium iterations                   |
C______________________________________________________________________|
C
        DO 20 IITER=1,MITER
C
C Select solution algorithm variable KRESL
          CALL ALGOR(IINCS ,IITER ,KRESL ,KUNLD )
C
          IF(NALGO.LT.0)THEN
C Set up prescribed displacements for tangential solution for the
C arc-length method
            CALL TANGEN
          ENDIF
C          
C REGULAR ANALYSIS:
C Assemble stiffness matrix and solve for iterative displacements
C (tangential solution for the arc-length method) the linearised system
C of discretised equilibrium equations using the chosen algorithm
C ---------------------------------------------------------------------
          IF(NMULTI.EQ.1)THEN
            CALL SOLINT
     1(   DTIME      ,IITER      ,KRESL      ,IFNEG      ,KUNLD      ,
     2    MXFRON     ,UNSYM      ,INCCUT     ,NSOLVE     )
C
C RVE ANALYSIS:
C Assemble the reduced stiffness matrices and solve for iterative 
C displacements the linearised system of discretised equilibrium 
C equations based on the chosen kinematical constraint option
C ---------------------------------------------------------------
          ELSEIF(NMULTI.EQ.2)THEN
            CALL SOLRVE(DTIME ,KUNLD ,UNSYM ,INCCUT, RESVEC )
          ENDIF
C
          IF(INCCUT)THEN
C System solution failed due to zero pivot: break equilibrium iteration
C loop and activate increment cutting
            GOTO 30
          ENDIF
C
C For Arc-Length method only: Compute iterative displacement according
C to the arc-length constraint and update the incremental and total
C load factors
C --------------------------------------------------------------------
          IF(NALGO.LT.0)THEN
            CALL ARCLEN
     1(   DFACT      ,DLAMD      ,DLENG      ,DLENM      ,DLENP      ,
     2    IFNEG      ,IINCS      ,IITER      ,INCCUT     ,TFACT      )
C
            IF(INCCUT)THEN
C No real roots for arc-length constraint equation: break equilibrium
C iteration loop and activate increment cutting
              GOTO 30
            ENDIF
          ENDIF
C
C Update incremental and total displacements. Also update nodal
C coordinates for large deformation analyses
C -------------------------------------------------------------
          CALL UPCONF
C
C Re-set converged load factors and print out increment information
C -----------------------------------------------------------------
          IF(IITER.EQ.1)THEN
            IF(IINCS.EQ.1)THEN
C Re-set previous converged load factors/arc-length
              IF(NALGO.LT.0)DLENGO=DLENG
              DFACTO=DFACT
              TFACTO=R0
            ENDIF
            IF(NALGO.GT.0)THEN
C Fixed increments option: print current total load factor
               WRITE(*,1050) IINCS,TTIME,TFACT
c               WRITE(16,1050)IINCS,TTIME,TFACT
            ELSE
C Arc-length: print current arc-length
               WRITE(*,1055) IINCS,DLENG
c               WRITE(16,1055)IINCS,DLENG
            ENDIF
          ENDIF
C
C Re-set relevant problem variables to last converged solution
C ------------------------------------------------------------
          CALL SWITCH( 2 )
C
C Update problem variables (stress and other state variables) and
C evaluate internal force vectors of all elements
C ---------------------------------------------------------------
          CALL INTFOR( DTIME, INCCUT, TTIME )
C
          IF(INCCUT)THEN
C Internal force calculation failed: break equilibrium iteration loop
C and activate load increment cutting
            GOTO 30
          ENDIF
C
C Assemble internal and external global force vectors, reactions,
C compute residual and check for convergence
C ---------------------------------------------------------------
          IF(NMULTI.EQ.1)THEN
            CALL CONVER(CONVRG,DIVERG,IITER,TOLER,TFACT)
          ELSEIF(NMULTI.EQ.2)THEN
            CALL CONVRV(CONVRG,DIVERG,IITER,TOLER,TFACT,RESVEC
     1                        ,FBD)
      
          ENDIF
C
          ITACT=IITER
C
          IF(CONVRG)THEN
C Iterations have converged: break equilibrium iteration loop and go to
C next load increment
             WRITE(*,1060)
c             WRITE(16,1060)
            IF(NALGO.GT.0)THEN
               WRITE(*,1063) DTIME,DFACT
c               WRITE(16,1063)DTIME,DFACT
            ELSE
               WRITE(*,1065) TFACT
               WRITE(*,1067) DFACT
c               WRITE(16,1065)TFACT
c               WRITE(16,1067)DFACT
            ENDIF
             WRITE(*,1060)
c             WRITE(16,1060)
            GOTO 40
          ELSEIF(DIVERG)THEN
C Iterations are diverging: break equilibrium iteration loop and
C activate load increment cutting
c             WRITE(16,1100)
             WRITE(*,1100)
            GOTO 30
          ENDIF
C
   20   CONTINUE
C______________________________________________________________________ 
C                                                                      |
C                End loop over equilibrium iterations                  |
C______________________________________________________________________|
C
C Newton-Raphson procedure did not converge within the prescribed
C maximum number of iterations !!
C Print corresponding message and proceed to increment cutting
C
        WRITE(16,1040)
        WRITE(*,1040)
C
C
C
   30   CONTINUE
CRESVEC
C
C Activate increment cutting
C
C REFERENCE: Section 5.4.3
C --------------------------
C
        IF(NALGO.GT.0)THEN
C For fixed increments option: split current load increment into two
C equally sized sub-increments
          WRITE(16,1110)
          WRITE(*,1110)
          IF(IPSUB.EQ.MSUBIN)THEN
C abort program if maximum permissible number of consecutive increment
C cuts has been exceeded (i.e. sub-increment stack array DFSUB is full)
            CALL ERRPRT('EE0002')
          ENDIF
          DFSUB(IPSUB)  =DFSUB(IPSUB)*RP5
          DFSUB(IPSUB+1)=DFSUB(IPSUB)
          DTSUB(IPSUB)  =DTSUB(IPSUB)*RP5
          DTSUB(IPSUB+1)=DTSUB(IPSUB)
          IPSUB=IPSUB+1
        ELSE
C For arc-length method: reduce the arc-length
          WRITE(16,1120)
          WRITE(*,1120)
          IF(IINCS.EQ.1)THEN
            DFACT=DFACTO*RP7
            DFACTO=DFACT
          ELSE
            DLENG=DLENGO*RP7
            DLENGO=DLENG
          ENDIF
        ENDIF
C Switch relevant variables to last converged values (in load increment
C cutting mode) before re-trying with reduced load increment/arc-length
        TFACT=TFACTO
        TTIME=TTIMEO
        CALL SWITCH( 3 )
        IINCS=IINCS-1
        GOTO 10
C
C
C
   40   CONTINUE
C
C Newton-Raphson iterations converged for the current load increment
C ------------------------------------------------------------------
C Reset some converged parameters
        IF(NALGO.LT.0)DLENGO=DLENG
        TFACTO=TFACT
        TTIMEO=TTIME
        IF(NALGO.GT.0)THEN
C Fixed increments option: update pointer to sub-increments stack array
          IPSUB=IPSUB-1
        ELSE
C Arc-length method: update arc-length according to the desired number
C of iterations for convergence and the actual number of iterations
C needed for convergence in the previous load step
          CALL LENGTH(DLENG ,DLENM ,ITACT ,ITDES )
        ENDIF
C
C Output results if required
C REFERENCE: Section 5.4.7
C
C
C Compute and display homogenised strains and stresses
        IF(NMULTI.EQ.2)THEN
          CALL RVOUTP( IINCS, FBD, TFACT )
          CALL OUTPUT(TFACT,TTIME,IINCS,IITER,NOUTP,FBD)
        ENDIF
C Produce vtk files
        CALL VTKOUT ( IINCS, TFACT, DATFIL, DUTYLR )
c        CALL OUTPUT(TFACT,TTIME,IINCS,IITER,NOUTP)
        IF((NALGO.GT.0.AND.IPSUB.EQ.0).OR.(NALGO.LT.0))THEN
          CALL RSTART
     1(  DFOLD      ,DLENG      ,DLENGO     ,DLENM      ,DLAMD      ,
     2   IFNEG      ,IINCS      ,MXFRON     ,NOUTP      ,TFACT      ,
     3   TFACTO     ,TTIME      ,TTIMEO     ,UNSYM      ,RSTINP     ,
     4   RSTOUT     ,1          ,INCRST     )
        ELSEIF(NALGO.GT.0.AND.IPSUB.NE.0)THEN
          CALL SWITCH( 1 )
          GOTO 10
        ENDIF
        IF(NALGO.LT.0.AND.FSTOP.NE.R0.AND.TFACT.GT.FSTOP)THEN
C Arc-length only: Break loop over increments and stop if maximum
C prescribed load factor has been exceeded
          GOTO 60
        ENDIF
C
C
   50 CONTINUE
C
C=======================================================================
C                                                                      |
C                   End loop over load increments                      |
C                                                                      |
C=======================================================================
C
   60 CONTINUE
C
C
C
C Exit HYPLAS
C ***********
C
C Close files before exit
      CALL FCLOSE
C Echo file names back to standard output and stop
       WRITE(*,1070)
      I=INDEX(RESFIL,' ')-1
       WRITE(*,1080)DATFIL(1:I)
       WRITE(*,1090)RESFIL(1:I)
      IF(IRV_RVEOPT.EQ.1)THEN
          WRITE(*,1092) 'Linear'
      ELSE IF(IRV_RVEOPT.EQ.2)THEN 
          WRITE(*,1092) 'Periodic'
      ELSE IF(IRV_RVEOPT.EQ.3)THEN
      	   WRITE(*,1092) 'Uniform traction'
      ENDIF	
      IF(INCRST.NE.0)THEN
         WRITE(*,1095)RSTOUT(1:I),INCRST
      ENDIF
      STOP ' '
      END
CDOC END_PROGRAM HYPLAS
