CDOC BEGIN_SUBROUTINE ININCR
CDOC Reads input data for load increments
CDOC
CDOC This routine reads from the data file and echos to the results file
CDOC the data required for load incrementation.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DFACT  <  Initial incremental load factor (for
CDOC C                          Arc-Length control).
CDOC DOUBLE_PRECISION DLENP  <  Maximum arc-length parameter (for
CDOC C                          Arc-Length control).
CDOC DOUBLE_PRECISION FSTOP  <  Maximum total load factor, above which
CDOC C                          the analysis will stop (for
CDOC C                          Arc-Length control).
CDOC INTEGER          ITDES  <  Desired number of equilibrium
CDOC C                          iterations (for Arc-Length control).
CDOC INTEGER          MINCS  >  Dimensioning parameter: maximum
CDOC C                          permissible number of prescribed load
CDOC C                          increments under fixed load increments
CDOC C                          option. Defines size of
CDOC C                          DFACTV, DTIMEV, MITERV and NOUTPV.
CDOC INTEGER          MITER  <  Maximum allowed number of equilibrium
CDOC C                          iterations in any increment (for
CDOC C                          Arc-Length control).
CDOC INTEGER          NALGO  >  Solution algorithm flag (positive for
CDOC C                          fixed increments option, negative for
CDOC C                          Arc-Length control).
CDOC INTEGER          NINCS  <  For fixed increments option: total
CDOC C                          number of specified load increments.
CDOC C                          For Arc-Length control: maximum
CDOC C                          specified number of load increments.
CDOC DOUBLE_PRECISION TOLER  <  Equilibrium convergence for all load
CDOC C                          increments (for Arc-Length control).
CDOC DOUBLE_PRECISION DFACTV <  Array with the incremental load factor
CDOC C                          of each specified load increment (for
CDOC C                          fixed load increments option).
CDOC DOUBLE_PRECISION DTIMEV <  Array with the time increment for each
CDOC C                          specified load increment (for fixed
CDOC C                          load increments option).
CDOC INTEGER          MITERV <  Array with the maximum allowed number of
CDOC C                          global (equilibrium) iterations for each
CDOC C                          specified load increment (for
CDOC C                          fixed load increments option).
CDOC INTEGER          NOUTP  <  Output code array with output frequency
CDOC C                          flags (for Arc-Length control).
CDOC INTEGER          NOUTPV <  Array with the output code for each
CDOC C                          specified load increment (for fixed load
CDOC C                          increments option).
CDOC DOUBLE_PRECISION TOLERV <  Array with the equilibrium convergence
CDOC C                          tolerance for each specified load
CDOC C                          increment (for fixed load increments
CDOC C                          option).
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, July 1996:  Initial coding
CHST
CHST E.de Souza Neto, June 1998:  Arc-length data reading included
CHST
      SUBROUTINE ININCR
     1(   DFACT      ,DLENP      ,FSTOP      ,ITDES      ,MINCS      ,
     2    MITER      ,NALGO      ,NINCS      ,TOLER      ,
     3    DFACTV     ,DTIMEV     ,MITERV     ,NOUTP      ,NOUTPV     ,
     4    TOLERV     )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION
     1    DFACTV(MINCS)      ,DTIMEV(MINCS)      ,MITERV(MINCS)      ,
     2    NOUTP(5)           ,NOUTPV(5,MINCS)    ,TOLERV(MINCS)
C
      LOGICAL FOUND
      CHARACTER*80 INLINE
      DIMENSION
     1    IWBEG(40), IWEND(40)
      PARAMETER
     1(   R0=0.0D0   )
C***********************************************************************
C READS INPUT DATA FOR LOAD INCREMENTATION
C
C REFERENCE: Figure 5.1
C***********************************************************************
 1000 FORMAT(///' Increment control with fixed load increments selected'
     1         /' ====================================================='
     2        //'       Number of proportional load increments =',I5/)
C
 1010 FORMAT(///' Increment control by the Arc-Length method selected'/
     1          ' ==================================================='//
     2          ' Arc length data'/' ---------------'//
     3' Maximum allowed number of increments ........= ',I5)
 1020 FORMAT(
     1' Initial load increment factor .............. = ',G15.6/
     2' Convergence tolerence ...................... = ',G15.6/
     3' Max. No. of iterations ..................... = ',I5)
 1030 FORMAT(/' Output control parameter for results'/,
     1        '       ( Output frequencies )'/
     2    ' Displacements......................... = ',I3/
     3    ' Reactions............................. = ',I3/
     4    ' State variables at gauss points....... = ',I3/
     5    ' State variables at nodes.............. = ',I3/
     6    ' Output to re-start file............... = ',I3)
 1040 FORMAT(/
     1' Desired number of iterations per increment .. =',I5/
     2' Maximum load factor ......................... =',G15.6/
     3' Maximum arc length parameter ................ =',G15.6)
C
      CALL FNDKEY
     1(   FOUND      ,IWBEG      ,IWEND      ,'INCREMENTS',
     2    INLINE     ,15         ,NWRD       )
      IF(.NOT.FOUND)CALL ERRPRT('ED0079')
C
      IF(NALGO.GT.0)THEN
C Fixed increments option
C -----------------------
        IF(NWRD.EQ.1)CALL ERRPRT('ED0034')
        NINCS=INTNUM(INLINE(IWBEG(2):IWEND(2)))
c        WRITE(16,1000)NINCS
        IF(NINCS.LE.0)    CALL ERRPRT('ED0013')
        IF(NINCS.GT.MINCS)CALL ERRPRT('ED0035')
        DO 10 IINCS=1,NINCS
          READ(15,*,ERR=997,END=997)DTIMEV(IINCS),DFACTV(IINCS),
     1                              TOLERV(IINCS),MITERV(IINCS),
     2                              (NOUTPV(I,IINCS),I=1,5)
          IF(TOLERV(IINCS).LE.R0)CALL ERRPRT('ED0142')
          IF(DTIMEV(IINCS).LT.R0)CALL ERRPRT('ED0227')
          IF(MITERV(IINCS).LT.1)CALL ERRPRT('ED0146')
   10   CONTINUE
      ELSE
C Arc-length control
C ------------------
        IF(NWRD.EQ.1) CALL ERRPRT('ED0097')
        NINCS=INTNUM(INLINE(IWBEG(2):IWEND(2)))
c        WRITE(16,1010)NINCS
        IF(NINCS.LE.0)CALL ERRPRT('ED0098')
C
        READ(15,*,ERR=998,END=998)DFACT,TOLER,MITER,(NOUTP(I),I=1,5),
     1                            ITDES,FSTOP,DLENP
        IF(TOLER.LT.R0)CALL ERRPRT('ED0142')
        IF(ITDES.LT.0)CALL ERRPRT('ED0143')
c        WRITE(16,1020)DFACT,TOLER,MITER
c        WRITE(16,1030)(NOUTP(I),I=1,5)
c        WRITE(16,1040)ITDES,FSTOP,DLENP
      ENDIF
C Send error messages in case of I/O error while reading increment data
      GOTO 999
  997 CALL ERRPRT('ED0178')
      GOTO 999
  998 CALL ERRPRT('ED0179')
  999 CONTINUE
      RETURN
      END
CDOC END_SUBROUTINE ININCR
