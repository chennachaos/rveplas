CDOC BEGIN_SUBROUTINE ERRPRT
CDOC Prints error message and can abort the program if requested
CDOC
CDOC This routine prints error/warning messages to the standard output
CDOC and results file and may abort the program depending on the entry
CDOC value of its argument. The character string passed as its argument
CDOC is an error code which will be searched for in the file ERROR.RUN,
CDOC assumed to be kept in the directory defined by the HYPLASHOME
CDOC environment variable. If the correponding character string (error
CDOC code) is found, the associated error/warning message is printed
CDOC in the standard output and results file.
CDOC The are 5 types of errors/warnings: Input data error, input data
CDOC warning, internal error, execution error and execution warning.
CDOC These are characterised, respectively, by error codes of the
CDOC types: 'ED????', 'WD????', 'EI????', 'EE????' and 'WE????'.
CDOC This routine aborts the program in case of input data error
CDOC (ED????), internal error ('EI????') or execution error (EE????).
CDOC The execution of the program is not interrupted in case of
CDOC warnings ('WE????' or 'WD????').
CDOC See file ERROR.RUN for more details.
CDOC 
CDOC BEGIN_PARAMETERS
CDOC CHARACTER        ERRCOD >  Character string containing the error
CDOC C                          code.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, July 1996: Initial coding
CHST
CHST E.de Souza Neto, October 2008: Handling of HYPLASHOME environment
CHST                                variable modified to cope with blank
CHST                                spaces in directory path name
CHST                                (this is common in WINDOWS).
CHST
CHST E.de Souza Neto, April 2011: Error message added in case of I/O
CHST                              error while reading file ERROR.RUN
CHST
      SUBROUTINE ERRPRT(ERRCOD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL       AVAIL ,FOUND
      CHARACTER*6   ERRCOD
      CHARACTER*22  HEADER
      CHARACTER*80  INLINE
      CHARACTER*256 HYPLASHOME
      DIMENSION IWBEG(128), IWEND(128)
C***********************************************************************
C PRINTS (TO THE STANDARD OUTPUT AND RESULTS FILE) MOST ERROR/WARNING
C MESSAGES OF HYPLAS.
C THE ARGUMENT "ERRCOD" CONTAINS THE ERROR/WARNING MESSAGE CODE. THE
C TEXT OF THE MESSAGE ASSOCIATED WITH THE GIVEN CODE MUST BE IN FILE
C "ERROR.RUN". THIS FILE IS KEPT IN THE FILE SYSTEM DIRECTORY WHOSE
C NAME IS STORED IN THE OPERATING SYSTEM ENVIRONMENT VARIABLE
C "HYPLASHOME" (SEE COMMENTS BELOW).
C***********************************************************************
 1000 FORMAT(///' ',74('*')/' *',72X,'*'/' *',25X,A22,25X,'*'/' *',72X,
     1      '*'/' ',74('*')/' *',72X,'*'/' * Code:         ',A6,51X,'*')
 1010 FORMAT(' *',A72,'*')
 1020 FORMAT(' *',72X,'*'/' ',74('*'))
 1030 FORMAT(//'   ASSOCIATED MESSAGES WILL NOT BE PRINTED !',//,
     1  '   The above error code was not found in file ERROR.RUN.',//)
 1040 FORMAT(//'   ASSOCIATED MESSAGES WILL NOT BE PRINTED !',//,
     1   '   File ERROR.RUN does not exist in HYPLASHOME directory.',//)
 1050 FORMAT(//'   THERE WAS AN I/O ERROR WHILE READING ERROR.RUN',//,
     1   '   The full error message will not be printed.',//)
C
      IF(ERRCOD(1:2).EQ.'ED')THEN
        HEADER='  INPUT  DATA  ERROR  '
      ELSEIF(ERRCOD(1:2).EQ.'WD')THEN
        HEADER=' INPUT  DATA  WARNING '
      ELSEIF(ERRCOD(1:2).EQ.'EI')THEN
        HEADER='   INTERNAL   ERROR   '
      ELSEIF(ERRCOD(1:2).EQ.'EE')THEN
        HEADER='    EXECUTION ERROR   '
      ELSEIF(ERRCOD(1:2).EQ.'WE')THEN
        HEADER='   EXECUTION WARNING  '
      ELSE
        HEADER='  UNKNOWN ERROR TYPE  '
      ENDIF
      WRITE(*,1000)HEADER,ERRCOD
      WRITE(16,1000)HEADER,ERRCOD
C
C
C WARNING: GETENV is a non-standard FORTRAN 77 instruction, used here to
C          obtain the value of the operating system environment variable
C          HYPLASHOME - A character string containig the name of the
C          directory where the file ERROR.RUN (containing the error/
C          warning messages of HYPLAS) is kept in the file system. You
C          may need to change this if your FORTRAN compiler does not
C          support the instruction GETENV.
      CALL GETENV('HYPLASHOME',HYPLASHOME)
C
C Opens file ERROR.RUN
      NWRD=NWORD(HYPLASHOME,IWBEG,IWEND)
      !INQUIRE(FILE=HYPLASHOME(1:IWEND(NWRD))//'/ERROR.RUN',EXIST=AVAIL)
      AVAIL=.FALSE.
      IF(.NOT.AVAIL)THEN
        WRITE(*,1040)
        WRITE(16,1040)
        GOTO 999
      ENDIF
      OPEN(23,FILE=HYPLASHOME(1:IWEND(NWRD))//'/ERROR.RUN',STATUS='OLD')
C
C Finds the character string with the given error code in file ERROR.RUN
      CALL FNDKEY
     1(   FOUND      ,IWBEG      ,IWEND      ,ERRCOD     ,INLINE     ,
     2    23         ,NWRD       )
      IF(.NOT.FOUND)THEN
        WRITE(*,1030)
        WRITE(16,1030)
        GOTO 998
      ENDIF
C
C Reads and echoes the associated error/warning message
      NLINES=INTNUM(INLINE(IWBEG(2):IWEND(2)))
      DO 10 I=1,NLINES
        READ(23,'(A72)',ERR=900,END=900)INLINE
        WRITE(*,1010)INLINE
        WRITE(16,1010)INLINE
   10 CONTINUE
      WRITE(*,1020)
      WRITE(16,1020)
C I/O error while reading file ERROR.RUN
      GOTO 998
  900 WRITE(*,1050) 
      WRITE(16,1050) 
C
  998 CONTINUE
      CLOSE(23,STATUS='KEEP')
  999 CONTINUE
C
C All codes WE???? and WD???? are WARNING codes (execution warnings and
C input data warnings, respectively) and do not stop the execution
C of HYPLAS. Any other type of message code is seen as an ERROR code and
C will cause HYPLAS to stop its execution. Refer to comments in the
C file ERROR.RUN.
      IF(ERRCOD(1:2).NE.'WE'.AND.ERRCOD(1:2).NE.'WD')CALL PEXIT
      RETURN
      END
CDOC END_SUBROUTINE ERRPRT
