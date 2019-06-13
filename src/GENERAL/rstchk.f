CDOC BEGIN_SUBROUTINE RSTCHK
CDOC Checks wether the main data is to be read from a re-start file
CDOC
CDOC This routine reads the data file and checks wether the main data
CDOC is to be read from it or from a re-start file.
CDOC
CDOC BEGIN_PARAMETERS
CDOC CHARACTER        RSTINP >  Character string containing the input
CDOC C                          re-start file name.
CDOC LOGICAL          RSTRT  <  Restart flag. Return value is set to
CDOC C                          .FALSE. if the main data is to
CDOC C                          be read from the data file. Set to
CDOC C                          .TRUE. if the main data is to
CDOC C                          be read from a re-start file.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, September 1996: Initial coding
CHST E.de Souza Neto, September 2009: Added check for missing re-start
CHST                                  file name in data file. HYPLAS now
CHST                                  stops if name is missing.
CHST                                  The absence of this check was
CHST                                  causing HYPLAS to crash
CHST                                  (segmentation fault) when the
CHST                                  re-start file name was missing.
CHST
      SUBROUTINE RSTCHK(  RSTINP     ,RSTRT      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL RSTRT
      CHARACTER*256 RSTINP
C
      LOGICAL AVAIL,FOUND
      CHARACTER*80 INLINE
      DIMENSION IWBEG(40),IWEND(40)
C***********************************************************************
C CHECKS WETHER MAIN DATA IS TO BE READ FROM INPUT RE-START FILE
C AND SET INPUT RE-START FILE NAME IF REQUIRED
C***********************************************************************
 1000 FORMAT(////,
     1'               Main input data read from re-start file'/
     2'               ======================================='///
     3'               Input re-start file name ----> ',A)
C
C Checks whether the input data file contains the keyword RESTART
C
        CALL FNDKEY
     1(   FOUND     ,IWBEG    ,IWEND    ,'RESTART',
     2    INLINE    ,15       ,NWRD     )
        IF(FOUND)THEN
c checks if the re-start file name root is missing from the data file
          IF(NWRD.LT.2) CALL ERRPRT('ED0199')
C sets re-start flag and name of input re-start file
          RSTRT=.TRUE.
          RSTINP=INLINE(IWBEG(2):IWEND(2))//'.rst'
          WRITE(16,1000)INLINE(IWBEG(2):IWEND(2))//'.rst'
          WRITE(18,1000)INLINE(IWBEG(2):IWEND(2))//'.rst'
          WRITE( *,1000)INLINE(IWBEG(2):IWEND(2))//'.rst'
C checks existence of the input re-start file
          INQUIRE(FILE=RSTINP,EXIST=AVAIL)
          IF(.NOT.AVAIL)CALL ERRPRT('ED0096')
        ELSE
          RSTRT=.FALSE.
        ENDIF
C
      RETURN
      END
CDOC END_SUBROUTINE RSTCHK
