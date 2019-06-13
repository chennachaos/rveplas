CDOC BEGIN_SUBROUTINE FOPEN
CDOC Opens/sets the names of the input data, results and re-start files
CDOC
CDOC This routine reads the input data file name from the standard
CDOC input, sets the names and opens the corresponding data and results
CDOC files, and then sets the re-start file name.
CDOC
CDOC BEGIN_PARAMETERS
CDOC CHARACTER         DATFIL <  Character string containing the input
CDOC C                           data file name.
CDOC CHARACTER         RESFIL <  Character string containing the results
CDOC C                           output file name.
CDOC CHARACTER         RSTOUT <  Character string containing the
CDOC C                           re-start output file name.
CDOC INTEGER           NMULTI <  Multi-scale analysis type flag. Set
CDOC                             according to input file extension
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, June 1996: initial coding
CHST
CHST D. de Bortoli,  March 2014: reading of RVE input data files added 
CHST
      SUBROUTINE FOPEN
     1(   DATFIL     ,RESFIL     ,RSTOUT     ,
     2    NMULTI     ,PLOTD )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER DATFIL*256,RESFIL*256,RSTOUT*256,PLOTD*256,MYDATA*256
      INTEGER NMULTI
      CHARACTER BELL*1
      LOGICAL AVAIL
C***********************************************************************
C READ DATA FILE NAME FROM STANDARD INPUT, SET NAMES FOR AND OPEN DATA
C AND RESULTS FILES AND SET RE-START FILE NAME
C
C REFERENCE: Figure 5.1
C***********************************************************************
C Read data file name from standard input
C ---------------------------------------
c       WRITE(*,'(////15X,A,/15X,A)')
c     1'Data file name must have extension .dat/.DAT or .rve/.RVE',
c     2'and must not contain blank spaces.'
c      WRITE(*,'(/15X,A)')'(Type EXIT or QUIT to terminate)'
!   10 WRITE(*,'(//15X,A$)')'Input data file name --------> '
!      READ(*,'(A)',ERR=10)DATFIL
C Sort out data, result and re-start file names and open data and
C result files
C ---------------------------------------------------------------
      BELL=CHAR(7)
      DATFIL = "8hexaele.rve"
C      DATFIL2='MASH.rve'
      MYDATA="MAHSHID.rve"
C Find end of data file name
      I=INDEX(DATFIL,' ')-1
      IF(I.EQ.0)THEN
        WRITE(*,'(/15X,A)')
     1          'Data file name must NOT begin with blank space !'
!        GOTO 10
      ELSEIF(I.EQ.4.AND.
     1         (DATFIL(1:4).EQ.'EXIT'.OR.DATFIL(1:4).EQ.'exit'
     2      .OR.DATFIL(1:4).EQ.'QUIT'.OR.DATFIL(1:4).EQ.'quit'))THEN
        WRITE(*,'(///15X,A,///)')'Program HYPLAS terminated by user.'
        STOP ' '
      ENDIF
C Check data file name extension...
C     ... regular analysis: it must be either .dat or .DAT
      IF(DATFIL(I-3:I).EQ.'.dat'.OR.DATFIL(I-3:I).EQ.'.DAT')THEN
        NMULTI=1
C     ... multi-scale analysis: it must be either .rve or .RVE
      ELSEIF(DATFIL(I-3:I).EQ.'.rve'.OR.DATFIL(I-3:I).EQ.'.RVE')THEN
        NMULTI=2
      ELSE
         WRITE(*,'(/15X,A,/15X,A)')
     1          'Data file name does not have a valid extension!',
     2          '(.dat, .DAT, .rve or .RVE) Please try again'
         WRITE(*,'(1X,A$)')BELL
!        GOTO 10
      ENDIF
C Check existence of data file
      INQUIRE(FILE=DATFIL,EXIST=AVAIL)
      IF(.NOT.AVAIL)THEN
        WRITE(*,'(/A,A,A,A)')
     1'               File "',DATFIL(1:I),'" not found !  ',
     2              ' Please try again'
        WRITE(*,'(1X,A$)')BELL
!        GOTO 10
      ENDIF
C Give name to results file...
C     ... regular analysis: extension .res
      IF(NMULTI.EQ.1)THEN
        RESFIL(1:I-3)=DATFIL(1:I-3)
        RESFIL(I-3:I)='.res'
        RESFIL(I+1:256)=DATFIL(I+1:256)

        PLOTD(1:I-3)=MYDATA(1:I-3)
        PLOTD(I-3:I)='.rvr'
        PLOTD(I+1:256)=MYDATA(I+1:256)
C     ... multi-scale analysis: extension .rvr
      ELSEIF(NMULTI.EQ.2)THEN
        RESFIL(1:I-3)=DATFIL(1:I-3)
        RESFIL(I-3:I)='.rvr'
        RESFIL(I+1:256)=DATFIL(I+1:256)
        PLOTD(1:I-3)=MYDATA(1:I-3)
        PLOTD(I-3:I)='.rvr'
        PLOTD(I+1:256)=MYDATA(I+1:256)
      ENDIF
C Give name to output re-start file (extension .rst)
        RSTOUT(1:I-3)=DATFIL(1:I-3)
        RSTOUT(I-3:I)='.rst'
        RSTOUT(I+1:256)=DATFIL(I+1:256)
        RSTOUT(1:I-3)=MYDATA(1:I-3)
        RSTOUT(I-3:I)='.rst'
        RSTOUT(I+1:256)=MYDATA(I+1:256)
C Open data and results file
      OPEN(UNIT=15,FILE=DATFIL,STATUS='OLD')
      OPEN(UNIT=16,FILE=RESFIL,STATUS='UNKNOWN')
      OPEN(UNIT=18,FILE=PLOTD,STATUS='UNKNOWN')
C      OPEN(UNIT=13,FILE=FNAME,STATUS='OLD',BLANK='ZERO',ERR=100)
C
      RETURN
      END
CDOC END_SUBROUTINE FOPEN
