CDOC BEGIN_SUBROUTINE FNDKEY
CDOC Finds and reads a line containing a specified keyword from a file.
CDOC
CDOC This routine searches for a given keyword positioned as the first
CDOC word of a line in a file.
CDOC If the given keyword is found then the corresponding line is read
CDOC and returned together with the number of words in the line and two
CDOC integer arrays containing the position of the beginning and end of
CDOC each word.
CDOC If there is more than one line in the given file containing the
CDOC the specified keyword as its first word, then this routine will
CDOC return the first line found in the search which, in general, does
CDOC not coincide with the first occurrence of the keyword in the file.
CDOC In this case, which line is returned depends on which line of the
CDOC file the search started. Hence, it is recommended that files to be
CDOC searched contain only one occurrence of the keyword as first word
CDOC of a line (occurrences of the keyword in positions other than as
CDOC the first word of a line are not detected by this routine and can
CDOC safely be present in the file).
CDOC
CDOC BEGIN_PARAMETERS
CDOC LOGICAL          FOUND  <  Logical flag. Its return value is set to
CDOC C                          .TRUE. if the specified keyword
CDOC C                          is found at the beginning of a line in
CDOC C                          the specified file and set to
CDOC C                          .FALSE. otherwise.
CDOC INTEGER          IWBEG  <  Integer array containg the position of
CDOC C                          the beginning of each word in the line
CDOC C                          containing the specified keyword.
CDOC C                          IWBEG(K) is the position of the
CDOC C                          beginning of the Kth word in the line.
CDOC INTEGER          IWEND  <  Integer array containg the position of
CDOC C                          the end of each word in the line
CDOC C                          containing the specified keyword.
CDOC C                          IWEND(K) is the position of the end of
CDOC C                          the Kth word in the line.
CDOC CHARACTER        KEYWRD >  Keyword.
CDOC CHARACTER        INLINE <  Contents of the line whose first word is
CDOC C                          the specified  keyword.
CDOC INTEGER          NFILE  >  Unit identifier of the file to be
CDOC C                          searched.
CDOC INTEGER          NWRD   <  Total number of words in the line
CDOC C                          containing the specified keyword.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, July 1996: Initial coding
CHST
      SUBROUTINE FNDKEY
     1(   FOUND      ,IWBEG      ,IWEND      ,KEYWRD     ,INLINE     ,
     2    NFILE      ,NWRD       )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL FOUND
      CHARACTER*80 INLINE
      CHARACTER*(*) KEYWRD
      DIMENSION
     1    IWBEG(40), IWEND(40)
C***********************************************************************
C FINDS AND READS A LINE CONTAINING A SPECIFIED KEYWORD FROM A FILE.
C THIS ROUTINE SEARCHES FOR A GIVEN KEYWORD POSITIONED AS THE FIRST
C WORD OF A LINE IN A FILE.
C IF THE GIVEN KEYWORD IS FOUND THEN THE CORRESPONDING LINE IS READ AND
C RETURNED TOGETHER WITH THE NUMBER OF WORDS IN THE LINE AND TWO INTEGER
C ARRAYS CONTAINING THE POSITION OF THE BEGINNING AND END OF EACH WORD.
C***********************************************************************
 1000 FORMAT(A80)
C
      FOUND=.TRUE.
      IEND=0
   10 READ(NFILE,1000,END=20)INLINE
      NWRD=NWORD(INLINE,IWBEG,IWEND)
      IF(NWRD.NE.0)THEN
        IF(INLINE(IWBEG(1):IWEND(1)).EQ.KEYWRD)THEN
          GOTO 999
        ENDIF
      ENDIF
      GOTO 10
   20 IF(IEND.EQ.0)THEN
        IEND=1
        REWIND NFILE
        GOTO 10
      ELSE
        REWIND NFILE
        FOUND=.FALSE.
      ENDIF
  999 RETURN
      END
CDOC END_SUBROUTINE FNDKEY
