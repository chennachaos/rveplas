CDOC BEGIN_INTEGER_FUNCTION NWORD
CDOC Returns the number of words contained in a character string
CDOC
CDOC The return value of this function is the number of words contained
CDOC in the character string passed in its argument list. The function
CDOC also sets the pointers to the beginning and end of each word. The
CDOC pointers are returned via argument list.
CDOC
CDOC BEGIN_PARAMETERS
CDOC CHARACTER        CHRSTR >  Character string.
CDOC INTEGER          IWBEG  <  Array of pointers to the beginning of
CDOC C                          the words contained in CHRSTR.
CDOC C                          For the Nth word, beginning at
CDOC C                          CHRSTR(I:I), the function sets
CDOC C                          IWBEG(N)=I.
CDOC INTEGER          IWEND  <  Array of pointers to the end of
CDOC C                          the words contained in CHRSTR.
CDOC C                          For the Nth word, ending at
CDOC C                          CHRSTR(I:I), the function sets
CDOC C                          IWEND(N)=I.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, July 1996: Initial coding
CHST
      INTEGER FUNCTION NWORD(CHRSTR,IWBEG,IWEND)
      IMPLICIT NONE
      CHARACTER*(*) CHRSTR
      INTEGER IWBEG(*), IWEND(*)
      LOGICAL OUT
      INTEGER I, LEN, LENGTH
C***********************************************************************
C FIND NUMBER OF WORDS CONTAINED IN A CHARACTER STRING AND SET POINTERS
C TO BEGINNING AND END OF EACH WORD
C***********************************************************************
      LENGTH=LEN(CHRSTR)
      NWORD=0
      OUT=.TRUE.
      DO 10 I=1,LENGTH
        IF(OUT)THEN
          IF(CHRSTR(I:I).NE.' ')THEN
            OUT=.FALSE.
            NWORD=NWORD+1
            IWBEG(NWORD)=I
          ENDIF
        ELSE
          IF(CHRSTR(I:I).EQ.' ')THEN
            OUT=.TRUE.
            IWEND(NWORD)=I-1
          ELSEIF(I.EQ.LENGTH)THEN
            IWEND(NWORD)=I
          ENDIF
        ENDIF
   10 CONTINUE
      RETURN
      END
CDOC END_INTEGER_FUNCTION NWORD
