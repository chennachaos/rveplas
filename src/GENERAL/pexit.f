CDOC BEGIN_SUBROUTINE PEXIT
CDOC Aborts the execution of HYPLAS
CDOC
CDOC This routine closes the open files and stops the execution of
CDOC HYPLAS, sending a message to the results file and to the standard
CDOC output. It is called in emergency situations when a irrecoverable
CDOC error occurs.
CDOC
      SUBROUTINE PEXIT
C Print message
      WRITE(*,'(///15X,A,///)')'Program HYPLAS aborted.'
      WRITE(16,'(///15X,A,///)')'Program HYPLAS aborted.'
C Close files
      CALL FCLOSE
C and exit program
      STOP ' '
      END
CDOC END_SUBROUTINE PEXIT
