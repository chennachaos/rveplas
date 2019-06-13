C  @(#)   Module:<rminve.f>   Version:1.6   Date:05/03/94
      SUBROUTINE RMINVE
     1(   A     ,AI     ,NSIZE      ,ERROR )
C$DP,1
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (MSIZE=50)
      CHARACTER*6 NAME
      LOGICAL
     1    ERROR
      DIMENSION
     1    A(NSIZE,NSIZE)     ,AI(NSIZE,NSIZE)
      DIMENSION
     2    INDX(MSIZE)
C$DP,1
      DATA R0,R1/0.0D0,1.0D0/
C$SP,1
C     DATA R0,R1/0.0  ,1.0  /
      DATA NAME/'RMINVE'/
C***********************************************************************
C Evaluate inverse matrix without determinant
C*Coded by
C G.C.Huang, Oct.,1991
C*Arrays
C A     - A matrix
C=AI    - inversed A matrix
C*Variables
C NSIZE - Size of the matrix
C=ERROR - Error flag
C***********************************************************************
cccccccccccccccD     CALL SENTRY(NAME,MODEDB)
C System checks
ccccccccccccccc      IF(NSIZE.LT.1.OR.NSIZE.GT.MSIZE)
ccccccccccccccc     1    CALL SYSCHK(NAME  ,10    ,NSIZE ,0     )
C Set up identity matrix
C ----------------------
      DO 20 ISIZE=1,NSIZE
        DO 10 JSIZE=1,NSIZE
          AI(ISIZE,JSIZE)=R0
   10   CONTINUE
        AI(ISIZE,ISIZE)=R1
   20 CONTINUE
C LU decompose the matrix
      CALL LUDCMP(A    ,INDX  ,DUMMY,NSIZE ,NSIZE ,ERROR )
      IF(ERROR)THEN
C Error. Singular matrix encountered.
        GOTO 999
      ENDIF
C Find inverse by columns
      DO 30 ISIZE=1,NSIZE
        CALL LUBKSB(A     ,AI(1,ISIZE),INDX  ,NSIZE ,NSIZE )
   30 CONTINUE
  999 CONTINUE
cccccccccccccccD     CALL SEXIT(MODEDB)
      RETURN
      END
