C  @(#)   Module:<lubksb.f>   Version:1.4   Date:05/03/94
      SUBROUTINE LUBKSB
     1(   A          ,B          ,INDX       ,
     2    N          ,NP         )
C$DP,1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER
     1    NAME*6
      DIMENSION
     1    A(NP,NP)   ,INDX(N)    ,B(N)
      DATA NAME/'LUBKSB'/
C***********************************************************************
C Routine to solve the set of N linear equations AX=B
C See NUMERICAL RECIPES p36.
C*ACRONYM
C LU_BacK_SuBustitutions
C*DESCRIPTION
C*HISTORY
C Name          Date         Comment
C G.C.Huang    Oct,92      initial coding
C*EXTERNAL
C Arrays
C A      - LU decomposed matrix
C=B      - Right hand side matrix as input and stored solutions as output
C INDX   - Permutation vector
C Variables
C N      - Size of the problem
C NP     - Physical size of A matrix
C (c) Copyright 1992, Rockfield Software Limited, Swansea, UK
C***********************************************************************
cccccccccccccccD     CALL SENTRY(NAME,MODEDB)
      II=0
      DO 12 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO 11 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
   11     CONTINUE
        ELSE IF (SUM.NE.0.) THEN
          II=I
        ENDIF
        B(I)=SUM
   12 CONTINUE
      DO 14 I=N,1,-1
        SUM=B(I)
        IF(I.LT.N)THEN
          DO 13 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
   13     CONTINUE
        ENDIF
        B(I)=SUM/A(I,I)
   14 CONTINUE
cccccccccccccccD     CALL SEXIT(MODEDB)
      RETURN
      END
