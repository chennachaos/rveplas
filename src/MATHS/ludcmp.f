C  @(#)   Module:<ludcmp.f>   Version:1.6   Date:10/02/95
      SUBROUTINE LUDCMP
     1(   A          ,INDX       ,
     2    D          ,N          ,NP         ,ERROR      )
C$DP,1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER
     1    NAME*6
      PARAMETER
     1(   NMAX=50    )
      LOGICAL
     1    ERROR
      DIMENSION
     1    A(NP,NP)   ,INDX(N)    ,VV(NMAX)
C$DP,3
      DATA
     1    R0         ,R1         ,TINY       /
     2    0.0D0      ,1.0D0      ,1.0D-19    /
C$SP,3
C     DATA
C    1    R0         ,R1         ,TINY       /
C    2    0.0        ,1.0        ,1.0E-19    /
      DATA NAME/'LUDCMP'/
C***********************************************************************
C Routine to do LU decomposition
C See NUMERICAL RECIPES p35.
C*ACRONYM
C LU_DeCOmPosition
C*DESCRIPTION
C*HISTORY
C Name          Date         Comment
C G.C.Huang    Oct,92      initial coding
C*EXTERNAL
C Arrays
C=A      - LU decomposed matrix
C=INDX   - Permutation vector
C Variables
C=D      - Row interchange indicator
C N      - Size of the problem
C NP     - Physical size of A matrix
C ERROR  - Error flag
C (c) Copyright 1992, Rockfield Software Limited, Swansea, UK
C***********************************************************************
ccccccccccccccD     CALL SENTRY(NAME,MODEDB)
      ERROR=.FALSE.
      D=R1
      DO 12 I=1,N
        AAMAX=R0
        DO 11 J=1,N
          IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
   11   CONTINUE
        IF(AAMAX.EQ.R0)THEN
C Error. Singular matrix encountered.
cccccccccccccc          CALL WRTER('A0176E',NAME,0 )
C Flag the error
          ERROR=.TRUE.
          GOTO 999
        ENDIF
        VV(I)=R1/AAMAX
   12 CONTINUE
      DO 19 J=1,N
        IF (J.GT.1) THEN
          DO 14 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO 13 K=1,I-1
                IF(ABS(A(I,K)).LT.TINY.AND.ABS(A(K,J)).LT.TINY)GOTO 13
C If A(I,K) and A(K,J) are very small, skip the following line. Otherwise
C we have problem here. (G.C.)
                SUM=SUM-A(I,K)*A(K,J)
   13         CONTINUE
              A(I,J)=SUM
            ENDIF
   14     CONTINUE
        ENDIF
        AAMAX=R0
        DO 16 I=J,N
          SUM=A(I,J)
          IF (J.GT.1)THEN
            DO 15 K=1,J-1
              IF(ABS(A(I,K)).LT.TINY.AND.ABS(A(K,J)).LT.TINY)GOTO 15
C If A(I,K) and A(K,J) are very small, skip the following line. Otherwise
C we have problem here. (G.C.)
              SUM=SUM-A(I,K)*A(K,J)
   15       CONTINUE
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*ABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
   16   CONTINUE
        IF (J.NE.IMAX)THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
   17     CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(J.NE.N)THEN
          IF(A(J,J).EQ.TINY)A(J,J)=TINY
          DUM=R1/A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
   18     CONTINUE
        ENDIF
   19 CONTINUE
      IF(A(N,N).EQ.R0)A(N,N)=TINY
 999  CONTINUE
ccccccccccccccD     CALL SEXIT(MODEDB)
      RETURN
      END
