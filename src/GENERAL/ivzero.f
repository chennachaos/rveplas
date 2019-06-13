CDOC BEGIN_SUBROUTINE IVZERO
CDOC Zero an integer array
CDOC
CDOC This routine initialises to zero the N components of the integer
CDOC array argument IV.
CDOC 
CDOC BEGIN_PARAMETERS 
CDOC INTEGER          IV     <  Zeroed integer array.
CDOC INTEGER          N      >  Dimension of IV.
CDOC END_PARAMETERS
CDOC
      SUBROUTINE IVZERO
     1(   IV         ,N          )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IV(N)
C***********************************************************************
C INITIALISES TO ZERO AN INTEGER ARRAY OF DIMENSION N
C***********************************************************************
      DO 10 I=1,N
        IV(I)=0
   10 CONTINUE
      RETURN
      END
CDOC END_SUBROUTINE IVZERO
