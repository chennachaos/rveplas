CDOC BEGIN_SUBROUTINE RVSUB
CDOC Subtracts two double precision vectors
CDOC
CDOC This function subtracts two double precision vectors passed as
CDOC arguments and stores the result in another vector U=V-W.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION U      <  Double precision vector with the result
CDOC C                          V-W.
CDOC DOUBLE_PRECISION V      >  Double precision vector.
CDOC DOUBLE_PRECISION W      >  Double precision vector.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, September 1996: Initial coding
CHST
      SUBROUTINE RVSUB
     1(   U          ,V          ,W          ,N          )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION
     1    U(N)               ,V(N)               ,W(N)
C***********************************************************************
C SUBTRACTS THE VECTOR 'W' FROM THE VECTOR 'V' AND STORE THE RESULT
C IN 'U'. U ,V AND W ARE DOUBLE PRECISION VECTORS OF DIMENSION N.
C***********************************************************************
      DO 10 I=1,N
        U(I)=V(I)-W(I)
   10 CONTINUE
      RETURN
      END
CDOC END_SUBROUTINE RVSUB
