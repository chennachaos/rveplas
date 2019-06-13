CDOC BEGIN_DOUBLE_PRECISION_FUNCTION EXP2X
CDOC f(x)=exp(2x)
CDOC
CDOC This function relates the eigenvalues of the Cauchy-Green strain
CDOC tensors to the principal logarithmic stretches.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION X      >  Point at which the function will be
CDOC C                          evaluated.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, August 1996: Initial coding
CHST
      DOUBLE PRECISION FUNCTION EXP2X(X)
C
      DOUBLE PRECISION X
      DOUBLE PRECISION R2
      DATA   R2  /2.0D0/
C***********************************************************************
C SCALAR FUNCTION THAT RELATES EIGENVALUES OF THE CAUCHY-GREEN
C TENSOR TO THE PRINCIPAL LOGARITHMIC STRECTHES
C
C REFERENCE: Section 3.1.7
C***********************************************************************
      EXP2X=DEXP(R2*X)
C
      RETURN
      END
CDOC END_DOUBLE_PRECISION_FUNCTION EXP2X
