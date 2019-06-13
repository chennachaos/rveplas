CDOC BEGIN_DOUBLE_PRECISION_FUNCTION DLGD2
CDOC dlg2(x)=log(x)/2
CDOC
CDOC This function relates the principal logarithmic stretches and the
CDOC eigenvalues of the Cauchy-Green strain tensor.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION X      >  Point at which the function will be
CDOC C                          evaluated.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, August 1996; Initial coding
CHST
      DOUBLE PRECISION FUNCTION DLGD2(X)
C
      DOUBLE PRECISION X
      DOUBLE PRECISION RP5
      DATA   RP5 /0.5D0/
C***********************************************************************
C SCALAR FUNCTION THAT RELATES PRINCIPAL LOGARITHMIC STRECTHES AND
C EIGENVALUES OF THE CAUCHY-GREEN TENSOR
C
C REFERENCE: Section 3.1.7
C***********************************************************************
      DLGD2=RP5*DLOG(X)
C
      RETURN
      END
CDOC END_DOUBLE_PRECISION_FUNCTION DLGD2
