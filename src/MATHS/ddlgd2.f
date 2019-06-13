CDOC BEGIN_DOUBLE_PRECISION_FUNCTION DDLGD2
CDOC ddlg2(x)=1/(2x)
CDOC
CDOC This is the derivative of the function defined in DLGD2,
CDOC that relates the principal logarithmic stretches and the
CDOC eigenvalues of the Cauchy-Green strain tensor.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION X      >  Point at which the function will be
CDOC C                          evaluated.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, August 1996: Initial coding
CHST
      DOUBLE PRECISION FUNCTION DDLGD2(X)
C
      DOUBLE PRECISION X
      DOUBLE PRECISION RP5
      DATA   RP5 /0.5D0/
C***********************************************************************
C DERIVATIVE OF THE SCALAR FUNCTION 'DLGD2' THAT RELATES PRINCIPAL
C LOGARITHMIC STRECTHES AND EIGENVALUES OF THE CAUCHY-GREEN TENSOR
C***********************************************************************
      DDLGD2=RP5/X
C
      RETURN
      END
CDOC END_DOUBLE_PRECISION_FUNCTION DDLGD2
