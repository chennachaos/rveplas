CDOC BEGIN_INTEGER_FUNCTION NFUNC
CDOC Simple integer calculation used in frontal slover
CDOC
      INTEGER FUNCTION NFUNC(N1,N2)
      I = N1
      J = N2
      NF = (J*J-J)/2+I
      NFUNC=NF
      RETURN
      END
CDOC END_INTEGER_FUNCTION NFUNC
