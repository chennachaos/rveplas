CDOC BEGIN_SUBROUTINE RDDAMA
CDOC Reads data for Lemaitre's ductile damage elasto-plastic model.
CDOC
CDOC This routine reads from the data file and echos to the results file
CDOC the material parameters necessary for Lemaitre's ductile damage
CDOC elasto-plastic model with piece-wise linear isotropic hardening.
CDOC It also sets the array of real properties and some components of
CDOC the array of integer material properties.
CDOC These arrays are used by subroutines SUDAMA and CTDAMA.
CDOC This routine also checks whether the number of integer and real
CDOC material properties, the number of real state variables and the
CDOC number of logical algorithmic variables required by the present
CDOC model is compatible with the dimension of the corresponding global
CDOC arrays of HYPLAS.
CDOC It also sets the unsymmetric tangent stiffness flag.
CDOC
CDOC BEGIN_PARAMETERS
CDOC INTEGER          IPROPS <  Array of integer material properties.
CDOC INTEGER          MIPROP >  Dimension of the global array of integer
CDOC C                          material properties.
CDOC INTEGER          MLALGV >  Dimension of the global array of logical
CDOC C                          algorithmic variables.
CDOC INTEGER          MRPROP >  Dimension of the global array of real
CDOC C                          material variables.
CDOC INTEGER          MRSTAV >  Dimension of the global array of real
CDOC C                          state variables.
CDOC INTEGER          NTYPE  >  Stress state type.
CDOC DOUBLE_PRECISION RPROPS <  Array of real material properties.
CDOC LOGICAL          UNSYM  <  Unsymmetric tangent stiffness flag.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, January 2000: Initial coding
CHST
CHST E.de Souza Neto, April 2011: I/O error message added
CHST
      SUBROUTINE RDDAMA
     1(   IPROPS     ,MIPROP     ,MLALGV     ,MRPROP     ,MRSTAV     ,
     2    NTYPE      ,RPROPS     ,UNSYM      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL UNSYM
      PARAMETER( IPHARD=6  ,NLALGV=2  ,NRSTAV=6 )
      DIMENSION
     1    IPROPS(*)          ,RPROPS(*)
      DATA R0   /0.0D0/
C***********************************************************************
C READS AND ECHOES MATERIAL PROPERTIES FOR LEMAITRE'S DUCTILE DAMAGE
C ELASTO-PLASTIC MATERIAL MODEL WITH NON-LINEAR (PIECEWISE LINEAR)
C ISOTROPIC HARDENING
C
C REFERENCE Box: 12.3
C***********************************************************************
 1000 FORMAT(' LEMAITRE''S DUCTILE DAMAGE elasto-plastic model'/)
 1100 FORMAT(
     1' Mass density ...................................... =',G15.6/
     2' Young''s modulus ................................... =',G15.6/
     3' Poisson''s ratio ................................... =',G15.6/
     4' Damage evolution law exponent ...................... =',G15.6/
     5' Damage evolution law denominator ................... =',G15.6)
 1200 FORMAT(/
     1' Number of points on hardening curve ............... =',I3//
     2'             R          uniaxial yield stress '/)
 1300 FORMAT(2(5X,G15.6))
C
      IF(NTYPE.NE.2.AND.NTYPE.NE.3)CALL ERRPRT('ED0172')
C Set unsymmetric tangent stiffness flag
      UNSYM=.TRUE.
C
C Read and echo some of the real properties
      WRITE(16,1000)
      READ(15,*,ERR=100,END=100)DENSE
      READ(15,*,ERR=100,END=100)YOUNG,POISS,DAMEXP,DAMDEN
      WRITE(16,1100)DENSE,YOUNG,POISS,DAMEXP,DAMDEN
      IF(YOUNG.LE.R0)CALL ERRPRT('ED0164')
      IF(DAMEXP.LE.R0)CALL ERRPRT('ED0170')
      IF(DAMDEN.LE.R0)CALL ERRPRT('ED0171')
C number of points on hardening curve
      READ(15,*,ERR=100,END=100)NHARD
      WRITE(16,1200)NHARD
      IF(NHARD.LT.2) CALL ERRPRT('ED0165')
C check dimensions of IPROPS
      IF(MIPROP.LT.3)CALL ERRPRT('ED0166')
      IPROPS(3)=NHARD
C check dimensions of RPROPS
      NRPROP=IPHARD+NHARD*2-1
      IF(NRPROP.GT.MRPROP)CALL ERRPRT('ED0167')
C
      RPROPS(1)=DENSE
      RPROPS(2)=YOUNG
      RPROPS(3)=POISS
      RPROPS(4)=DAMEXP
      RPROPS(5)=DAMDEN
C Read and set hardening curve
      DO 10 IHARD=1,NHARD
        READ(15,*,ERR=100,END=100)RPROPS(IPHARD+IHARD*2-2),
     1                            RPROPS(IPHARD+IHARD*2-1)
        WRITE(16,1300)RPROPS(IPHARD+IHARD*2-2),
     1                RPROPS(IPHARD+IHARD*2-1)
   10 CONTINUE
C Check dimension of RSTAVA and LALGVA
      IF(NRSTAV.GT.MRSTAV)CALL ERRPRT('ED0168')
      IF(NLALGV.GT.MLALGV)CALL ERRPRT('ED0169')
C
      GOTO 200
C Issue error message and abort program execution in case of I/O error
  100 CALL ERRPRT('ED0201')
  200 CONTINUE
C
      RETURN
      END
CDOC END_SUBROUTINE RDDAMA
